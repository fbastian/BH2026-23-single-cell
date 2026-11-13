"""Shared setup for the reproduction notebooks (Guilliams 2022 mouse atlas rebuild + HLiCA x mouse SAMap, 2026-11-12).

Every step's code is written verbatim to WORK by the notebook (%%writefile) and run here as a subprocess in the conda env it
ran in originally (the kernel only orchestrates). Paths in the code point to the REPRO root, a mirror of the RECORD root.

Safety rules (outputs of record are never written):
  * link(rel)  - an input that this notebook does not rebuild is symlinked from RECORD into REPRO.
  * fresh(rel) - an output folder is made real before a step writes it: any symlinked ancestor inside REPRO is converted to
                 a real folder whose children are linked back to RECORD, and symlinks inside the folder itself are removed.
  * run(), launch()/wait_for() - snapshot RECORD before each step; raise if any record file changed, appeared or disappeared.
Load with %run -i (shared namespace, so RUN_HEAVY set in the notebook is seen by launch/wait_for).
"""
import os, shutil, subprocess, time, glob
import pandas as pd

RECORD = "/mnt/masld_vault/cross_species"
REPRO = "/mnt/masld_vault/cross_species_repro"
NB = os.path.expanduser("~/liver_cross_species/notebooks")
WORK = f"{NB}/work"; LOGS = f"{WORK}/logs"
ENVS = os.path.expanduser("~/miniforge3/envs")
PY = {e: f"{ENVS}/{e}/bin/python" for e in ("rapids-sc", "masld_integration", "sc-preprocessing", "samap", "masld_atlas")}
RSCRIPT = {e: f"{ENVS}/{e}/bin/Rscript" for e in ("masld_tables",)}
RUN_HEAVY = False   # heavy notebooks (09-12) do nothing until this is set to True in the notebook (setup is loaded with %run -i)
_LAUNCHED = {}
os.makedirs(LOGS, exist_ok=True); os.makedirs(REPRO, exist_ok=True)
os.umask(0o002)


def rec(rel): return f"{RECORD}/{rel}"
def rep(rel): return f"{REPRO}/{rel}"


def free_gb(path=REPRO):
    s = os.statvfs(path); return s.f_bavail * s.f_frsize / 1e9


def guard(gb, path=REPRO):
    f = free_gb(path)
    if f < gb: raise RuntimeError(f"only {f:.0f} GB free on {path}, this step needs >= {gb} GB")
    print(f"free space on {path}: {f:.0f} GB (needs {gb})")


def link(*rels):
    """Symlink RECORD inputs into REPRO unless REPRO already has them (rebuilt by an earlier notebook)."""
    for rel in rels:
        dst = rep(rel)
        if os.path.lexists(dst):
            print(("rebuilt   " if not os.path.islink(dst) else "linked    ") + rel); continue
        _materialise(os.path.dirname(dst))
        if os.path.lexists(dst):   # the parent was just created with its record children linked (fix 2026-11-13)
            print("linked    " + rel + "  (record)"); continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        assert os.path.exists(rec(rel)), f"missing in record: {rel}"
        os.symlink(rec(rel), dst); print("linked    " + rel + "  (record)")


def _materialise(path):
    """Make every component of path inside REPRO a real folder; children of converted links point back to RECORD."""
    rel = os.path.relpath(path, REPRO)
    if rel.startswith(".."): raise ValueError(path)
    cur = REPRO
    for part in [p for p in rel.split(os.sep) if p and p != "."]:
        cur = os.path.join(cur, part)
        recp = os.path.join(RECORD, os.path.relpath(cur, REPRO))
        if os.path.islink(cur):
            target = os.readlink(cur); os.unlink(cur); os.makedirs(cur)
            for ch in os.listdir(target): os.symlink(os.path.join(target, ch), os.path.join(cur, ch))
        elif not os.path.exists(cur):
            # review fix 2026-11-12: a folder not yet in REPRO but present in RECORD is created with its record children
            # linked (e.g. 08_totalvi_clean/embeddings.h5ad when only streaks/ is rebuilt)
            os.makedirs(cur)
            if os.path.isdir(recp):
                for ch in os.listdir(recp): os.symlink(os.path.join(recp, ch), os.path.join(cur, ch))


def fresh(*rels):
    """Output folders: real folders, no symlinked entries inside (so writes never follow a link into RECORD)."""
    for rel in rels:
        d = rep(rel); _materialise(d)
        for e in os.scandir(d):
            if e.is_symlink(): os.unlink(e.path)
        print("output    " + rel)


def is_real(rel):
    p = rep(rel); rel_ = os.path.relpath(p, REPRO); cur = REPRO
    for part in rel_.split(os.sep):
        cur = os.path.join(cur, part)
        if os.path.islink(cur) or not os.path.exists(cur): return False
    return True


def _snapshot():
    out = {}
    for root, _, files in os.walk(RECORD):
        for f in files:
            p = os.path.join(root, f)
            try: out[p] = os.stat(p).st_mtime_ns
            except FileNotFoundError: pass
    return out


def _tail(log, n=25):
    lines = [l.rstrip() for l in open(log, errors="replace") if "warn" not in l.lower() and l.strip()]
    print("\n".join(lines[-n:]))


def _record_changes(before):
    after = _snapshot()
    changed = [k for k, v in after.items() if k in before and before[k] != v]
    new = [k for k in after if k not in before]; gone = [k for k in before if k not in after]
    return changed + new + gone


def run(script, env, args=(), gpu=None, threads=None, path_prefix=None, cwd=None, tail=25):
    """Run WORK/<script> in its conda env; log to WORK/logs/<script>.log; fail loudly; never touch RECORD."""
    envv = os.environ.copy(); envv["MPLBACKEND"] = "Agg"   # the kernel's inline backend is not valid in a subprocess
    if gpu is not None: envv["CUDA_VISIBLE_DEVICES"] = str(gpu)
    if threads: envv.update({k: str(threads) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})
    if path_prefix: envv["PATH"] = path_prefix + ":" + envv["PATH"]
    s = f"{WORK}/{script}"
    cmd = ([RSCRIPT[env], s] if script.endswith(".R") else ["bash", s] if script.endswith(".sh") else [PY[env], "-I", s]) + [str(a) for a in args]
    log = f"{LOGS}/{os.path.basename(script)}.log"
    before = _snapshot(); t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd or WORK, env=envv)
    changed = _record_changes(before)
    _tail(log, tail)
    print(f"--- {script} [{env}] exit {p.returncode}, {(time.time() - t0) / 60:.1f} min, log {log}")
    if changed: raise RuntimeError(f"RECORD files changed (must not happen): {changed[:5]}")
    if p.returncode: raise RuntimeError(f"{script} failed, see {log}")


def launch(name, lines):
    """Heavy steps: write a bash wrapper and start it in a detached screen session with an exit sentinel."""
    if not RUN_HEAVY: print("RUN_HEAVY is False: not launched"); return
    sh = f"{WORK}/{name}.sh"; sent = f"{LOGS}/{name}.exit"
    if os.path.exists(sent): os.remove(sent)
    _LAUNCHED[name] = _snapshot()
    with open(sh, "w") as fh:   # body in a subshell, so 'set -e' inside still reaches the sentinel line (review fix)
        fh.write("#!/bin/bash\numask 002\nexport MPLBACKEND=Agg\n(\nset -o pipefail\n" + "\n".join(lines) + f"\n)\necho $? > {sent}\n")
    subprocess.run(["screen", "-dmS", f"nb_{name}", "bash", sh], check=True)
    print(f"launched screen nb_{name}; sentinel {sent}")


def wait_for(name, log=None, every=120):
    if not RUN_HEAVY: print("RUN_HEAVY is False: nothing to wait for"); return
    sent = f"{LOGS}/{name}.exit"
    while not os.path.exists(sent):
        if log and os.path.exists(log):
            l = [x for x in open(log, errors="replace") if x.strip()]
            print(time.strftime("%H:%M"), l[-1].rstrip()[:150] if l else "", flush=True)
        time.sleep(every)
    code = open(sent).read().strip(); print(name, "exit", code)
    if name in _LAUNCHED:
        ch = _record_changes(_LAUNCHED[name])
        if ch: raise RuntimeError(f"RECORD files changed, new or deleted (must not happen): {ch[:5]}")
    if code != "0": raise RuntimeError(f"{name} failed")


def compare_tsv(rel, index_col=0, atol=1e-6):
    """Compare a REPRO table with the RECORD table of the same name."""
    a = pd.read_csv(rec(rel), sep="\t", index_col=index_col); b = pd.read_csv(rep(rel), sep="\t", index_col=index_col)
    same_shape = a.shape == b.shape and list(a.columns) == list(b.columns)
    num = a.select_dtypes("number").columns
    diff = (a[num] - b[num].reindex_like(a[num])).abs().max().max() if same_shape and len(num) else float("nan")
    txt = (a.drop(columns=num).astype(str) == b.drop(columns=num).reindex_like(a.drop(columns=num)).astype(str)).all().all() if same_shape else False
    print(f"{rel}: shape record {a.shape} / repro {b.shape}; max numeric diff {diff:.3g}; text columns identical: {txt}")
    return same_shape and (not len(num) or diff <= atol) and txt


def show(rel, width=900):
    from IPython.display import Image, display
    display(Image(filename=rep(rel), width=width))
