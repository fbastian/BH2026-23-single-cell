"""Execute the light notebooks 01-08 in order (kernel masld_atlas); executed copies -> notebooks/executed/. Stops at the first failure."""
import os, sys, time, nbformat
from nbclient import NotebookClient
NB = os.path.expanduser("~/liver_cross_species/notebooks"); os.makedirs(f"{NB}/executed", exist_ok=True)
names = sorted(f for f in os.listdir(NB) if f.endswith(".ipynb") and f[:2] in {"01", "02", "03", "04", "05", "06", "07", "08"})
# usage: execute_light.py [first] or execute_light.py 02 03 ... 01 (explicit order)
order = sys.argv[1:] if len(sys.argv) > 2 else None
start = sys.argv[1] if len(sys.argv) == 2 else "01"
todo = [next(x for x in names if x[:2] == k) for k in order] if order else [x for x in names if x[:2] >= start]
for n in todo:
    t0 = time.time(); nb = nbformat.read(f"{NB}/{n}", as_version=4)
    try:
        NotebookClient(nb, timeout=None, kernel_name="masld_atlas", resources={"metadata": {"path": f"{NB}/work"}}).execute()
        ok = True
    except Exception as e:
        ok = False; print(n, "FAILED:", repr(e)[:2000], flush=True)
    nbformat.write(nb, f"{NB}/executed/{n}")
    print(n, "ok" if ok else "FAILED", f"{(time.time() - t0) / 60:.1f} min", flush=True)
    if not ok: sys.exit(1)
print("all light notebooks executed")
