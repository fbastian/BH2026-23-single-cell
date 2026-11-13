#!/bin/bash
# Step 06: split -> scDblFinder + decontX per sample (12 parallel) -> assemble. Screen g22_dbl_dx; sentinel logs/06.exit.
set -uo pipefail
cd ~/liver_cross_species/guilliams2022
[ "$(df --output=avail / | tail -1)" -gt 150000000 ] || { echo "ABORT: <150 GB free on /"; echo 1 > logs/06.exit; exit 1; }
OUT=/mnt/masld_vault/cross_species/guilliams2022
SPLIT=$HOME/liver_cross_species/guilliams2022/tmp_split
RES=$OUT/06_doublets_ambient/per_sample
rm -f logs/06.exit logs/06b.fail; mkdir -p logs/06b
[ "$(cat logs/06a.exit 2>/dev/null)" = 0 ] || { ~/miniforge3/envs/rapids-sc/bin/python -I code/06a_split_samples.py $OUT $SPLIT > logs/06a.log 2>&1; echo $? > logs/06a.exit; [ "$(cat logs/06a.exit)" = 0 ] || { echo 1 > logs/06.exit; exit 1; }; }
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PY=~/miniforge3/envs/sc-preprocessing/bin/python
ls -S $SPLIT/*.h5ad | xargs -n1 basename | sed 's/\.h5ad$//' |
  xargs -P 12 -I{} sh -c "$PY -I code/06b_doublets_decontx_sample.py $SPLIT $RES {} > logs/06b/{}.log 2>&1 || echo {} >> logs/06b.fail"
n=$(ls $RES/*.cells.parquet 2>/dev/null | wc -l); echo "per-sample done: $n/46"
[ "$n" -eq 46 ] && [ ! -s logs/06b.fail ] || { echo 1 > logs/06.exit; exit 1; }
~/miniforge3/envs/rapids-sc/bin/python -I code/06c_assemble.py $OUT > logs/06c.log 2>&1
echo $? > logs/06.exit
