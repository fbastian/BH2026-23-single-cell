#!/bin/bash
cd ~/liver_cross_species/guilliams2022
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
R=/mnt/masld_vault/cross_species/guilliams2022/06_doublets_ambient/final_dx
mkdir -p $R logs/06e; rm -f logs/06e.fail
echo CS114 CS137 CS138 CS139 CS144 CS48 CS52 CS88 CS89 CS92 CS93 CS96 CS97 | tr ' ' '\n' |
  xargs -P 13 -I{} sh -c "~/miniforge3/envs/sc-preprocessing/bin/python -I code/06e_decontx_final_sample.py \$PWD/tmp_split $R {} > logs/06e/{}.log 2>&1 || echo {} >> logs/06e.fail"
n=$(ls $R/*.dx300.npz | wc -l); echo "dx done $n/13"
[ "$n" -eq 13 ] && [ ! -s logs/06e.fail ] || { echo 1 > logs/06e.exit; exit 1; }
~/miniforge3/envs/rapids-sc/bin/python -I code/06f_clean_annotate_split.py > logs/06f.log 2>&1
echo $? > logs/06e.exit
