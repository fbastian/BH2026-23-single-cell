#!/bin/bash
cd ~/liver_cross_species/guilliams2022
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
echo CS140 CS139 CS93 CS144 CS89 CS115 ABU19 ABU20 ABU21 CS98 | tr ' ' '\n' | xargs -P 10 -I{} sh -c "~/miniforge3/envs/sc-preprocessing/bin/python -I code/06d_decontx_ladder.py $PWD/tmp_split /mnt/masld_vault/cross_species/guilliams2022/06_doublets_ambient/ladder {} > logs/06d/{}.log 2>&1"
echo $? > logs/06d.exit
