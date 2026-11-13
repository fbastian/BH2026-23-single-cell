#!/bin/bash
P=~/liver_cross_species/guilliams2022; L=$P/logs; E=~/miniforge3/envs; O=/mnt/masld_vault/cross_species/guilliams2022
CUDA_VISIBLE_DEVICES=0 $E/masld_integration/bin/python -I $P/code/03_integrate.py $O totalvi > $L/03_totalvi.log 2>&1; echo $? > $L/03t.exit
[ "$(cat $L/03t.exit)" = 0 ] || exit 1
CUDA_VISIBLE_DEVICES=0 $E/rapids-sc/bin/python -I $P/code/03p_post_gpu.py $O > $L/03p.log 2>&1; echo $? > $L/03p.exit
$E/rapids-sc/bin/python -I $P/code/04_umap_draft_panel.py $O > $L/04.log 2>&1; echo $? > $L/04.exit
