#!/bin/bash
# 07 overrides -> 08 TotalVI on clean counts -> 08p GPU post + panel. Screen g22_clean_tvi; sentinels logs/07.exit 08.exit 08p.exit
P=~/liver_cross_species/guilliams2022; L=$P/logs; E=~/miniforge3/envs
cd $P; rm -f $L/07.exit $L/08.exit $L/08p.exit
$E/rapids-sc/bin/python -I code/07_apply_overrides.py > $L/07.log 2>&1; echo $? > $L/07.exit
[ "$(cat $L/07.exit)" = 0 ] || exit 1
CUDA_VISIBLE_DEVICES=0 $E/masld_integration/bin/python -I code/08_totalvi_clean.py > $L/08_totalvi_clean.log 2>&1; echo $? > $L/08.exit
[ "$(cat $L/08.exit)" = 0 ] || exit 1
CUDA_VISIBLE_DEVICES=0 $E/rapids-sc/bin/python -I code/08p_post_gpu.py > $L/08p.log 2>&1; echo $? > $L/08p.exit
