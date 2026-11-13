#!/bin/bash
# 10 strict QC -> 11 TotalVI on strict cells -> 11p GPU post + panel. Screen g22_strict; sentinels logs/10.exit 11.exit 11p.exit
P=~/liver_cross_species/guilliams2022; L=$P/logs; E=~/miniforge3/envs
cd $P; rm -f $L/10.exit $L/11.exit $L/11p.exit
[ "$(df --output=avail / | tail -1)" -gt 150000000 ] || { echo 1 > $L/10.exit; exit 1; }
CUDA_VISIBLE_DEVICES=0 $E/rapids-sc/bin/python -I code/10_strict_qc.py > $L/10_strict.log 2>&1; echo $? > $L/10.exit
[ "$(cat $L/10.exit)" = 0 ] || exit 1
CUDA_VISIBLE_DEVICES=0 $E/masld_integration/bin/python -I code/11_totalvi_strict.py > $L/11_totalvi_strict.log 2>&1; echo $? > $L/11.exit
[ "$(cat $L/11.exit)" = 0 ] || exit 1
CUDA_VISIBLE_DEVICES=0 $E/rapids-sc/bin/python -I code/11p_post_gpu.py > $L/11p.log 2>&1; echo $? > $L/11p.exit
