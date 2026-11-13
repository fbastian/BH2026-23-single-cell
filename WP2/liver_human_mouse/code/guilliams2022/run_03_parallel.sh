#!/bin/bash
P=~/liver_cross_species/guilliams2022; L=$P/logs; E=~/miniforge3/envs; O=/mnt/masld_vault/cross_species/guilliams2022
until grep -q "^written" $L/02_qc.log; do grep -q "^Traceback" $L/02_qc.log && { echo "02 failed" > $L/03.exit; exit 1; }; sleep 15; done
echo 0 > $L/02.exit
CUDA_VISIBLE_DEVICES=1 $E/rapids-sc/bin/python -I $P/code/03h_harmony_gpu.py $O > $L/03_harmony_gpu.log 2>&1 && echo 0 > $L/03h.exit || echo 1 > $L/03h.exit &
CUDA_VISIBLE_DEVICES=0 $E/masld_integration/bin/python -I $P/code/03_integrate.py $O totalvi > $L/03_totalvi.log 2>&1; echo $? > $L/03t.exit
wait
