#!/bin/bash
P=~/liver_cross_species/guilliams2022
CUDA_VISIBLE_DEVICES=1 ~/miniforge3/envs/rapids-sc/bin/python -I $P/code/03h_harmony_gpu.py /mnt/masld_vault/cross_species/guilliams2022 > $P/logs/03_harmony_gpu.log 2>&1; echo $? > $P/logs/03h.exit
