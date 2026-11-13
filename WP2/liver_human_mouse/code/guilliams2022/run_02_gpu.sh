#!/bin/bash
P=~/liver_cross_species/guilliams2022
CUDA_VISIBLE_DEVICES=1 ~/miniforge3/envs/rapids-sc/bin/python -I $P/code/02_qc.py > $P/logs/02_qc.log 2>&1
