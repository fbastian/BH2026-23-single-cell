#!/bin/bash
P=~/liver_cross_species/guilliams2022
~/miniforge3/envs/rapids-sc/bin/python -I $P/code/02b_preqc_L0_metrics.py > $P/logs/02b.log 2>&1; echo $? > $P/logs/02b.exit
