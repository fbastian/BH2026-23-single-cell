#!/bin/bash
LOG=$1; SENT=$2
~/miniforge3/envs/masld_integration/bin/python -I ~/liver_cross_species/guilliams2022/code/01_load_raw.py > "$LOG" 2>&1
echo $? > "$SENT"
