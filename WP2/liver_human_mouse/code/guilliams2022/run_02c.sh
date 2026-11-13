#!/bin/bash
P=~/liver_cross_species/guilliams2022
while [ ! -f $P/logs/02b.exit ]; do sleep 20; done
[ "$(cat $P/logs/02b.exit)" = 0 ] || { echo 99 > $P/logs/02c.exit; exit 1; }
~/miniforge3/envs/rapids-sc/bin/python -I $P/code/02c_thresholds_ridges.py > $P/logs/02c.log 2>&1; echo $? > $P/logs/02c.exit
