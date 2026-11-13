#!/bin/bash
cd ~/liver_cross_species/guilliams2022
~/miniforge3/envs/rapids-sc/bin/python -I code/05_draft_annotation.py > logs/05_annotation.log 2>&1
echo $? > logs/05.exit
