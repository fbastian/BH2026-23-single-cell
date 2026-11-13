#!/bin/bash
cd ~/liver_cross_species/samap_maps
PATH=$HOME/miniforge3/envs/mmseqs2/bin:$PATH ~/miniforge3/envs/rapids-sc/bin/python -I code/03_build_mmseqs_maps_liver.py > logs/03_mmseqs.log 2>&1; echo $? > logs/03.exit
