#!/bin/bash
# chain: prep -> SAMap; exit sentinels in logs/. Free-space guard (vault >= 40 GB).
umask 002
L=~/liver_cross_species/samap_run/logs
source ~/miniforge3/etc/profile.d/conda.sh; conda activate samap
free=$(df --output=avail -BG /mnt/masld_vault | tail -1 | tr -dc 0-9); [ "$free" -lt 40 ] && { echo "vault only ${free}G free" > $L/prep.log; echo 9 > $L/prep.exit; exit 9; }
cd ~/liver_cross_species/samap_run
python code/01_prep.py > $L/prep.log 2>&1; echo $? > $L/prep.exit; [ "$(cat $L/prep.exit)" = 0 ] || exit 1
python code/02_samap.py > $L/samap.log 2>&1; echo $? > $L/samap.exit
