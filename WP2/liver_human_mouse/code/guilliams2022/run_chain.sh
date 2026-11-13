#!/bin/bash
# waits for step 01, then 02 QC (rapids-sc env: has scikit-image for scrublet), then 03 Harmony + TotalVI
P=~/liver_cross_species/guilliams2022; L=$P/logs; E=~/miniforge3/envs
[ "$(df --output=avail /mnt/masld_vault | tail -1)" -gt 50000000 ] || { echo "ABORT: <50 GB free on vault" > $L/chain.exit; exit 1; }
while [ ! -f $L/01.exit ]; do sleep 30; done
[ "$(cat $L/01.exit)" = 0 ] || { echo "01 failed" > $L/chain.exit; exit 1; }
$E/rapids-sc/bin/python -I $P/code/02_qc.py > $L/02_qc.log 2>&1; echo $? > $L/02.exit
[ "$(cat $L/02.exit)" = 0 ] || { echo "02 failed" > $L/chain.exit; exit 1; }
$E/masld_integration/bin/python -I $P/code/03_integrate.py /mnt/masld_vault/cross_species/guilliams2022 harmony > $L/03_harmony.log 2>&1; echo $? > $L/03h.exit
$E/masld_integration/bin/python -I $P/code/03_integrate.py /mnt/masld_vault/cross_species/guilliams2022 totalvi > $L/03_totalvi.log 2>&1; echo $? > $L/03t.exit
echo done > $L/chain.exit
