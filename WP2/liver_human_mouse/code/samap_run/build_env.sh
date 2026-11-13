#!/bin/bash
# SAMap env following the Project 23 WP1/SETUP.md recipe (python 3.12, hnswlib from conda-forge, sc-samap 3.0.1). umask 002 (shared conda).
umask 002
source ~/miniforge3/etc/profile.d/conda.sh
conda create -y -n samap python=3.12 hnswlib -c conda-forge > ~/liver_cross_species/samap_run/logs/env.log 2>&1 || { echo 1 > ~/liver_cross_species/samap_run/logs/env.exit; exit 1; }
conda activate samap
pip install anndata scanpy sc-samap==3.0.1 igraph leidenalg pandas numpy scipy holoviews bokeh pyarrow matplotlib >> ~/liver_cross_species/samap_run/logs/env.log 2>&1
python -c "import samap, hnswlib, scanpy; from samap.mapping import SAMAP; from importlib.metadata import version; print('ok', version('sc-samap'), scanpy.__version__)" >> ~/liver_cross_species/samap_run/logs/env.log 2>&1
echo $? > ~/liver_cross_species/samap_run/logs/env.exit
