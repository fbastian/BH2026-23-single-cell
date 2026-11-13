#!/bin/bash
# cDNA references matching our annotations: human GENCODE v48 (HLiCA, CxG schema 7.1.0), mouse GENCODE vM18 (= Ensembl 93 = 10x mm10-3.0.0)
R=/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/refs/raw; cd $R
for u in https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz \
         https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.transcripts.fa.gz; do
  f=$(basename $u); [ -s $f ] || curl -s --retry 3 -o $f $u; gzip -t $f && echo "$f ok $(stat -c %s $f)" || echo "$f BAD"
done
sha256sum *.fa.gz > SHA256SUMS
