#!/bin/csh -ef


mkdir -p girskis

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180714/suppl/GSE180714%5FMPRA%5Fdata.txt.gz
mv GSE180714_MPRA_data.txt.gz girskis


echo DONE
