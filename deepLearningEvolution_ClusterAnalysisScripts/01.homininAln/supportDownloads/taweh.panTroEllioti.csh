#!/bin/csh -evf

#Taweh (panTro Ellioti)

#accession SRX360479
#wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748156/SRR748156
#wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748157/SRR748157
#wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748158/SRR748158
#wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748159/SRR748159
fastq-dump SRR748156
fastq-dump SRR748157
fastq-dump SRR748158
fastq-dump SRR748159

#accession SRX243518
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748160/SRR748160
fastq-dump SRR748160

cat *.fastq | gzip > taweh.panTroEllioti.fastq.gz
rm *.fastq
rm SRR*

echo DONE
