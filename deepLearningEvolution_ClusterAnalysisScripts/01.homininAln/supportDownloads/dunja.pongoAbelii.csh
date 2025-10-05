#!/bin/csh -evf

#Dunja (pongo Abelii)

#accession SRX243482
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748038/SRR748038 && fastq-dump SRR748038 && rm SRR748038 && gzip SRR748038.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748039/SRR748039 && fastq-dump SRR748039 && rm SRR748039 && gzip SRR748039.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748040/SRR748040 && fastq-dump SRR748040 && rm SRR748040 && gzip SRR748040.fastq"

#accession SRX243481
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748035/SRR748035 && fastq-dump SRR748035 && rm SRR748035 && gzip SRR748035.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748036/SRR748036 && fastq-dump SRR748036 && rm SRR748036 && gzip SRR748036.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748037/SRR748037 && fastq-dump SRR748037 && rm SRR748037 && gzip SRR748037.fastq"

echo DONE
