#!/bin/csh -evf

#Accession SRX243488
sbatch --wrap="wget https://trace.ncbi.nlm.nih.gov/Traces?run=SRR748053 && fastq-dump SRR748053 && rm SRR748053 && gzip SRR748053.fastq"
sbatch --wrap="wget https://trace.ncbi.nlm.nih.gov/Traces?run=SRR748054 && fastq-dump SRR748054 && rm SRR748054 && gzip SRR748054.fastq"
sbatch --wrap="wget https://trace.ncbi.nlm.nih.gov/Traces?run=SRR748055 && fastq-dump SRR748055 && rm SRR748055 && gzip SRR748055.fastq"

#Accession SRX243487
sbatch --wrap="wget https://trace.ncbi.nlm.nih.gov/Traces?run=SRR748051 && fastq-dump SRR748051 && rm SRR748051 && gzip SRR748051.fastq"
sbatch --wrap="wget https://trace.ncbi.nlm.nih.gov/Traces?run=SRR748052 && fastq-dump SRR748052 && rm SRR748052 && gzip SRR748052.fastq"

echo DONE
