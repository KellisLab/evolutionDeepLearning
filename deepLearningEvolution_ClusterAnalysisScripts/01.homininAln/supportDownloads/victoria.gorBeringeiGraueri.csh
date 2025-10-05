#!/bin/csh -evf

#Accession SRX243533
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748192/SRR748192 && fastq-dump SRR748192 && rm SRR748192 && gzip SRR748192.fastq"

#Accession SRX243532
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748191/SRR748191 && fastq-dump SRR748191 && rm SRR748191 && gzip SRR748191.fastq"

#Accession SRX243531
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748190/SRR748190 && fastq-dump SRR748190 && rm SRR748190 && gzip SRR748190.fastq"

#Accession SRX243530
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748189/SRR748189 && fastq-dump SRR748189 && rm SRR748189 && gzip SRR748189.fastq"

#Accession SRX243529
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748188/SRR748188 && fastq-dump SRR748188 && rm SRR748188 && gzip SRR748188.fastq"

#Accession SRX243528
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748187/SRR748187 && fastq-dump SRR748187 && rm SRR748187 && gzip SRR748187.fastq"

echo DONE
