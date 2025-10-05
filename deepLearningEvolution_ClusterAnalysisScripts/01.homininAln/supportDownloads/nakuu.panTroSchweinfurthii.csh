#!/bin/csh -evf

#Nakuu (panTro Schweinfurthii)

#accession SRX237583
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726412/SRR726412 && fastq-dump SRR726412 && rm SRR726412"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726413/SRR726413 && fastq-dump SRR726413 && rm SRR726413"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726415/SRR726415 && fastq-dump SRR726415 && rm SRR726415"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726416/SRR726416 && fastq-dump SRR726416 && rm SRR726416"

#accession SRX237541
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726408/SRR726408 && fastq-dump SRR726408 && rm SRR726408"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726409/SRR726409 && fastq-dump SRR726409 && rm SRR726409"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726410/SRR726410 && fastq-dump SRR726410 && rm SRR726410"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR726411/SRR726411 && fastq-dump SRR726411 && rm SRR726411"

echo DONE
