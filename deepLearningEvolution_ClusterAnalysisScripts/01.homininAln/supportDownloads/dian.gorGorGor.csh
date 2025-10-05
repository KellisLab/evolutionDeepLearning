#!/bin/csh -evf

#Dian (gorGorGor)
#SRX243460
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747981/SRR747981 && fastq-dump SRR747981 && rm SRR747981 && gzip SRR747981.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747982/SRR747982 && fastq-dump SRR747982 && rm SRR747982 && gzip SRR747982.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747983/SRR747983 && fastq-dump SRR747983 && rm SRR747983 && gzip SRR747983.fastq"

#SRX243459
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747978/SRR747978 && fastq-dump SRR747978 && rm SRR747978 && gzip SRR747978.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747979/SRR747979 && fastq-dump SRR747979 && rm SRR747979 && gzip SRR747979.fastq"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR747980/SRR747980 && fastq-dump SRR747980 && rm SRR747980 && gzip SRR747980.fastq"

#SRX243538
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748197/SRR748197 && fastq-dump SRR748197 && rm SRR748197 && gzip SRR748197.fastq"

#SRX243537
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748196/SRR748196 && fastq-dump SRR748196 && rm SRR748196 && gzip SRR748196.fastq"

echo DONE
