#!/bin/csh -evf

#Clara (panTro Troglodytes)

#accession SRX243496
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748073/SRR748073 && fastq-dump SRR748073 && rm SRR748073"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748074/SRR748074 && fastq-dump SRR748074 && rm SRR748074"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748075/SRR748075 && fastq-dump SRR748075 && rm SRR748075"

#accession SRX243495
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748071/SRR748071 && fastq-dump SRR748071 && rm SRR748071"
sbatch --wrap="wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR748072/SRR748072 && fastq-dump SRR748072 && rm SRR748072"

echo DONE
