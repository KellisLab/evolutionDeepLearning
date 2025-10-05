#!/bin/csh -evf

wget http://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam
bedTools bamtofastq -i T_hg19_1000g.bam -fq densiova.fastq

echo DONE
