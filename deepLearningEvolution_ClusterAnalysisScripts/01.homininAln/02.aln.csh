#!/bin/csh -ef

set ref = ../reference/hg38.fa

foreach fq (*.fastq.gz)
	echo $fq
	set fqPrefix = $fq:r:r
	sbatch --mem=16G -c 8 --wrap="bwa mem -t 8 $ref $fq | samtools view -b > $fqPrefix.hg38aln.bam"
end

echo DONE
