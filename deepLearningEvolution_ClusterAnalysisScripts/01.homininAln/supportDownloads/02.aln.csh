#!/bin/csh -ef

set ref = ../../reference/hg38.fa

foreach fq (*.fastq.gz)
	echo $fq
	set fqPrefix = $fq:r:r
	sbatch --mem=64G -c 8 --wrap="bwa mem -t 8 $ref $fq | samtools sort -o $fqPrefix.hg38aln.sorted.bam"
end

echo DONE
