#!/bin/csh -ef

mkdir -p priorOutput
set ref = ~/data/reference/hg38.fa

foreach b (*.bam)
	set bPrefix = $b:r
	sbatch --mem=32G --wrap="~/go/bin/samAssembler prior $b $ref priorOutput/$bPrefix.prior.txt"
end
