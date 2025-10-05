#!/bin/csh -ef

mkdir -p buildOutput
set ref = /home/rimangan/data/reference/hg38.fa

foreach b (*.bam)
	set bPrefix = $b:r
	mkdir -p buildOutput/$bPrefix.multiFa
	sbatch --wrap="~/go/bin/samAssembler build -qNameA $bPrefix.A -qNameB $bPrefix.B -empiricalPrior priorOutput/$bPrefix.prior.txt -multiFaDir buildOutput/$bPrefix.multiFa $b $ref buildOutput/$bPrefix.A.fa buildOutput/$bPrefix.B.fa"
end

echo DONE
