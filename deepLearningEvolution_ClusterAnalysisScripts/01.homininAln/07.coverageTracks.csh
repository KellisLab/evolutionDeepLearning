#!/bin/csh -evf

#OLD: first run module load deeptools and samtools

mkdir -p coverageWigs
set sizes = ~/data/reference/hg38.chrom.sizes
#OLD: set binSize = 10

foreach b (*.bam)
	set bPrefix = $b:r
	sbatch --mem=64G --wrap="~/go/bin/samToWig -defaultValue 0 $b $sizes coverageWigs/$bPrefix.coverage.wig.gz"
end

echo DONE
