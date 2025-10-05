#!/bin/csh -ef

foreach chrom (`cat chromNames.txt`)
	sbatch --wrap="./support.03.csh $chrom"
end

echo DONE
