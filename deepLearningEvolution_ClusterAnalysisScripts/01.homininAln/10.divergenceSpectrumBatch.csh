#!/bin/csh -ef

mkdir -p divergenceSpectra

foreach v (vcf/*.vcf.gz)
	set vPrefix = $v:t:r:r
	sbatch --wrap="python divergenceSpectrum.py $v divergenceSpectra/$vPrefix.txt"
end

echo DONE
