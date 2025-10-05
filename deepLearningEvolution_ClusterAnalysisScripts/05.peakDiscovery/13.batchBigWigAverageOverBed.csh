#!/bin/csh -ef

set chainDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/chainOutput

foreach chainFile ($chainDir/*A.merged.chain)
	set sample = $chainFile:t:r:r:r
	echo $sample
	sbatch --wrap="./support.13.bigWigAverageOverBed.csh $sample"
end

echo DONE
