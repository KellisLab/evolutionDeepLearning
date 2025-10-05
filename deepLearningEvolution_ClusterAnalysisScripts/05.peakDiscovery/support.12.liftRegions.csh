#!/bin/csh -ef

set liftOver = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/liftOver
set chainFile = $1

set haplotype = $chainFile:t:r:r
mkdir -p finalReproduciblePeaks/assemblyCoords/$haplotype
foreach biosampleFile (finalReproduciblePeaks/byBiosample/*.bed)
	set biosample = $biosampleFile:t:r
	$liftOver $biosampleFile $chainFile finalReproduciblePeaks/assemblyCoords/$haplotype/$biosample.bed /dev/null
end

echo DONE
