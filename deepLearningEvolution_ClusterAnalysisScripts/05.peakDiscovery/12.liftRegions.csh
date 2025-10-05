#!/bin/csh -ef

set chainDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/chainOutput
set liftOver = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/liftOver

mkdir -p finalReproduciblePeaks/assemblyCoords

foreach chainFile ($chainDir/*.chain)
	sbatch --wrap="./support.12.liftRegions.csh $chainFile"
end

echo DONE
