#!/bin/csh -ef

set sample = $1
set bigWigAverageOverBed = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/bigWigAverageOverBed
set bigWigDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/bigWigOutput
mkdir -p bigWigAverageOverBed
mkdir -p bigWigAverageOverBed/$sample.A
mkdir -p bigWigAverageOverBed/$sample.B

foreach biosampleFile (finalReproduciblePeaks/assemblyCoords/$sample.A/*.bed)	
	set biosample = $biosampleFile:t:r
	$bigWigAverageOverBed $bigWigDir/$sample.A/$sample.A.$biosample.bw $biosampleFile bigWigAverageOverBed/$sample.A/$biosample.txt
	$bigWigAverageOverBed $bigWigDir/$sample.B/$sample.B.$biosample.bw finalReproduciblePeaks/assemblyCoords/$sample.B/$biosample.bed bigWigAverageOverBed/$sample.B/$biosample.txt
end

echo DONE
