#!/bin/csh - ef

#set bedDir = ~/data/ancoraUpdated/enformerPredictions

set bedDir = ~/data/ancoraUpdated/enformerPredictionsB

#mkdir -p ~/data/ancoraUpdated/bigWigOutput

foreach b ($bedDir/*)
	set sample = $b:t
	sbatch -p kellis --wrap="./support.06.csh $sample"
end

echo DONE
