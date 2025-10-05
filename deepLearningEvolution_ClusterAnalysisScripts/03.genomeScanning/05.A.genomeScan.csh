#!/bin/csh -evf

set trackList = ~/data/ancoraUpdated/trackList.txt

mkdir -p enformerPredictions

set outDir = ~/data/ancoraUpdated/enformerPredictions
set rileyDir = /net/bmc-lab4/data/kellis/users/rimangan
set model = $rileyDir/enformer/enformer_1

foreach f (buildOutput/*.A.fa)
	set sample = $f:t:r
	mkdir -p $outDir/$sample
	sbatch -p kellis --time-min="20-0:0:0" --wrap="python ~/data/enformer/genomeScan/genomeScan.py -c chromSizesOutput/$sample.chrom.sizes -f $f -t $trackList -o $outDir/$sample -m $model"
end

echo DONE
