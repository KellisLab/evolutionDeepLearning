#!/bin/csh -ef

set sample = $1
#set bedDir = ~/data/ancoraUpdated/enformerPredictions
set bedDir = ~/data/ancoraUpdated/enformerPredictionsB
set bigWigfunc = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/wigToBigWig
set sizes = ~/data/ancoraUpdated/chromSizesOutput/$sample.chrom.sizes

mkdir -p ~/data/ancoraUpdated/bigWigOutput/$sample

foreach track ($bedDir/$sample/*.bed)
	set tPrefix = $track:t:r
	~/go/bin/bedToWig -useRange Annotation $track $sizes /dev/stdout | $bigWigfunc /dev/stdin $sizes ~/data/ancoraUpdated/bigWigOutput/$sample/$sample.$tPrefix.bw && rm $track
end

echo DONE

