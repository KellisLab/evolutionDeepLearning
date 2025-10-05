#!/bin/csh -ef

mkdir -p unionSetsHg38

set bwDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/bigWigOutput/altai.Neanderthal.A

foreach bw ($bwDir/*.bw)
	set biosample = $bw:t:r:e
	echo $biosample
	cat hg38Beds/*.$biosample.*.bed > unmerged.unionSet.bed
	~/go/bin/bedMerge unmerged.unionSet.bed unionSetsHg38/$biosample.unionSet.bed
	rm unmerged.unionSet.bed
end

cat unionSetsHg38/*.unionSet.bed > unmerged.unionSet.bed
~/go/bin/bedMerge unmerged.unionSet.bed allBiosamples.unionSet.bed
rm unmerged.unionSet.bed

echo DONE
