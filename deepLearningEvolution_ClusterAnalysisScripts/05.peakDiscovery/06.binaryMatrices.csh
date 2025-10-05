#!/bin/csh -ef

mkdir -p biosamplePeaksByGenome

foreach unionSet (unionSetsHg38/*.unionSet.bed)
	echo $unionSet
	set biosample = $unionSet:t:r:r
	set outputDir = biosamplePeaksByGenome
	echo $biosample
	mkdir -p biosamplePeaksByGenome/$biosample
	foreach bedFile (hg38Beds/*$biosample*.hg38.bed)
		echo $bedFile
		set genome = $bedFile:t:r:r:r:r:r
		echo $genome
		~/go/bin/intervalOverlap $bedFile $unionSet biosamplePeaksByGenome/$biosample/$genome.bed
	end
end

echo DONE
