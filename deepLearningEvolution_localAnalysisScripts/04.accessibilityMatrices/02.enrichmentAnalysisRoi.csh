#!/bin/csh -ef

set roiFile = ~/Scavenger/Datasets/Genomes/hg38/bed/HGE_merged.bed
set roiPrefix = $roiFile:t:r
mkdir -p enrichmentResults

foreach biosampleFile (divergentPeaks/upregulatedResults/*.bed)
	set biosample = $biosampleFile:t:r:r
	echo $biosample
	~/go/bin/overlapEnrichments -trimToRefGenome normalApproximate $roiFile $biosampleFile divergentPeaks/backgroundRegions/$biosample.background.bed enrichmentResults/$biosample.$roiPrefix.txt
end


echo DONE
