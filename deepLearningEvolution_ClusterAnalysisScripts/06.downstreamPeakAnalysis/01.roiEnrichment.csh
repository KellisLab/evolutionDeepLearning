#!/bin/csh -ef

mkdir -p enrichmentResults

foreach roiFile (regionsOfInterest/*.bed)
	set roi = $roiFile:t:r
	foreach biosampleBackgroundFile (backgroundRegions/*.bed)
		set biosample = $biosampleBackgroundFile:t:r:r
		echo $biosample
		~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate upregulatedResults/$biosample.upregulated.bed $roiFile $biosampleBackgroundFile enrichmentResults/$roi.$biosample.upregulated.txt
		~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate downregulatedResults/$biosample.downregulated.bed $roiFile $biosampleBackgroundFile enrichmentResults/$roi.$biosample.downregulated.txt
	end
end

echo DONE
