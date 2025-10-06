#!/bin/csh -ef

mkdir -p enrichmentResults

foreach roiFile (regionsOfInterest/*.bed)
	set roi = $roiFile:t:r
	foreach biosampleBackgroundFile (backgroundRegions/*.bed)
		set biosample = $biosampleBackgroundFile:t:r:r
		echo $biosample
		~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInHominins/$biosample.upInHominins.bed $roiFile $biosampleBackgroundFile enrichmentResults/$roi.$biosample.upInHominins.txt
		~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInGreatApes/$biosample.upInGreatApes.bed $roiFile $biosampleBackgroundFile enrichmentResults/$roi.$biosample.upInGreatApes.txt
	end
end

echo DONE
