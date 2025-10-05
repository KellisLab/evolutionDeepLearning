#!/bin/csh -evf

mkdir -p finalReproduciblePeaks/byBiosample

foreach biosampleFile (finalReproduciblePeaks/archaic/*.bed)
	set biosample = $biosampleFile:t:r:r
	cat finalReproduciblePeaks/archaic/$biosample.final.bed finalReproduciblePeaks/modern/$biosample.final.bed finalReproduciblePeaks/greatApe/$biosample.final.bed > finalReproduciblePeaks/$biosample.groupMerged.bed
	~/go/bin/bedMerge finalReproduciblePeaks/$biosample.groupMerged.bed tmp.bed
	~/go/bin/bedFilter -maxLength 2000 tmp.bed tmp2.bed
	~/go/bin/bedFormat -coordName tmp2.bed finalReproduciblePeaks/byBiosample/$biosample.bed
	rm finalReproduciblePeaks/$biosample.groupMerged.bed
end

echo DONE
