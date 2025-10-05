#!/bin/csh -ef

cat finalReproduciblePeaks/archaic/*.bed finalReproduciblePeaks/modern/*.bed finalReproduciblePeaks/greatApe/*.bed > finalReproduciblePeaks/union.unmerged.bed
~/go/bin/bedMerge finalReproduciblePeaks/union.unmerged.bed tmp.bed ; ~/go/bin/bedFormat -coordName tmp.bed unionSet.Idr.bed
rm finalReproduciblePeaks/union.unmerged.bed tmp.bed

mkdir -p finalReproduciblePeaks/byBiosample

foreach biosampleFile (finalReproduciblePeaks/archaic/*.bed)
	set biosample = $biosampleFile:t:r:r
	cat finalReproduciblePeaks/archaic/$biosample.final.bed finalReproduciblePeaks/greatApe/$biosample.final.bed finalReproduciblePeaks/modern/$biosample.final.bed > finalReproduciblePeaks/$biosample.cat.bed
	~/go/bin/intervalOverlap finalReproduciblePeaks/$biosample.cat.bed unionSet.Idr.bed finalReproduciblePeaks/byBiosample/$biosample.bed
	rm finalReproduciblePeaks/$biosample.cat.bed
end

echo DONE
