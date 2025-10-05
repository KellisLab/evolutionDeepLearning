#!/bin/csh -ef

set chromHmmDir = ~/data/reference/epimap/chromatinStates/groupByState/
set noGap = ~/data/reference/hg38.noGap.bed # if using whole genome background instead of background.merged.bed

mkdir -p epimap
foreach b (*.merged.bed)
	set bPrefix = $b:t:r:r
	mkdir -p epimap/$bPrefix
	foreach c ($chromHmmDir/*.bed)
		set cPrefix = $c:t:r
		~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate $c $b background.merged.bed epimap/$bPrefix/$cPrefix.bed
	end
end


echo DONE
