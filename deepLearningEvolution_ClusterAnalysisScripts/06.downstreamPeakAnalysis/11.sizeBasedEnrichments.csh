#!/bin/csh -evf

set chromHmmDir = /home/rimangan/data/reference/epimap/chromatinStates/groupByState
set noGap = /home/rimangan/data/reference/hg38.noGap.bed
mkdir -p sizeBasedEnrichments
~/go/bin/bedFilter -maxLength 1000 background.merged.bed sizeBasedEnrichments/background.leq1000.bed
~/go/bin/bedFilter -minLength 1001 background.merged.bed sizeBasedEnrichments/background.greaterThan1000.bed

mkdir -p sizeBasedEnrichments/enrichmentResults
foreach b ($chromHmmDir/*.bed)
	echo $b
	set bPrefix = $b:t:r
	~/go/bin/overlapEnrichments -trimToRefGenome normalApproximate $b sizeBasedEnrichments/background.leq1000.bed $noGap sizeBasedEnrichments/enrichmentResults/$bPrefix.leq1000.bed
	~/go/bin/overlapEnrichments -trimToRefGenome normalApproximate $b sizeBasedEnrichments/background.greaterThan1000.bed $noGap sizeBasedEnrichments/enrichmentResults/$bPrefix.greaterThan1000.bed
end

echo DONE
