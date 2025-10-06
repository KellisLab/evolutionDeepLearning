#!/bin/csh -ef

set bg = ../background.merged.bed
mkdir -p enrichmentResults

~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate ../upregulated.merged.bed subSets/allUp.bed $bg enrichmentResults/mpraUp.enformerUp.txt
~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate ../upregulated.merged.bed subSets/allDown.bed $bg enrichmentResults/mpraDown.enformerUp.txt
~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate ../downregulated.merged.bed subSets/allUp.bed $bg enrichmentResults/mpraUp.enformerDown.txt
~/go/bin/overlapEnrichments -trimToRefGenome -relationship any normalApproximate ../downregulated.merged.bed subSets/allDown.bed $bg enrichmentResults/mpraDown.enformerDown.txt

cat enrichmentResults/*.txt | grep -v 'Method' > enrichment.Summary.txt
echo DONE
