#!/bin/csh -ef

#foreach stateDir (epimapEnrichments/*)
#	set state = $stateDir:t
#	echo $state
#	cat $stateDir/*.txt | grep -v 'Method' > epimapEnrichments/$state.txt
#end

set regElements = /home/rimangan/data/reference/epimap/chromatinStates/regulatoryElements.bed
set noGap = /home/rimangan/data/reference/hg38.noGap.bed

~/go/bin/overlapEnrichments -trimToRefGenome -secondFileList peakSets.txt normalApproximate $regElements peakSets.txt $noGap regElementsEnrichments.txt

echo DONE
