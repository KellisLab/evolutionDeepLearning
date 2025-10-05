#!/bin/csh -ef

set state = $1
set epimapDir = /home/rimangan/data/reference/epimap/chromatinStates/groupByState
set noGap = /home/rimangan/data/reference/hg38.noGap.bed
mkdir -p epimapEnrichments/$state

foreach groupFile ($epimapDir/*.$state.bed)
	set groupPrefix = $groupFile:t:r:r
	~/go/bin/overlapEnrichments -secondFileList peakSets.txt normalApproximate $groupFile peakSets.txt $noGap epimapEnrichments/$state/$groupPrefix.$state.txt
end

echo DONE
