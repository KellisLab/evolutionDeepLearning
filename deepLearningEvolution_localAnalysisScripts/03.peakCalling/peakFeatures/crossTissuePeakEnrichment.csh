#!/bin/csh -ef

touch enrichmentResults.txt
rm enrichmentResults.txt

set noGap = /Users/dickinsonia/Scavenger/Datasets/Genomes/hg38/hg38.noGap.bed

foreach aFile (byBiosample/*.bed)
	set aPrefix = $aFile:t:r
	foreach bFile(byBiosample/*.bed)
		set bPrefix = $bFile:t:r
		echo $aPrefix $bPrefix
		set result = `~/go/bin/overlapEnrichments -trimToRefGenome normalApproximate $aFile $bFile $noGap stdout | grep -v Method | awk '{print $9}'`
		echo $aPrefix $bPrefix $result >> enrichmentResults.txt
	end
end


echo DONE
