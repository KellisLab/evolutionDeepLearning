#!/bin/csh -ef

touch jaccard_result.txt
rm jaccard_result.txt

foreach aFile (byBiosample/*.bed)
	set aPrefix = $aFile:t:r
	echo $aPrefix
	foreach bFile(byBiosample/*.bed)
		set bPrefix = $bFile:t:r
		set result = `/Users/dickinsonia/Software/bedtools2/bin/bedtools jaccard -a $aFile -b $bFile | grep -v 'intersection' | awk '{print $3}'`
		echo $aPrefix $bPrefix $result >> jaccard_result.txt
	end
end


echo DONE
