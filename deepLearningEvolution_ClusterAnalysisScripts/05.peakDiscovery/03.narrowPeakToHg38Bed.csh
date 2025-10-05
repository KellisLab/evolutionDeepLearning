#!/bin/csh -ef

mkdir -p assemblyCoordsBeds
mkdir -p hg38Beds
set liftOver = /home/rimangan/data/Software/x86_kentBin/liftOver

foreach chain (toHg38Chains/*.B.merged.To.Hg38.chain)
	set chainPrefix = $chain:t:r:r:r:r
	echo $chainPrefix
	foreach narrowPeak (narrowPeak/$chainPrefix.*.narrowPeak)
		echo $narrowPeak
		set narrowPeakPrefix = $narrowPeak:t:r
		cat $narrowPeak | grep -v 'nextItemButton' | cut -f1-4 > assemblyCoordsBeds/$narrowPeakPrefix.assemblyCoords.bed
		cat $narrowPeak | grep -v 'nextItemButton' | cut -f1-4 | $liftOver /dev/stdin $chain hg38Beds/$narrowPeakPrefix.hg38.bed /dev/null
	end
end

echo DONE
