#!/bin/csh -ef

set hops = /home/rimangan/data/reference/hops/hopsPn_pValue.hg38.bed

mkdir -p hopsOverlapsPValue
foreach b (enhancerPromoterSplits/*.bed)
	set bPrefix = $b:t:r
	~/go/bin/intervalOverlap $b $hops hopsOverlapsPValue/$bPrefix.hopsPValue.txt
end

echo DONE
