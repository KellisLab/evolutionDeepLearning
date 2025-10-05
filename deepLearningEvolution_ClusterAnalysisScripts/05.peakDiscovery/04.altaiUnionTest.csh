#!/bin/csh -ef

#cat hg38Beds/altai.* > unmerged.union.bed 
#~/go/bin/bedMerge unmerged.union.bed altai.unionSet.bed
#rm unmerged.union.bed

mkdir -p altaiActive
#foreach f (hg38Beds/altai.*)
#	set fPrefix = $f:t:r
#	~/go/bin/intervalOverlap $f altai.unionSet.bed altaiActive/$fPrefix.active.txt
#end



echo DONE
