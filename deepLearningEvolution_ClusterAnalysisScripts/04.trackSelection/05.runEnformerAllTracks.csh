#!/bin/csh 

foreach split (harSplits/*.hg38.bed)
	sbatch -p kellis --wrap="./support.05.runEnformerAllTracks.csh $split"
end

echo DONE
