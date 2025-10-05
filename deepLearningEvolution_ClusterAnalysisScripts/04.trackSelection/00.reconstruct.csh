#!/bin/csh -evf

mkdir -p reconOut
set tree = 4d.mod

foreach f (*.fa.gz)
	set fPrefix = $f:r:r
	sbatch --mem=16G --wrap="~/go/bin/reconstructSeq $tree $f reconOut/$fPrefix.fa.gz"
end

echo DONE
