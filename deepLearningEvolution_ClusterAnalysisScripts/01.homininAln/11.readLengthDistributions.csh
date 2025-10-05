#!/bin/csh -ef

mkdir -p readLengthDistributions
set humanBamDir = /net/bmc-lab4/data/kellis/users/tnikitha/humanAssemblies/finishedHumanBamFiles

foreach b (finishedBamFiles/*.bam)
	echo $b
	set bPrefix = $b:t:r
	sbatch --wrap="~/go/bin/samInfo readLength $b readLengthDistributions/$bPrefix.readLengthDistribution.txt"
end

foreach b ($humanBamDir/*.sorted.bam)
	echo $b
	set bPrefix = $b:t:r:r
	sbatch --wrap="~/go/bin/samInfo readLength $b readLengthDistributions/$bPrefix.readLengthDistribution.txt"
end

echo DONE
