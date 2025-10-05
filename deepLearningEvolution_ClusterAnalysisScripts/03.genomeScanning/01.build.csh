#!/bin/csh -ef

mkdir -p buildOutput
set ref = ~/data/reference/hg38.simple.fa

set humanBamDir = ~/data/humanAssemblies/finishedHumanBamFiles
set rileyBamDir = /net/bmc-lab4/data/kellis/users/rimangan/homininAln/finishedBamFiles

set buildDir = ~/data/ancoraUpdated/buildOutput
set humanPriorDir = ~/data/humanAssemblies/priorOutput
set rileyPriorDir = /net/bmc-lab4/data/kellis/users/rimangan/homininAln/priorOutput
set coveragePeaks = ~/data/samCoverageTrials/coveragePeaks


foreach b ($humanBamDir/*.sorted.bam)
        set bPrefix = $b:t:r:r
        mkdir -p buildOutput/$bPrefix.multiFa
	sbatch -p kellis --wrap="~/go/bin/ancora build -qNameA $bPrefix.A -qNameB $bPrefix.B -empiricalPrior $humanPriorDir/$bPrefix.prior.txt -multiFaDir $buildDir/$bPrefix.multiFa -problematicRegionsBed $coveragePeaks/$bPrefix.peaks.bed.gz $b $ref $buildDir/$bPrefix.A.fa $buildDir/$bPrefix.B.fa"
end

foreach b ($rileyBamDir/*.bam)
	set bPrefix = $b:t:r
	mkdir -p buildOutput/$bPrefix.multiFa
	sbatch -p kellis --wrap="~/go/bin/ancora build -qNameA $bPrefix.A -qNameB $bPrefix.B -empiricalPrior $rileyPriorDir/$bPrefix.prior.txt -multiFaDir $buildDir/$bPrefix.multiFa -problematicRegionsBed $coveragePeaks/$bPrefix.peaks.bed.gz $b $ref $buildDir/$bPrefix.A.fa $buildDir/$bPrefix.B.fa"
end


echo DONE
