#!/bin/csh -evf

# to run, add paths or install cluster modules for the following required software
set bwa = bwa
set samtools = samtools
set gatk = gatk # this was used in some debugging, but didn't make it to the submitted manuscript
set bedtools = bedtools
#
set branchLength = 0.01
set propIndel = 0 #for now, we're just looking at substitution
set chromSize = 1000000
set chromNum = 2
set transitionBias = 3
set simulateSamEpsilon = 0.01
set simulateSamDeaminationRate = 0.05
set lambda = 0.01
#
mkdir -p multiFa
mkdir -p empiricalMultiFa
mkdir -p conventionalMultiFa
mkdir -p lambdaMultiFa
mkdir -p output

# first, we set up the reference genome
~/go/bin/randSeq -numSeq $chromNum -lenSeq $chromSize -setSeed $seed ref.fa
cp ~/data/reference/hg38.fa ref.fa
set seed = `expr $seed + 1`
$bwa index ref.fa
$samtools faidx ref.fa
touch ref.dict
rm ref.dict
#$gatk CreateSequenceDictionary -R ref.fa # no longer needed, we dropped the gatk analysis in the paper

set seeds = (25 50 75 100 125)
set coverages = (1 5 10 20 30)
set readLengths = (40 75 150)

foreach s ($seeds)
    foreach cov ($coverages)
        foreach readLength ($readLengths)
        ./support.validateWithPrior.csh $s $cov $readLength
        end
    end
end


echo DONE
