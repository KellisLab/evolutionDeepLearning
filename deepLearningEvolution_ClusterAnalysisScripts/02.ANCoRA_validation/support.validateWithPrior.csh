#!/bin/csh -evf

set seed = $1
set startingSeed = $seed
set coverage = $2
set readLength = $3
echo $seed
echo $coverage
echo $readLength
set fragmentLength = `expr $readLength \* 2`

# on Luria, add to path and install modules
set bwa = bwa
set samtools = samtools
#set gatk = gatk
set bedtools = bedtools

set branchLength = 0.01
set propIndel = 0 #for now, we're just looking at substitution
set chromSize = 1000000
set chromNum = 2
set transitionBias = 3
set simulateSamEpsilon = 0.01
set simulateSamDeaminationRate = 0.05
set lambda = 0.01

mkdir -p multiFa
mkdir -p empiricalMultiFa
mkdir -p conventionalMultiFa
mkdir -p lambdaMultiFa
mkdir -p output

#first, a reference genome. We'll do a genome with two chromsomes, called Sequence_0 and Sequence_1
#~/go/bin/randSeq -numSeq $chromNum -lenSeq $chromSize -setSeed $seed ref.fa
cp ~/data/reference/hg38.fa ref.fa
set seed = `expr $seed + 1`
$bwa index ref.fa
$samtools faidx ref.fa
touch ref.dict
rm ref.dict
#$gatk CreateSequenceDictionary -R ref.fa

#Next, simulated forward eolution. 
~/go/bin/faFilter -name Sequence_0 ref.fa Sequence_0.fa
~/go/bin/simulateEvol withIndels -qName evol1 -branchLength $branchLength -propIndel $propIndel -setSeed $seed -transitionBias $transitionBias Sequence_0.fa Sequence_0.evol1.fa
set seed = `expr $seed + 1`
~/go/bin/simulateEvol withIndels -qName evol2 -branchLength $branchLength -propIndel $propIndel -setSeed $seed -transitionBias $transitionBias Sequence_0.fa Sequence_0.evol2.fa
~/go/bin/faFilter -name Sequence_1 ref.fa Sequence_1.fa
set seed = `expr $seed + 1`
~/go/bin/simulateEvol withIndels -qName evol1 -branchLength 0.05 -propIndel 0.1 -setSeed $seed -transitionBias 3 Sequence_1.fa Sequence_1.evol1.fa
set seed = `expr $seed + 1`
~/go/bin/simulateEvol withIndels -qName evol2 -branchLength 0.05 -propIndel 0.1 -setSeed $seed -transitionBias 3 Sequence_1.fa Sequence_1.evol2.fa

#I then filtered out the evolved sequences from the resulting multiFa alignments and trimmmed out gaps:

~/go/bin/faFilter -notName Sequence_0 Sequence_0.evol1.fa tmp.Sequence_0.evol1.fa
~/go/bin/faFormat -noGaps tmp.Sequence_0.evol1.fa Sequence_0.evol1.diverge.noGaps.fa
~/go/bin/faFilter -notName Sequence_1 Sequence_1.evol1.fa tmp.Sequence_1.evol1.fa
~/go/bin/faFormat -noGaps tmp.Sequence_1.evol1.fa Sequence_1.evol1.diverge.noGaps.fa

~/go/bin/faFilter -notName Sequence_0 Sequence_0.evol2.fa tmp.Sequence_0.evol2.fa
~/go/bin/faFormat -noGaps tmp.Sequence_0.evol2.fa Sequence_0.evol2.diverge.noGaps.fa
~/go/bin/faFilter -notName Sequence_1 Sequence_1.evol2.fa tmp.Sequence_1.evol2.fa
~/go/bin/faFormat -noGaps tmp.Sequence_1.evol2.fa Sequence_1.evol2.diverge.noGaps.fa
rm tmp*

cat Sequence_0.evol1.diverge.noGaps.fa Sequence_1.evol1.diverge.noGaps.fa > evol1.diverge.noGaps.fa
cat Sequence_0.evol2.diverge.noGaps.fa Sequence_1.evol2.diverge.noGaps.fa > evol2.diverge.noGaps.fa
rm Sequence_0.evol1.diverge.noGaps.fa Sequence_1.evol1.diverge.noGaps.fa Sequence_0.evol2.diverge.noGaps.fa Sequence_1.evol2.diverge.noGaps.fa

# We'll merge the diverged sequences and place them in the multiFa folder for validation later
~/go/bin/mergeMultiFa Sequence_0.evol1.fa Sequence_0.evol2.fa multiFa/Sequence_0.evol.fa
~/go/bin/mergeMultiFa Sequence_1.evol1.fa Sequence_1.evol2.fa multiFa/Sequence_1.evol.fa
rm Sequence_0.evol1.fa Sequence_0.evol2.fa Sequence_1.evol1.fa Sequence_1.evol2.fa
rm Sequence_0.fa Sequence_1.fa

# Now we test samAssembler for a series of sequencing libraries

set seed = `expr $seed + 1`
~/go/bin/simulateSam -coverage $coverage -flatErrorRate $simulateSamEpsilon -ancientErrorRate $simulateSamDeaminationRate -setSeed $seed -readLength $readLength -fragmentLength $fragmentLength evol1.diverge.noGaps.fa evol1.$readLength.$coverage.Reads.bam
set seed = `expr $seed + 1`
~/go/bin/simulateSam -coverage $coverage -flatErrorRate $simulateSamEpsilon -ancientErrorRate $simulateSamDeaminationRate -setSeed $seed -readLength $readLength -fragmentLength $fragmentLength evol2.diverge.noGaps.fa evol2.$readLength.$coverage.Reads.bam

$bedtools bamtofastq -i evol1.$readLength.$coverage.Reads.bam -fq evol1.$readLength.$coverage.readA.fastq -fq2 evol1.$readLength.$coverage.readB.fastq
$bedtools bamtofastq -i evol2.$readLength.$coverage.Reads.bam -fq evol2.$readLength.$coverage.readA.fastq -fq2 evol2.$readLength.$coverage.readB.fastq

cat evol1.$readLength.$coverage.readA.fastq evol2.$readLength.$coverage.readA.fastq > merged.readA.fastq
cat evol1.$readLength.$coverage.readB.fastq evol2.$readLength.$coverage.readB.fastq > merged.readB.fastq

$bwa mem -t 8 ref.fa merged.readA.fastq merged.readB.fastq | $samtools view -bh - > diverged.RefAln.bam
$samtools sort diverged.RefAln.bam -o diverged.RefAln.sorted.bam

## score gatk prediction
#mkdir -p output/gatk
#mkdir -p output/gatk/haplotypes
#$gatk AddOrReplaceReadGroups -I diverged.RefAln.sorted.bam -O diverged.RefAln.sorted.rg.bam -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM sample1 -CREATE_INDEX True
#$gatk HaplotypeCaller -I diverged.RefAln.sorted.rg.bam -R ref.fa -O output/gatk/$readLength.$startingSeed.$coverage.vcf
#~/go/bin/faFormat -noGapBed ref.noGap.bed ref.fa /dev/null
#~/go/bin/haplotypeGenerator -includeRef ref.fa output/gatk/$readLength.$startingSeed.$coverage.vcf  ref.noGap.bed output/gatk/haplotypes
#~/go/bin/mergeMultiFa output/gatk/haplotypes/Sequence_0.0.$chromSize.fa multiFa/Sequence_0.evol.fa Sequence_0.$readLength.$coverage.gatk.validate.fa
#~/go/bin/mergeMultiFa output/gatk/haplotypes/Sequence_1.0.$chromSize.fa multiFa/Sequence_1.evol.fa Sequence_1.$readLength.$coverage.gatk.validate.fa

~/go/bin/samAssembler prior -minCoverage 0 diverged.RefAln.sorted.bam ref.fa diverged.prior.txt
~/go/bin/samAssembler build -delta 0.05 -empiricalPrior diverged.prior.txt -multiFaDir empiricalMultiFa diverged.RefAln.sorted.bam ref.fa /dev/null /dev/null
~/go/bin/samAssembler build -delta 0.05 -multiFaDir conventionalMultiFa diverged.RefAln.sorted.bam ref.fa /dev/null /dev/null
~/go/bin/samAssembler build -delta 0.05 -multiFaDir lambdaMultiFa -lambda $lambda diverged.RefAln.sorted.bam ref.fa /dev/null /dev/null

~/go/bin/mergeMultiFa empiricalMultiFa/Sequence_0.fa multiFa/Sequence_0.evol.fa Sequence_0.$readLength.$coverage.empirical.validate.fa
~/go/bin/mergeMultiFa empiricalMultiFa/Sequence_1.fa multiFa/Sequence_1.evol.fa Sequence_1.$readLength.$coverage.empirical.validate.fa
~/go/bin/mergeMultiFa conventionalMultiFa/Sequence_0.fa multiFa/Sequence_0.evol.fa Sequence_0.$readLength.$coverage.conventional.validate.fa
~/go/bin/mergeMultiFa conventionalMultiFa/Sequence_1.fa multiFa/Sequence_1.evol.fa Sequence_1.$readLength.$coverage.conventional.validate.fa
~/go/bin/mergeMultiFa lambdaMultiFa/Sequence_0.fa multiFa/Sequence_0.evol.fa Sequence_0.$readLength.$coverage.lambda.validate.fa
~/go/bin/mergeMultiFa lambdaMultiFa/Sequence_1.fa multiFa/Sequence_1.evol.fa Sequence_1.$readLength.$coverage.lambda.validate.fa

ls -1 *.$readLength.$coverage.gatk.validate.fa > fileList.$readLength.$coverage.gatk.txt
ls -1 *.$readLength.$coverage.empirical.validate.fa > fileList.$readLength.$coverage.empirical.txt
ls -1 *.$readLength.$coverage.conventional.validate.fa > fileList.$readLength.$coverage.conventional.txt
ls -1 *.$readLength.$coverage.lambda.validate.fa > fileList.$readLength.$coverage.lambda.txt

~/go/bin/samAssembler score baseMatrix fileList.$readLength.$coverage.gatk.txt output/$readLength.$startingSeed.$coverage.gatk.txt
~/go/bin/samAssembler score baseMatrix fileList.$readLength.$coverage.empirical.txt output/$readLength.$startingSeed.$coverage.empirical.txt
~/go/bin/samAssembler score baseMatrix fileList.$readLength.$coverage.conventional.txt output/$readLength.$startingSeed.$coverage.conventional.txt
~/go/bin/samAssembler score baseMatrix fileList.$readLength.$coverage.lambda.txt output/$readLength.$startingSeed.$coverage.lambda.txt

rm evol1.$readLength.$coverage.Reads.bam evol2.$readLength.$coverage.Reads.bam evol1.$readLength.$coverage.readA.fastq evol1.$readLength.$coverage.readB.fastq
rm evol2.$readLength.$coverage.readA.fastq evol2.$readLength.$coverage.readB.fastq merged.readA.fastq merged.readB.fastq
rm fileList.$readLength.$coverage.*.txt
rm ref.fa.*
rm *.fa
rm *.bam
rm *.bai
#rm -r output/gatk

echo DONE
