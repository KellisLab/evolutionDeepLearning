#!/bin/csh -evf

mkdir -p divergenceAnalysis
set bigWigAverageOverBed = /home/rimangan/data/Software/x86_kentBin/bigWigAverageOverBed
set bw = /home/rimangan/data/reference/Conservation/haqer.divergenceDensity.bw
mkdir -p divergenceVcf
cat ~/data/homininAln/5wayPrimate/reconOut/substitutionsVcf/*.vcf.gz > divergenceVcf/allHumanHcaVariants.vcf.gz
set vcfFile = divergenceVcf/allHumanHcaVariants.vcf.gz
set intervalOverlap = ~/go/bin/intervalOverlap

mkdir -p divergenceAnalysis/inputBeds
mkdir -p divergenceAnalysis/outputTables
mkdir -p divergenceAnalysis/mutationCounts
foreach b (*.merged.bed)
	set bPrefix = $b:t:r
	awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $1"."$2"."$3}' $b > divergenceAnalysis/inputBeds/$bPrefix.bed
	#$bigWigAverageOverBed $bw divergenceAnalysis/inputBeds/$bPrefix.bed divergenceAnalysis/outputTables/$bPrefix.txt
	$intervalOverlap -mergedOutput divergenceAnalysis/inputBeds/$bPrefix.bed $vcfFile divergenceAnalysis/mutationCounts/$bPrefix.txt
end

echo DONE
