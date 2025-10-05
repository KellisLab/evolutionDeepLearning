#!/bin/csh

set split = $1
set splitPrefix = $split:t:r:r
set trackFile = allDnaseTracks.FixedNames.txt
set humanGenome = ~/data/reference/hg38.fa
set ancestralizedGenome = hg38.ancestralizedHars.fa
set hcaGenome = ~/data/reference/hca/hca.fa
set harFile = $split
set liftOver = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/liftOver
mkdir -p enformerResults

# Enformer predictions on hg38
python amplitudeAnalysis.py -f $humanGenome -c human -r $harFile -t $trackFile -o enformerResults/hg38.$splitPrefix.enformerResults.txt

# Enformer predictions on hg38.Ancestral_HARs
set chain = hg38.liftOver.hg38_ancestralizedHars.chain
$liftOver $harFile $chain harSplits/$splitPrefix.positionsOnAncestralizedAssembly.bed unmapped.txt
python amplitudeAnalysis.py -f $ancestralizedGenome -c ancestralized -r harSplits/$splitPrefix.positionsOnAncestralizedAssembly.bed -t $trackFile -o enformerResults/ancestralizedHars.$splitPrefix.enformerResults.txt

# Enformer predictions on hca
set chain = /home/rimangan/data/reference/liftOver/hg38ToHca.over.chain.gz
$liftOver $harFile $chain harSplits/$splitPrefix.positionsOnHca.bed unmapped.txt
python amplitudeAnalysis.py -f $hcaGenome -c hca -r harSplits/$splitPrefix.positionsOnHca.bed -t $trackFile -o enformerResults/hcaHars.$splitPrefix.enformerResults.txt

echo DONE
