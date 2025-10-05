#!/bin/csh -evf

set trackFile = brainDnaseTracks.txt
set humanGenome = ~/data/reference/hg38.fa
set ancestralizedGenome = hg38.ancestralizedHars.fa
set harFile = /net/bmc-lab4/data/kellis/users/rimangan/reference/regionsOfInterest/HARs.GSE180714.filtered.bed
set liftOver = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/liftOver 

#python amplitudeAnalysis.py -f $humanGenome -c human -r $harFile -t $trackFile -o hg38.enformerResults.txt

set chain = hg38.liftOver.hg38_ancestralizedHars.chain
$liftOver $harFile $chain hars.positionsOnAncestralizedAssembly.bed unmapped.txt
python amplitudeAnalysis.py -f $ancestralizedGenome -c ancestralized -r hars.positionsOnAncestralizedAssembly.bed -t $trackFile -o ancestralizedHars.enformerResults.txt

echo DONE
