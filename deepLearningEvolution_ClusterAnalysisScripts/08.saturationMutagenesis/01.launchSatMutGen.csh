#!/bin/csh -ef

#set bed = chr20.45095781.45096191.bed
#set track = 370 # brain_male_embryo_105_days

#set bed = chr1.201381602.201382791.bed
#set track = 570 # heart embryo 80 days

#set bed = chr1.201381602.201382791.bed
#set track = 347 # hepatocyte from h9



set bPrefix = $bed:t:r
mkdir -p satMutResults
python support.01.satMutGen.py -t $track -e ../enformer_1/ -f /home/rimangan/data/reference/hg38.fa -b $bed -o satMutResults/$bPrefix.track.$track.satMut.txt

echo DONE
