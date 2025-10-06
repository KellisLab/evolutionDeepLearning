#!/bin/csh -ef

# === INPUT VARIABLES ===
#set track = 570 #heart_embryo_80_days
#set bedFile = heart_embryo_80_days.chr1.201381090.201381626.bed

#set track = 347 #hepatocyte from H9
#set bedFile = chr1.201381602.201382791.bed

#set track = 0 #cerebellum_male_adult_27_years_35_years
#set bedFile = bed/chr19.1795882.1796497.bed #GRIN3B

#set track = 179 # brain_female_embryo_142_days
#set bedFile = bed/chr6.163756292.163756452.bed #QKI

#set track = 370 # brain_male_embryo_105_days
#set bedFile = bed/chr6.11236447.11236874.bed #NEDD9

#set track = 77 # astrocyte of the spinal cord
#set bedFile = bed/chr15.40266842.40267256.bed #EIF2AK4

#set track = 370
#set bedFile = bed/chr12.1933407.1934070.bed #CACNA1C

set track = 370
set bedFile = bed/chr10.86261934.86262389.bed #GRID1

# === CONSTANT PATHS ===
set enformerPath = ~/KellisProjects/Enformer/enformer_1/
set genomeFasta = /Users/dickinsonia/Scavenger/Datasets/Genomes/hg38/hg38.fa
set outPrefix = `basename $bedFile .bed`
set outFile = satMutResults/${outPrefix}.track.${track}.txt

# === COMMAND ===
python 01.satMutGen.py \
  -e $enformerPath \
  -f $genomeFasta \
  -b $bedFile \
  -o $outFile \
  -t $track

echo DONE
