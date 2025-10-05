#!/bin/csh 

mkdir -p harSplits
set harBed = /net/bmc-lab4/data/kellis/users/rimangan/reference/regionsOfInterest/HARs.GSE180714.filtered.bed
split -l 160 $harBed harSplits/partition_

@ i = 1
foreach file (harSplits/partition_*)
    set new_name = "harSplits/partition_$i.hg38.bed"
    mv $file $new_name
    @ i++
end

echo DONE
