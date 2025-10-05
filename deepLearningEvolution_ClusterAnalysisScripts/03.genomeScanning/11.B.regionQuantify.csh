#!/bin/csh -ef

set chainDir = ~/data/ancoraUpdated/chainOutput
set bigWigDir = ~/data/ancoraUpdated/bigWigOutput
set regionsBed = /net/bmc-lab4/data/kellis/users/rimangan/reference/regionsOfInterest/haqer.hg38.bed
set liftedRegionsDir = ~/data/ancoraUpdated/liftedRegions
set bigWigAverageOverBed = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/bigWigAverageOverBed

mkdir -p regionQuantificationOutput
set regionsPrefix = $regionsBed:t:r:r
echo $regionsPrefix

foreach chain ($chainDir/*.B.merged.chain)
        set assembly = $chain:t:r:r
        foreach bw ($bigWigDir/$assembly/*.bw)
                echo $bw
                set bwPrefix = $bw:t:r
                $bigWigAverageOverBed $bw liftedRegions/$regionsPrefix.$assembly.bed regionQuantificationOutput/$bwPrefix.$regionsPrefix.$assembly.txt
        end
end


echo DONE

