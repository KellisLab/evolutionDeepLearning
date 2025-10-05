#!/bin/csh -ef

set chainDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/chainOutput
set regionsBed = /net/bmc-lab4/data/kellis/users/rimangan/reference/regionsOfInterest/haqer.hg38.bed
set liftOver = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/liftOver

mkdir -p liftedRegions

set regionsPrefix = $regionsBed:t:r:r
echo $regionsPrefix

foreach chain ($chainDir/*.B.merged.chain)
        set assembly = $chain:t:r:r
        $liftOver $regionsBed $chain liftedRegions/$regionsPrefix.$assembly.bed /dev/null
end

echo DONE
