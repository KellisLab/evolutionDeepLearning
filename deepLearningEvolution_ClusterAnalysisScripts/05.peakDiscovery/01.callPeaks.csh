#!/bin/csh -ef


set bigWigDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/bigWigOutput

foreach dir ($bigWigDir/*.B)
    if (-d $dir) then
        echo $dir
        sbatch --wrap="./support.01.callPeaks.csh $dir"
    endif
end

echo DONE
