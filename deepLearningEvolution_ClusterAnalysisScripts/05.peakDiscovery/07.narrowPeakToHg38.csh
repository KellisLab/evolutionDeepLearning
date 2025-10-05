#!/bin/csh -ef

set liftOver = /home/rimangan/data/Software/x86_kentBin/liftOver

mkdir -p narrowPeakHg38

foreach chain (toHg38Chains/*.merged.To.Hg38.chain)
        set haplotype = $chain:t:r:r:r:r
        echo $haplotype
        foreach narrow (narrowPeak/$haplotype.*.narrowPeak)
                echo $narrow
                set narrowPrefix = $narrow:t:r
                
                # Strip header line before running liftOver
                grep -v "^track" $narrow | sort -k5,5nr > temp_noHeader.bed
                $liftOver temp_noHeader.bed $chain narrowPeakHg38/$narrowPrefix.hg38.narrowPeak /dev/null
                rm temp_noHeader.bed
	end
end

echo DONE
