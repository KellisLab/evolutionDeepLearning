#!/bin/csh -ef

set faSize = /net/bmc-lab4/data/kellis/users/rimangan/Software/x86_kentBin/faSize

mkdir -p chromSizesOutput

foreach b (buildOutput/*.fa)
	set bPrefix = $b:t:r
	#$faSize -detailed $b > chromSizesOutput/$bPrefix.chrom.sizes
	samtools faidx $b
end

echo DONE
