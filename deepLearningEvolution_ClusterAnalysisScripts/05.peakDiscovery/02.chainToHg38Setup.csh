#!/bin/csh -ef

set chainSwap = /home/rimangan/data/Software/x86_kentBin/chainSwap
set chainDir = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/chainOutput

mkdir -p toHg38Chains

foreach chainFile ($chainDir/*.chain)
	set chainPrefix = $chainFile:t:r
	$chainSwap $chainFile toHg38Chains/$chainPrefix.To.Hg38.chain
end

echo DONE
