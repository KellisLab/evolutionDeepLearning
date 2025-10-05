#!/bin/csh -ef

set buildOutput = /net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/buildOutput
mkdir -p vcf
foreach multiFaDir ($buildOutput/*.multiFa)
	echo $multiFaDir
	set sample = $multiFaDir:t:r
	mkdir -p vcf/$sample
	foreach chromFile ($multiFaDir/*.fa)
		set chrom = $chromFile:t:r
		~/go/bin/multiFaToVcf -substitutionsOnly $chromFile $chrom vcf/$sample/$chrom.vcf.gz
	end
	cat vcf/$sample/*.vcf.gz > vcf/$sample.merged.vcf.gz
end

echo DONE
