#!/bin/csh -ef

set chrom = $1

cp ~/data/ancoraUpdated/buildOutput/HG02615.multiFa/$chrom.fa ~/data/ancoraUpdated/$chrom.merged.fa

foreach seq (`cat seqNames.txt`)
	set multiFaDir = ~/data/ancoraUpdated/buildOutput/$seq.multiFa
	~/go/bin/mergeMultiFa $chrom.merged.fa $multiFaDir/$1.fa $chrom.tmp.fa
	mv $chrom.tmp.fa $chrom.merged.fa
end

echo DONE
	
