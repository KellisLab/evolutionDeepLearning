#!/bin/csh -evf

set sizes = /home/rimangan/data/reference/hg38.chrom.sizes

mkdir -p mergedMultiFa

foreach chrom (`cat $sizes | cut -f1`)
	echo $chrom
	set firstTime = 1
	foreach currDir (buildOutput/*)
		echo $currDir
		if (-d $currDir) then
			if ($firstTime == 1) then
				set currDirPrefix = $currDir:t:r
				echo $currDirPrefix
				set firstTime = 0
				cp $currDir/$chrom.fa mergedMultiFa/$chrom.merged.fa
			else
				~/go/bin/mergeMultiFa mergedMultiFa/$chrom.merged.fa $currDir/$chrom.fa mergedMultiFa/tmp.$chrom.fa.gz
				mv mergedMultiFa/tmp.$chrom.fa.gz mergedMultiFa/$chrom.merged.fa.gz
			endif
		endif
	end
end

echo DONE
