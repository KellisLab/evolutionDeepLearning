#!/bin/csh -evf

# here we'll insert Olivia's multiFaToVcf code when complete


# now we generate a chain

set target = hg38
mkdir -p chainOutput

foreach currDir (buildOutput/*)
	if (-d $currDir) then
		set dirPrefix = $currDir:t:r
		echo $dirPrefix
		mkdir -p $currDir/split
		mkdir -p chainOutput/$dirPrefix
		foreach chromFile ($currDir/*.fa)
			echo $chromFile
			set chrom = $chromFile:t:r
			echo $chrom
			~/go/bin/faFilter -notName $dirPrefix.B $chromFile $currDir/split/$chrom.$dirPrefix.QueryA.fa
			~/go/bin/faFilter -notName $dirPrefix.A $chromFile $currDir/split/$chrom.$dirPrefix.QueryB.fa
			~/go/bin/multiFaToChain $currDir/split/$chrom.$dirPrefix.QueryA.fa $chrom $chrom chainOutput/$dirPrefix/$dirPrefix.A.$chrom.chain
			~/go/bin/multiFaToChain $currDir/split/$chrom.$dirPrefix.QueryB.fa $chrom $chrom chainOutput/$dirPrefix/$dirPrefix.B.$chrom.chain
		end
		cat chainOutput/$dirPrefix/$dirPrefix.A.*.chain > chainOutput/$dirPrefix.A.merged.chain
		cat chainOutput/$dirPrefix/$dirPrefix.B.*.chain > chainOutput/$dirPrefix.B.merged.chain
	endif
end

echo DONE
