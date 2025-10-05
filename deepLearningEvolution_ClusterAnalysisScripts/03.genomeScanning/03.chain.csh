#!/bin/csh -ef

mkdir -p chainOutput

foreach seq (`cat haploNames.txt`)
	set name = $seq:r
	foreach chromFile (*.merged.fa)
		set chrom = $chromFile:t:r:r
		~/go/bin/multiFaToChain -querySeqName $seq $chromFile $chrom $chrom chainOutput/$chrom.$seq.chain
	end
	cat chainOutput/*.$seq.chain > chainOutput/$seq.merged.chain
	rm chainOutput/*.$seq.chain
end

echo DONE
