#!/bin/csh -ef
mkdir -p ancestralizedByChrom
mkdir -p chainTmp
mkdir -p chainOutput
foreach f (output/*fa.gz)
	set chrom = $f:t:r:r
	echo $chrom
	# first we pull out the sequence made by multiFaSequenceSwap and put it in another directory
	~/go/bin/faFilter -name $chrom.ancestral_HARs $f ancestralizedByChrom/$chrom.fa.gz

	# next we want to make a chain file to interoperate between hg38 coordinates and coordinates in the 
	# swapped sequence. First ,we make a multiFa that is pairwise human / ancestral_HARs, then we run
	# multiFaToChain	
	~/go/bin/faFilter -name hg38 $f chainTmp/$chrom.hg38.fa
	~/go/bin/faFilter -name $chrom.ancestral_HARs $f chaimpTmp/$chrom.ancestral_HARs.fa
	cat chainTmp/$chrom.hg38.fa chaimpTmp/$chrom.ancestral_HARs.fa > chainTmp/$chrom.pairwise.fa
	~/go/bin/multiFaToChain chainTmp/$chrom.pairwise.fa $chrom $chrom.ancestral_HARs chainOutput/$chrom.chain
end
zcat ancestralizedByChrom/*.fa.gz > hg38.ancestralizedHars.fa
zcat chainOutput/*.chain > hg38.liftOver.hg38_ancestralizedHars.chain

samtools faidx hg38.ancestralizedHars.fa


echo done
