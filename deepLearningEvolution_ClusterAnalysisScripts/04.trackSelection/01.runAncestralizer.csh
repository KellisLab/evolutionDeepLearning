#!/bin/csh -ef

# note: these paths are from our cluster. You'll need to update with your own. 
# HARs are not provided, but they are available from GEO:GSE180714

set reconDir = /net/bmc-lab4/data/kellis/users/rimangan/homininAln/5wayPrimate/reconOut
set HarFile  = /net/bmc-lab4/data/kellis/users/rimangan/reference/regionsOfInterest/HARs.GSE180714.filtered.bed

#reconDir consists of a directory with reconstructed genomes, HarFile consists of HARs (bed file)


mkdir -p harsByChrom #-p is used in mkdir to ensure that no error occurs if such a director already exists
mkdir -p output

~/go/bin/bedSplit byChrom $HarFile harsByChrom #the "bedSplit" command is called, byChrom is an input to the funcion
#which tells it to split the incoming "HarFile" by chromosome, these results (new file for every chromosome) are then put in a directory called harsByChrom
set background = hg38
set foreground = hca
foreach chromFile (harsByChrom/*.bed)
	set chrom = $chromFile:t:r #tail, path removed and root, file extension name removed
	echo $chrom
	sbatch --wrap="~/go/bin/multiFaSequenceSwap $reconDir/$chrom.fa.gz $chromFile $background $foreground output/$chrom.fa.gz $chrom.ancestral_HARs"
end #several files are made, corresponding to each chromosome swapped
	









echo done
