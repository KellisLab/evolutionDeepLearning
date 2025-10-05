#!/bin/csh -evf

set sizes = ~/data/reference/hg38.chrom.sizes
set wigToBigWig = /home/rimangan/data/Software/x86_kentBin/wigToBigWig
mkdir -p coverageTrackHub

foreach w (coverageWigs/*.wig.gz)
	set wPrefix = $w:t:r:r
	~/go/bin/wigFilter -chrom chr10 $w $sizes coverageWigs/$wPrefix.chr10.wig.gz
	$wigToBigWig coverageWigs/$wPrefix.chr10.wig.gz $sizes coverageWigs/$wPrefix.chr10.bw
end

# now we make a trackhub to display coverageWigs
set thDir = chr10CoverageTrackHub
mkdir -p $thDir
set urlPrefix = "https://vertgenlab.cloud.duke.edu/riley/$thDir/hg38" 

#make hub.txt file
echo "hub coverageWigs" > $thDir/hub.txt
echo "shortLabel coverageWigs" >> $thDir/hub.txt
echo "longLabel Wigs of coverage depth" >> $thDir/hub.txt
echo "email rimangan@mit.edu" >> $thDir/hub.txt
echo "genomesFile genomes.txt" >> $thDir/hub.txt

#make genomes.txt
echo "genome hg38" > $thDir/genomes.txt
echo "trackDb trackDb.txt" >> $thDir/genomes.txt
mkdir -p $thDir/hg38

#make trackDb.txt
touch $thDir/trackDb.txt
rm $thDir/trackDb.txt

foreach bw (coverageWigs/*.bw)
	set track = $bw:t:r
	echo "track $track" >> $thDir/trackDb.txt
	echo "type bigWig" >> $thDir/trackDb.txt
	echo "bigDataUrl $urlPrefix/$track.bw" >> $thDir/trackDb.txt
	echo "shortLabel $track" >> $thDir/trackDb.txt
	echo "longLabel $track" >> $thDir/trackDb.txt
	echo "visibiility full" >> $thDir/trackDb.txt
	echo "autoScale on" >> $thDir/trackDb.txt
	echo "alwaysZero on" >> $thDir/trackDb.txt
	echo "" >> $thDir/trackDb.txt
end

echo DONE
