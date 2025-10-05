#!/bin/csh -ef

set bigWigDir = $1

set bigWigToBedGraph = /home/rimangan/data/Software/x86_kentBin/bigWigToBedGraph
mkdir -p narrowPeak
set outDir = narrowPeak
mkdir -p cutoffAnalysis

foreach bigWigFile ($bigWigDir/*.bw)
	set bwPrefix = $bigWigFile:t:r
	$bigWigToBedGraph $bigWigFile $bwPrefix.bedGraph
	set peakOutput = $outDir/$bwPrefix.narrowPeak
	singularity exec -B /net/bmc-lab4/data/kellis/users/rimangan/enformer/macs3/ /home/software/macs3/macs3.sif \
		 macs3 bdgpeakcall -i $bwPrefix.bedGraph --cutoff-analysis --cutoff-analysis-steps 100 -o cutoffAnalysis/$bwPrefix.cutoff.txt
	set cutoff = `python findElbow.py cutoffAnalysis/$bwPrefix.cutoff.txt`
	singularity exec -B /net/bmc-lab4/data/kellis/users/rimangan/enformer/macs3/ /home/software/macs3/macs3.sif macs3 bdgpeakcall -i $bwPrefix.bedGraph -o narrowPeak/$bwPrefix.$cutoff.narrowPeak  -c $cutoff
	rm $bwPrefix.bedGraph
end

echo DONE
