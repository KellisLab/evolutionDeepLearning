#!/bin/csh -ef

set liftOver = "/Users/dickinsonia/Software/KentBin/liftOver"
set chainFile = "/Users/dickinsonia/Scavenger/Datasets/Genomes/ChainFiles/hg19ToHg38.over.chain"

mkdir -p formattedBeds

# script expects raw to contain the files from: https://personal.broadinstitute.org/cboix/epimap/links/pergroup/

python support.01.convertToBed.py --input_dir raw --gtf ~/Scavenger/Datasets/Genomes/hg38/genes/gencode.v44.basic.annotation.gtf.gz --output_dir formattedBeds/

foreach b (formattedBeds/*.bed)
    set bPrefix = $b:t:r
    $liftOver $b $chainFile formattedBeds/$bPrefix.hg38.bed /dev/null
    rm $b
end