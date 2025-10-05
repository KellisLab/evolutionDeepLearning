# Genome Scanning

This folder contains analysis scripts related to scanning the genome to construct genome-wide accessibility prediction profiles. 

Note that there is some overlap with 01.homininAln, which provides information on how to construct a multiple alignment from reads.
In this case, we are assuming we are starting with a directory (or two directories) of BAM files, which have been uniformly aligned to hg38. This is an ideal starting point for this pipeline. We note that in the script, we have preserved our local paths from our clusters, which won't match your machine. If you are building your own alignments with our pipeline, you can simply refactor 01.build.csh to point to your directory with finished bam files, one file per diploid organism.

As ever, please reach out with any questions.

Best,
Riley

rimangan@mit.edu
