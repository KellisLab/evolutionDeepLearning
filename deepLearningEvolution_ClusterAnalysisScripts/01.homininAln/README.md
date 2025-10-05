# Hominin Alignment Directory

This directory contains all scripts for generating the 35-way alignment of great apes used throughout the study. 
This includes:
(1) accessing raw reads from the web
(2) aligning reads to hg38
(3) Using ANCoRA to generate reference guided assemblies from BAMs
(4) Aligning the assemblies together

These sequence files are considerably large, and web-accessible. So while analysis scripts are present here for public access, the raw libraries are not. Instead, download links are provided.

## Important Note

The raw data files used in this study were a bit heterogenous. Some files were in CRAM format, others BAM to hg19, others FASTQ. I'm assuming that for each file, you are unalizning the FASTQ and then realigning the hg38. This is time-consuming, but best practice for uniform processing. 

Given this heterogeneity, I haven't written a 'plug-and-play' script that can process all of them. In reality, I processed each file one at a time. So for example:

For jimmi.panTroVerus (a chimp sample, presumably from an individual named Jimmi), there were several files all sequenced from the same individual. To maximize coverage, I used all files to construct a final BAM. These files were in sra, so I used fastq-dump to get final fastq files before progressing to hg38 alignments.

If libraries are spread across multiple fastq files, I found it was convenient to use a slurm array to run bwa.

While I believe there is enough detail here to reconstruct our process, if you are either (1) trying to replicate or (2) using this for your own alignments, let me know if you have questions.

Sincerely,
Riley

rimangan@mit.edu
