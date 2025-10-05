# ANCoRA Validation Directory

This directory contains scripts used in Figure S1 and Figure 1B for validating the performance of 
ANCoRA, which we used to construct reference-guided personalized assemblies.

The critical working script is 01.validateWithPrior.csh, which calls the helper script support.validateWithPrior.csh.
There's some required software (bwa, samtools, bedtools), and variables in the 01 script to specify paths (right at the top of the script).
Also, this script assumes you have gonomics executables installed at ~/go/bin. You can change relevant portions if you store gonomics executables elsewhere. Gonomics is available at https://github.com/vertgenlab/gonomics




