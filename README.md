# evolutionDeepLearning
Analysis scripts associated with "Deep learning predicts cis-regulatory turnover in human evolution". Mangan et al, BioRxiv 2025.


Much of this work was done back and forth between local machines and the Luria compute cluster here at MIT, and analysis files are represented in the two directories, respectively.

I've added in a fair bit of annotation, but certainly not complete. 


Some required software:
* Gonomics: a suite of bioinformatics software in the go programming language. Available from: 
	* https://github.com/vertgenlab/gonomics
* Enformer: a pre-trained sequence-to-function model used for accessibility predictions. 
	* https://www.nature.com/articles/s41592-021-01252-x
* Bedtools: For manipulating genomic intervals.
	* https://bedtools.readthedocs.io/en/latest/
* Samtools: For alignment, BAM/CRAM handling.
	* https://www.htslib.org/ 
* KentUtils: UCSC utilities.
	* https://github.com/ENCODE-DCC/kentUtils

Please reach out with any questions:
Riley J. Mangan
rimangan@mit.edu
