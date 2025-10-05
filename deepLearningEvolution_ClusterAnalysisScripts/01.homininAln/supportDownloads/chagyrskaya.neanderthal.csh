#!/bin/csh -evf

# note that these are aligned against hg19

sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr1-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr2-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr3-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr4-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr5-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr6-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr7-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr8-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr9-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr10-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr11-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr12-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr13-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr14-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr15-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr16-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr17-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr18-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr19-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr20-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr21-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr22-reali.bam"
sbatch --wrap="wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chrX-reali.bam"

echo DONE
