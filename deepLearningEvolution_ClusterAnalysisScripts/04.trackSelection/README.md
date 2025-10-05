# Track Selection

This directory contains scripts related to selecting Enformer accessibility tracks for analysis, analysis displayed in Figure S4.
Note that you'll need to start with a 5-way alignment of great apes, available as a part of the data release from 
Mangan et al, 2022 (pmid36423581) on vertgenlab.org. https://www.vertgenlab.org/dataPages/2022_HAQERs.html

Once you have downloaded the multiple alignment, place 00.reconstruct.csh in the folder with the .fa.gz files to run reconstructSeq (gonomics) which will producce the necessary alignment in a subdir called reconOut.

We again assume gonomics is installed in ~/go/bin, and gonomics is available from https://github.com/vertgenlab/gonomics

Please reach out if you have questions.

Best,
Riley

rimangan@mit.edu




