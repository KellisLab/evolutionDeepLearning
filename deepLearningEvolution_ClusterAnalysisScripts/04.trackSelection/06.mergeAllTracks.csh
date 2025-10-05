#!/bin/csh

cat enformerResults/*.txt | grep -v 'Amplitude' | gzip > allTracks.enformerResults.merged.txt.gz

echo DONE
