# BigWig VIZ

This directory contains a few files for visualizing bigWig format genomic data. These leverage the pyBigWig package for python.
Again, paths are based on our cluster, so won't work out of the box. Essentially, you'll need one folder with prediction bigWigs, and another with assembly coordinate mappings. These are critical as ANCoRA assemblies each have their own assembly system, so we use one set of coords (hg38) as global coords, and convert to the correct coords for each genome.

Happy to answer any questions at rimangan@mit.edu

Best,
Riley
