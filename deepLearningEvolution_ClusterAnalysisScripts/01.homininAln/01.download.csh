#!/bin/csh -ef

#make sure to run these first!
#module load perl/5.24.1
#module load sra/2.10.0

#reads from primate diversity project are found here: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR748156&display=data-access


echo "Stopping the program early, as it can start a massive, multigigabyte download. If you know what you are doing, take a look in the file"
exit 0

# I've organized this so each genome has its own download script, which are organized in the supportDownloads directory. You can add them like this:
#./supportDownloads/taweh.panTroEllioti.csh


# IMPORTANT NOTE:
# Please read the README before making changes here


echo DONE
