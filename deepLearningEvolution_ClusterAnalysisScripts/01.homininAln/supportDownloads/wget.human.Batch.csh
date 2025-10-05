#!/bin/csh -evf

# HG02610 - Gambian Mandinka male
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242417/HG02610.final.cram"

# HG02615 - Gambian Mandinka female
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989037/HG02615.final.cram"

# HG00513 - Han female
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241684/HG00513.final.cram"

# HG03055 - Mende female
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242541/HG03055.final.cram"

# NA18860 - Yoruba male
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989358/NA18860.final.cram"

# NA18858 - Yoruba female
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239384/NA18858.final.cram"

# NA20505 - Toscani female
#sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239788/NA20505.final.cram"

# HG00525 - Han Female (#2)
sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241686/HG00525.final.cram"

# HG03079 - Mende Female (#2)
sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242958/HG03079.final.cram"

# NA20517 - Toscani Female (#2)
sbatch --wrap="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239800/NA20517.final.cram"

echo DONE
