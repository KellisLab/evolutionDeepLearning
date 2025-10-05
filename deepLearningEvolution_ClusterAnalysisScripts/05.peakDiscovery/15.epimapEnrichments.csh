#!/bin/csh -ef

set states = (EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF_Rpts)

foreach state ($states)
	sbatch --wrap="./support.15.epimapEnrichments.csh $state"
end

echo DONE
