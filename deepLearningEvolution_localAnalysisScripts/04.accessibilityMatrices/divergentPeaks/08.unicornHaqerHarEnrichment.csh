#!/bin/csh -ef

mkdir -p unicornEnrichment
~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInHominins.merged.bed regionsOfInterest/UNICORNs.bed all.homininGreatApe.merged.bed unicornEnrichment/upInHominins.unicorn.bed
~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInGreatApes.merged.bed regionsOfInterest/UNICORNs.bed all.homininGreatApe.merged.bed unicornEnrichment/upInGreatApes.unicorn.bed

~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInHominins.merged.bed regionsOfInterest/HAR_Walsh_List_hg38.bed all.homininGreatApe.merged.bed unicornEnrichment/upInHominins.har.bed
~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInGreatApes.merged.bed regionsOfInterest/HAR_Walsh_List_hg38.bed all.homininGreatApe.merged.bed unicornEnrichment/upInGreatApes.har.bed

~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInHominins.merged.bed regionsOfInterest/haqer.ordered.bed all.homininGreatApe.merged.bed unicornEnrichment/upInHominins.haqer.bed
~/go/bin/overlapEnrichments -trimToSearchSpace -relationship any normalApproximate upInGreatApes.merged.bed regionsOfInterest/haqer.ordered.bed all.homininGreatApe.merged.bed unicornEnrichment/upInGreatApes.haqer.bed

cat unicornEnrichment/*.bed | grep -v 'Method' > unicornEnrichmentSummary.txt


echo DONE
