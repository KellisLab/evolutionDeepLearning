#!/bin/csh -ef


mkdir -p enhancerPromoterSplits

foreach b (*.merged.bed)
	set bPrefix = $b:t:r
	~/go/bin/intervalOverlap -nonOverlap ~/data/reference/epimap/chromatinStates/stateUnionSets/TssA.bed $b enhancerPromoterSplits/$bPrefix.TssANonOverlap.bed
	~/go/bin/intervalOverlap ~/data/reference/epimap/chromatinStates/stateUnionSets/TssA.bed $b enhancerPromoterSplits/$bPrefix.TssA.bed
	~/go/bin/intervalOverlap ~/data/reference/epimap/chromatinStates/regulatoryElements.bed enhancerPromoterSplits/$bPrefix.TssANonOverlap.bed enhancerPromoterSplits/$bPrefix.OtherRegulatoryStates.bed
	~/go/bin/intervalOverlap -nonOverlap ~/data/reference/epimap/chromatinStates/regulatoryElements.bed enhancerPromoterSplits/$bPrefix.TssANonOverlap.bed enhancerPromoterSplits/$bPrefix.nonRegulatory.bed
	rm enhancerPromoterSplits/$bPrefix.TssANonOverlap.bed
end
