#!/bin/csh -ef

set regionDir = "../accessibilityMatrices/divergentPeaks/upInHominins"
cat $regionDir/astrocyte_of_the_spinal_cord.upInHominins.bed \
$regionDir/bipolar_neuron_from_GM23338_doxycycline_4_days.upInHominins.bed \
$regionDir/brain_female_embryo_142_days.upInHominins.bed $regionDir/brain_male_embryo_105_days.upInHominins.bed \
$regionDir/brain_pericyte.upInHominins.bed $regionDir/cerebellum_male_adult_27_years_35_years.upInHominins.bed \
$regionDir/inferior_parietal_cortex_male_adult_84_years.upInHominins.bed $regionDir/medulloblastoma.upInHominins.bed \
$regionDir/neural_stem_progenitor_cell_from_H1_hESC.upInHominins.bed > tmp.bed
~/go/bin/bedMerge tmp.bed brain.upInHominins.merged.bed
rm tmp.bed

set modernDir = "../accessibilityMatrices/divergentPeaks/upInModern"
cat $modernDir/astrocyte_of_the_spinal_cord.upInModern.bed \
$modernDir/bipolar_neuron_from_GM23338_doxycycline_4_days.upInModern.bed \
$modernDir/brain_female_embryo_142_days.upInModern.bed $modernDir/brain_male_embryo_105_days.upInModern.bed \
$modernDir/brain_pericyte.upInModern.bed $modernDir/cerebellum_male_adult_27_years_35_years.upInModern.bed \
$modernDir/inferior_parietal_cortex_male_adult_84_years.upInModern.bed $modernDir/medulloblastoma.upInModern.bed \
$modernDir/neural_stem_progenitor_cell_from_H1_hESC.upInModern.bed > tmp.bed
~/go/bin/bedMerge tmp.bed brain.upInModern.merged.bed
rm tmp.bed


set homininRegions = "brain.upInHominins.merged.bed"
set modernRegions = "brain.upInModern.merged.bed"
set discordantRegions = "../accessibilityMatrices/divergentPeaks/discordant.bed"

mkdir -p mergedOutput

foreach b (formattedBeds/*.bed)
	set bPrefix = $b:t:r
	~/go/bin/intervalOverlap -mergedOutput $homininRegions $b mergedOutput/brain.upInHominin.$bPrefix.links.txt
	~/go/bin/intervalOverlap -mergedOutput $discordantRegions $b mergedOutput/discordant.$bPrefix.links.txt
	~/go/bin/intervalOverlap -mergedOutput $modernRegions $b mergedOutput/brain.upInModern.$bPrefix.links.txt
end

echo DONE
