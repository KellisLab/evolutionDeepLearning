#!/bin/csh -ef

cat upInArchaic/*.bed > upInArchaic.unmerged.bed
~/go/bin/bedMerge upInArchaic.unmerged.bed upInArchaic.merged.bed
cat upInModern/*.bed > upInModern.unmerged.bed
~/go/bin/bedMerge upInModern.unmerged.bed upInModern.merged.bed
cat upInHominins/*.bed > upInHominins.unmerged.bed
~/go/bin/bedMerge upInHominins.unmerged.bed upInHominins.merged.bed
cat upInGreatApes/*.bed > upInGreatApes.unmerged.bed
~/go/bin/bedMerge upInGreatApes.unmerged.bed upInGreatApes.merged.bed

cat upInArchaic.unmerged.bed upInModern.unmerged.bed > all.modernArchaic.unmerged.bed
~/go/bin/bedMerge all.modernArchaic.unmerged.bed all.modernArchaic.merged.bed
cat upInHominins.unmerged.bed upInGreatApes.unmerged.bed > all.homininGreatApe.unmerged.bed
~/go/bin/bedMerge all.homininGreatApe.unmerged.bed all.homininGreatApe.merged.bed

echo DONE