library(ggplot2)
library(ggpubr)
library(dplyr)
data <- read.csv("rawMpraResults.txt", header=TRUE, sep="\t")
metadata <- read.csv("metadata.txt", header=TRUE, sep="\t")
metadata <- select(metadata, c("hCONDEL_ID", "hg38_mpra_chr", "hg38_mpra_start_pos", "hg38_mpra_end_pos"))
data <- merge(data, metadata, by="hCONDEL_ID")

data$minusLog10pAdjSkewK562 <- -log10(data$padj_Skew_K562)
sum(data$padj_Skew_K562 < 0.05 | 
      data$padj_Skew_GM12878 < 0.05 | 
      data$padj_Skew_HEK293 < 0.05 | 
      data$padj_Skew_NPC < 0.05 | 
      data$padj_Skew_SKNSH < 0.05 | data$padj_Skew_HEPG2 < 0.05, na.rm=TRUE)
ggscatter(data, "log2FoldChange_Skew_K562", "minusLog10pAdjSkewK562")

# K562
upK562 <- subset(data, data$padj_Skew_K562 < 0.05 & 
                   (data$padj_Human_K562 < 0.1 | data$padj_Chimp_K562 < 0.1) & 
                   data$log2FoldChange_Skew_K562 > 0)
downK562 <- subset(data, data$padj_Skew_K562 < 0.05 & 
                     (data$padj_Human_K562 < 0.1 | data$padj_Chimp_K562 < 0.1) & 
                     data$log2FoldChange_Skew_K562 < 0)
# GM12878
upGM12878 <- subset(data, data$padj_Skew_GM12878 < 0.05 & 
                      (data$padj_Human_GM12878 < 0.1 | data$padj_Chimp_GM12878 < 0.1) & 
                      data$log2FoldChange_Skew_GM12878 > 0)
downGM12878 <- subset(data, data$padj_Skew_GM12878 < 0.05 & 
                        (data$padj_Human_GM12878 < 0.1 | data$padj_Chimp_GM12878 < 0.1) & 
                        data$log2FoldChange_Skew_GM12878 < 0)
# HEK293
upHEK293 <- subset(data, data$padj_Skew_HEK293 < 0.05 & 
                     (data$padj_Human_HEK293 < 0.1 | data$padj_Chimp_HEK293 < 0.1) & 
                     data$log2FoldChange_Skew_HEK293 > 0)
downHEK293 <- subset(data, data$padj_Skew_HEK293 < 0.05 & 
                       (data$padj_Human_HEK293 < 0.1 | data$padj_Chimp_HEK293 < 0.1) & 
                       data$log2FoldChange_Skew_HEK293 < 0)
# NPC
upNPC <- subset(data, data$padj_Skew_NPC < 0.05 & 
                  (data$padj_Human_NPC < 0.1 | data$padj_Chimp_NPC < 0.1) & 
                  data$log2FoldChange_Skew_NPC > 0)
downNPC <- subset(data, data$padj_Skew_NPC < 0.05 & 
                    (data$padj_Human_NPC < 0.1 | data$padj_Chimp_NPC < 0.1) & 
                    data$log2FoldChange_Skew_NPC < 0)
# SKNSH
upSKNSH <- subset(data, data$padj_Skew_SKNSH < 0.05 & 
                    (data$padj_Human_SKNSH < 0.1 | data$padj_Chimp_SKNSH < 0.1) & 
                    data$log2FoldChange_Skew_SKNSH > 0)
downSKNSH <- subset(data, data$padj_Skew_SKNSH < 0.05 & 
                      (data$padj_Human_SKNSH < 0.1 | data$padj_Chimp_SKNSH < 0.1) & 
                      data$log2FoldChange_Skew_SKNSH < 0)
# HEPG2
upHEPG2 <- subset(data, data$padj_Skew_HEPG2 < 0.05 & 
                    (data$padj_Human_HEPG2 < 0.1 | data$padj_Chimp_HEPG2 < 0.1) & 
                    data$log2FoldChange_Skew_HEPG2 > 0)
downHEPG2 <- subset(data, data$padj_Skew_HEPG2 < 0.05 & 
                      (data$padj_Human_HEPG2 < 0.1 | data$padj_Chimp_HEPG2 < 0.1) & 
                      data$log2FoldChange_Skew_HEPG2 < 0)


writeBed <- function(data, filename, nameField) {
  bed_data <- data.frame(
    chr = data$hg38_mpra_chr,
    start = data$hg38_mpra_start_pos,
    end = data$hg38_mpra_end_pos,
    name = nameField
  )
  write.table(bed_data, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

writeBed(upK562, "subSets/upK562.bed", "upK562")
writeBed(downK562, "subSets/downK562.bed", "downK562")
writeBed(upGM12878, "subSets/upGM12878.bed", "upGM12878")
writeBed(downGM12878, "subSets/downGM12878.bed", "downGM12878")
writeBed(upHEK293, "subSets/upHEK293.bed", "upHEK293")
writeBed(downHEK293, "subSets/downHEK293.bed", "downHEK293")
writeBed(upNPC, "subSets/upNPC.bed", "upNPC")
writeBed(downNPC, "subSets/downNPC.bed", "downNPC")
writeBed(upSKNSH, "subSets/upSKNSH.bed", "upSKNSH")
writeBed(downSKNSH, "subSets/downSKNSH.bed", "downSKNSH")
writeBed(upHEPG2, "subSets/upHEPG2.bed", "upHEPG2")
writeBed(downHEPG2, "subSets/downHEPG2.bed", "downHEPG2")

