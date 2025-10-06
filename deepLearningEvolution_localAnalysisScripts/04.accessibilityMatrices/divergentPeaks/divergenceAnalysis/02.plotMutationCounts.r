library(ggpubr)
library(ggridges)
library(dplyr)

processAndMergeData <- function(mutationCountFile, allRegionsFile, condition) {
  mutationCountData <- read.csv(mutationCountFile, header=TRUE, sep="\t")
  allData <- read.csv(allRegionsFile, header=FALSE, sep="\t")
  allData <- select(allData, c("V1", "V2", "V3"))
  colnames(allData) <- c("chr", "start", "end")
  allData$name <- paste(allData$chr, allData$start, allData$end, sep=".")
  mergedData <- merge(allData, mutationCountData[, c("name", "mutation_count")], by="name", all.x=TRUE)
  mergedData$mutation_count[is.na(mergedData$mutation_count)] <- 0
  mergedData$length <- mergedData$end - mergedData$start
  mergedData$divergence <- mergedData$mutation_count / mergedData$length * 100 # mutations per 100 bases
  mergedData$condition <- condition
  return(mergedData)
}

upInHomininData <- processAndMergeData("mutationCountsSummary/upInHominins.merged.Summary.txt", "../upInHominins.merged.bed", "upInHominins")
upInGreatApesData <- processAndMergeData("mutationCountsSummary/upInGreatApes.merged.Summary.txt", "../upInGreatApes.merged.bed", "upInGreatApes")
noChangeData <- processAndMergeData("mutationCountsSummary/noChangeHomininGreatApe.merged.Summary.txt", "../noChange.merged.bed", "noChange")
noChangeSizeMatchedData <- processAndMergeData("mutationCountsSummary/noChange.subsampled.merged.Summary.txt", "../noChange.subsampled.merged.bed", "noChangeSizeMatched")

mergedData <- rbind(upInHomininData, upInGreatApesData, noChangeData, noChangeSizeMatchedData)
mergedData$condition <- factor(mergedData$condition, levels = c("upInHominins", "upInGreatApes", "noChange", "noChangeSizeMatched"))

custom_colors <- c("upInHominins" = "#E96D4E",  
                   "upInGreatApes" = "#254857", 
                   "noChange" = "#6D6E71",     
                   "noChangeSizeMatched" = "#000000")

genomeWideAverageDivergence <- 100 * (18661171 / 3088269832) # number of mutations in human/HCA vcf file divided by total length of hg38, scaled per 100 base pairs

pdf("divergenceAnalysis.pdf", height=3, width=3)
ggplot(mergedData, aes(x = divergence, color = condition)) +
  stat_ecdf(geom = "step", size=1) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Divergence", y = "Cumulative Probability", color = "Condition") +
  theme_classic() + xlim(0, 3) + theme(legend.position="none") +
  labs(x="Mutations Per 100 Base Pairs", title="Greater Sequence Divergence in Functionally Divergent Elements") + geom_vline(xintercept=genomeWideAverageDivergence, linetype="dashed")
dev.off()

pairwise.t.test(mergedData$divergence, mergedData$condition, p.adjust.method = "bonferroni")

