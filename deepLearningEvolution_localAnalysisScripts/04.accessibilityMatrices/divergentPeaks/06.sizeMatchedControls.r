library(ggpubr)
library(dplyr)
library(tidyr)

processBeds <- function(bedFile, condition) {
  data <- read.csv(bedFile, header=FALSE, sep="\t")
  colnames(data) <- c("chr", "start", "end", "log2Fc", "score", "strand", "minusLogPAdj")
  data$name <- paste(data$chr, data$start, data$end, sep=".")
  data$length <- data$end - data$start
  data$Change <- condition
  return(data)
}

upInHomininData <- processBeds("upInHominins.merged.bed", "upInHominins")
upInGreatApesData <- processBeds("upInGreatApes.merged.bed", "upInGreatApes")
noChangeData <- processBeds("noChange.merged.bed", "noChange")

mergedData <- rbind(noChangeData, upInHomininData, upInGreatApesData)
table(mergedData$Change)
gghistogram(mergedData,  "length", y="..density..", facet.by="Change", fill="Change", bins=50) + xlim(0, 4000)

medians <- mergedData %>%
  group_by(Change) %>%
  summarize(median_length = median(length, na.rm = TRUE))

pdf("pdf/lengthDistributions.pdf", height=4, width=4)
ggplot(mergedData, aes(x = length, color = Change)) +
  stat_ecdf(geom = "step", size = 1) +  # Draws CDF curves
  xlim(400, 4000) +
  labs(title = "Shorter Open Chromatin Elements Are More Likely to Exhibit Divergent Accessibility", 
       x = "Element Length (bp)", 
       y = "Cumulative Probability") +
  theme_classic() + 
  scale_color_manual(values = c("#6D6E71", "#254857", "#E96D4E")) + 
  theme(legend.position="none")
dev.off()

pairwise_results <- pairwise.t.test(mergedData$length, mergedData$Change,
                                    p.adjust.method = "bonferroni")

divergentRegions <- rbind(upInHomininData, upInGreatApesData)
divMin <- min(divergentRegions$length, na.rm=TRUE)
divMax <- max(divergentRegions$length, na.rm=TRUE)
noChangeDataFiltered <- noChangeData %>%
  filter(length >= divMin & length <= divMax)
binBreaks <- quantile(divergentRegions$length, 
                      probs = seq(0, 1, length.out = 50),
                      na.rm = TRUE)
binBreaks <- unique(binBreaks)  # remove duplicates if any
noChangeDataFiltered <- noChangeDataFiltered %>%
  mutate(bin = cut(length, breaks = binBreaks, include.lowest = TRUE))
divergentRegions <- divergentRegions %>%
  mutate(bin = cut(length, breaks = binBreaks, include.lowest = TRUE))
binCountsDivergent <- divergentRegions %>%
  group_by(bin) %>%
  summarise(n_div = n(), .groups = "drop")
binCountsNoChange <- noChangeDataFiltered %>%
  group_by(bin) %>%
  summarise(n_noChange = n(), .groups = "drop")

binSampleSizes <- binCountsDivergent %>%
  left_join(binCountsNoChange, by = "bin") %>%
  mutate(
    n_noChange = ifelse(is.na(n_noChange), 0, n_noChange),
    sampleSize = pmin(n_div, n_noChange)  # example: 1-to-1 match
  ) %>%
  select(bin, sampleSize)


set.seed(123) 
noChangeSubsample <- data.frame()

for (i in seq_len(nrow(binSampleSizes))) {
  thisBin <- binSampleSizes$bin[i]
  thisN   <- binSampleSizes$sampleSize[i]
  df_candidates <- noChangeDataFiltered[noChangeDataFiltered$bin == thisBin, ]
  toSample <- min(nrow(df_candidates), thisN)
  sampledDF <- df_candidates %>%
    dplyr::slice_sample(n = toSample)
  noChangeSubsample <- rbind(noChangeSubsample, sampledDF)
}
noChangeSubsample <- noChangeSubsample %>% select(-bin)
noChangeSubsample$Change <- "noChangeSubsample"


mergedData <- rbind(mergedData, noChangeSubsample)
gghistogram(mergedData,  "length", y="density", facet.by="Change", fill="Change") + xlim(0, 4000)

write_bed <- function(df, filename) {
  bed <- df %>%
    mutate(
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    select(chr, start, end, name)
  write.table(bed, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

write_bed(noChangeSubsample, "noChange.subsampled.merged.bed")

