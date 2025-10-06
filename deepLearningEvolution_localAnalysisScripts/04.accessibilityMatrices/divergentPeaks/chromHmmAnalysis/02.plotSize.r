library(dplyr)
library(ggpubr)
library(gridExtra)
library(tidyr)


readBedFile <- function(filename) {
  bed <- read.csv(filename, header=FALSE, sep="\t")
  bed <- bed %>% select(c("V1", "V2", "V3"))
  colnames(bed) <- c("chr", "start", "end")
  bed$length <- bed$end - bed$start
  bed$filename <- basename(filename)
  fields <- unlist(strsplit(basename(filename), "\\."))
  bed$group <- fields[1]
  bed$chromHmmCategory <- fields[3]
  return(bed)
}

bed_files <- list.files("enhancerPromoterSplits", pattern="\\.bed$", full.names=TRUE)
mergedData <- do.call(rbind, lapply(bed_files, readBedFile))


line_counts_df <- mergedData %>%
  group_by(filename, group, chromHmmCategory) %>%
  summarise(line_count = n(), .groups = "drop")

proportions_df <- mergedData %>%
  group_by(group, chromHmmCategory) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(TssA_status = ifelse(chromHmmCategory == "TssA", "TssA", "nonTssA")) %>%
  mutate(regulatory_status = ifelse(chromHmmCategory == "nonRegulatory", "nonRegulatory", "regulatory"))

# Fisher's on TssA vs. nonTssA
contingency_data <- proportions_df %>%
  group_by(group, TssA_status) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = TssA_status, values_from = total, values_fill = 0)
run_fishers_test <- function(df, group1, group2) {
  subset_df <- df %>% filter(group %in% c(group1, group2))
  contingency_table <- as.matrix(subset_df[, c("TssA", "nonTssA")])
  rownames(contingency_table) <- subset_df$group
  fisher_result <- fisher.test(contingency_table)
  return(data.frame(
    Comparison = paste(group1, "vs", group2),
    p_value = fisher_result$p.value
  ))
}
comparisons <- list(
  c("noChangeHomininGreatApe", "upInHominins"),
  c("noChangeHomininGreatApe", "upInGreatApes"),
  c("upInHominins", "upInGreatApes")
)
fishers_results <- do.call(rbind, lapply(comparisons, function(x) run_fishers_test(contingency_data, x[1], x[2])))
print(fishers_results)

custom_colors <- c("#E96D4E", "#254857", "#6D6E71")

tssAPlotData <- subset(proportions_df, chromHmmCategory == "TssA")
tssAPlotData$group <- factor(tssAPlotData$group, levels=c("upInHominins", "upInGreatApes", "noChangeHomininGreatApe"))

# Fisher's noRegualtoryElement vs. regulatoryElement
regContingencyTable <- proportions_df %>%
  group_by(group, regulatory_status) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = regulatory_status, values_from = total, values_fill = 0)
run_fishers_test_reg <- function(df, group1, group2) {
  subset_df <- df %>% filter(group %in% c(group1, group2))
  contingency_table <- as.matrix(subset_df[, c("regulatory", "nonRegulatory")])
  rownames(contingency_table) <- subset_df$group
  fisher_result <- fisher.test(contingency_table)
  return(data.frame(
    Comparison = paste(group1, "vs", group2),
    p_value = fisher_result$p.value
  ))
}
regFishersResults <- do.call(rbind, lapply(comparisons, function(x) run_fishers_test_reg(regContingencyTable, x[1], x[2])))
nonRegulatoryPlotData <- subset(proportions_df, chromHmmCategory == "nonRegulatory")
nonRegulatoryPlotData$group <- factor(nonRegulatoryPlotData$group, levels=c("upInHominins", "upInGreatApes", "noChangeHomininGreatApe"))


a <- ggbarplot(tssAPlotData, x="group", y="proportion", fill="group", color="group", palette= custom_colors) +
  theme(axis.title.x=element_blank(), legend.position="none") +
  labs(y="Proportion of Elements", title="Proportion of Elements Overlapping Active TSS Regions")
b <- ggbarplot(nonRegulatoryPlotData, x="group", y="proportion", fill="group", color="group", palette= custom_colors) +
  theme(axis.title.x=element_blank(), legend.position="none") +
  labs(y="Proportion of Elements", title="Proportion of Elements Not Overlapping ChromHmm Regulatory States")

pdf("pdf/EnhancerPromoter.pdf", height=8, width=5)
grid.arrange(a,b)
dev.off()

### Length Distributions ###

ggplot(mergedData, aes(x=filename, y=length)) +
  geom_boxplot(outlier.shape=NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate labels for readability
  labs(x="Filename", y="Region Length", title="Distribution of Region Lengths by File") +
  ylim(0, 3000) + coord_flip()

pdf("pdf/TssA.lengthDistributions.pdf", height=6, width=6)
ggplot(subset(mergedData, chromHmmCategory == "TssA"), aes(x=length, color=group)) +
  stat_ecdf(geom="step", size=1) +
  scale_color_manual(values=rev(custom_colors)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(x="Region Length", y="Cumulative Probability",
       title="Cumulative Distribution of Region Lengths by Group (TssA)") +
  scale_x_continuous(limits=c(0, 3000))
dev.off()

tssA_data <- subset(mergedData, chromHmmCategory == "TssA")
pairwise_results <- pairwise.t.test(tssA_data$length, tssA_data$group,
                                    p.adjust.method = "bonferroni")
otherReg_data <- subset(mergedData, chromHmmCategory == "OtherRegulatoryStates")
pairwise_results <- pairwise.t.test(otherReg_data$length, otherReg_data$group,
                                    p.adjust.method = "bonferroni")

pdf("pdf/chromHmmSizeDifference.pdf", height=6, width=7)
ggplot(mergedData, aes(x=length, color=chromHmmCategory)) +
  stat_ecdf(geom="step", size=1) +
  scale_color_manual(values=c("#2a9d8f", "#e76f51", "#f4a261")) +
  theme_classic() +
  labs(x="Region Length", y="Cumulative Probability",
       title="Cumulative Distribution of Region Lengths by Group (TssA)") +
  scale_x_continuous(limits=c(0, 3000))
dev.off()

pairwise_results <- pairwise.t.test(mergedData$length, mergedData$chromHmmCategory,
                                    p.adjust.method = "bonferroni")


pdf("pdf/OtherRegulatoryStates.lengthDistributions.pdf", height=6, width=6)
ggplot(subset(mergedData, chromHmmCategory == "OtherRegulatoryStates"), aes(x=length, color=group)) +
  stat_ecdf(geom="step", size=1) +
  scale_color_manual(values=rev(custom_colors)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(x="Region Length", y="Cumulative Probability",
       title="Cumulative Distribution of Region Lengths by Group (Other Regulatory States)") +
  scale_x_continuous(limits=c(0, 3000))
dev.off()
