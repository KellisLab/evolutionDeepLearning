library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)

files <- list.files("hopsOverlapsPValue", pattern = "\\.txt$", full.names = TRUE)

data_list <- list()
for (file in files) {
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("Chr", "Start", "End", "Pn_pValue")
  filename <- tools::file_path_sans_ext(basename(file))
  df$Condition <- filename
  data_list[[filename]] <- df
}

mergedData <- do.call(rbind, data_list)
mergedData$minusLogP <- -log10(mergedData$Pn_pValue)

sig_threshold <- 3
# use Wald's method to estimate SE
prop_data <- mergedData %>%
  mutate(Significant = minusLogP > sig_threshold) %>%
  group_by(Condition) %>%
  summarise(
    Count_Significant = sum(Significant),
    Count_Not_Significant = sum(!Significant),
    Total = n(),
    PropSig = Count_Significant / Total,
    CI_lower = PropSig - 1.96 * sqrt((PropSig * (1 - PropSig)) / Total),
    CI_upper = PropSig + 1.96 * sqrt((PropSig * (1 - PropSig)) / Total),
    .groups = "drop"
  ) %>%
  mutate(
    # Split Condition string and extract State (2nd last) and Group (4th last)
    State = sapply(strsplit(Condition, "\\."), function(x) ifelse(length(x) >= 2, x[length(x)-1], NA)),
    Group = sapply(strsplit(Condition, "\\."), function(x) ifelse(length(x) >= 4, x[length(x)-3], NA))
  )

tssData <- subset(prop_data, State=="TssA")

a <- ggplot(tssData, aes(x = Condition, y = PropSig, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = rev(c("#E96D4E", "#254857", "#6D6E71", "#000000"))) +
  theme_classic() +
  labs(y = "Proportion of Overlapping Variants with Significant Horizonal Pleiotropy", x = "", title="TssA") +
  theme(legend.position = "none") + coord_flip()
pairwise_fisher_test <- function(row1, row2) {
  test_matrix <- matrix(
    c(row1$Count_Significant, row1$Count_Not_Significant,
      row2$Count_Significant, row2$Count_Not_Significant),
    nrow = 2, byrow = TRUE
  )
  test_result <- fisher.test(test_matrix)
  return(test_result$p.value)
}
pairwise_results <- expand.grid(Condition1 = tssData$Condition, 
                                Condition2 = tssData$Condition, 
                                stringsAsFactors = FALSE) %>%
  filter(Condition1 != Condition2) %>%  # Avoid self-comparisons
  arrange(Condition1, Condition2) %>%
  distinct() %>% 
  rowwise() %>%
  mutate(p_value = {
    row1 <- prop_data %>% filter(Condition == Condition1)
    row2 <- prop_data %>% filter(Condition == Condition2)
    if (nrow(row1) > 0 & nrow(row2) > 0) {
      pairwise_fisher_test(row1, row2)
    } else {
      NA
    }
  })
pairwise_results$pAdj <- p.adjust(pairwise_results$p_value, method="fdr")


otherRegData <- subset(prop_data, State=="OtherRegulatoryStates")
b <- ggplot(otherRegData, aes(x = Condition, y = PropSig, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = rev(c("#E96D4E", "#254857", "#6D6E71", "#000000"))) +
  theme_classic() +
  labs(y = "Proportion of Overlapping Variants with Significant Horizonal Pleiotropy", x = "", title="Other Regulatory States") +
  theme(legend.position = "none") + coord_flip()

otherReg_pairwise_results <- expand.grid(Condition1 = otherRegData$Condition, 
                                Condition2 = otherRegData$Condition, 
                                stringsAsFactors = FALSE) %>%
  filter(Condition1 != Condition2) %>%  # Avoid self-comparisons
  arrange(Condition1, Condition2) %>%
  distinct() %>% 
  rowwise() %>%
  mutate(p_value = {
    row1 <- prop_data %>% filter(Condition == Condition1)
    row2 <- prop_data %>% filter(Condition == Condition2)
    if (nrow(row1) > 0 & nrow(row2) > 0) {
      pairwise_fisher_test(row1, row2)
    } else {
      NA
    }
  })
otherReg_pairwise_results$pAdj <- p.adjust(otherReg_pairwise_results$p_value, method="fdr")

pdf("pdf/hopsWithLabels.pdf")
grid.arrange(a,b)
dev.off()

aNoLabel <- a + theme(axis.text.y=element_blank())
bNoLabel <- b + theme(axis.text.y=element_blank())
pdf("pdf/hopsWithoutLabels.pdf", height=8, width=5)
grid.arrange(aNoLabel, bNoLabel)
dev.off()
