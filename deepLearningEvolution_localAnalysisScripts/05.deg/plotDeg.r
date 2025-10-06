library(ggpubr)
library(dplyr)
library(tidyr)

data <- read.csv("summary.overlaps.txt", header=FALSE, sep="\t")

data <- data %>% select(c("V2", "V3", "V9", "V10", "V11"))

colnames(data) <- c("Filename1", "Filename2", "Enrich", "pEnrich", "pDeplete")
data$log2Enrich <- log2(data$Enrich)
parsed <- do.call(rbind, strsplit(as.character(data$Filename1), split = "\\."))
data_parsed <- data.frame(
  pCutoff   = parsed[, 2],
  direction = parsed[, 4],
  fcCutoff  = parsed[, 5],
  padLength = parsed[, 6],
  stringsAsFactors = FALSE
)
data_parsed$pCutoff   <- as.numeric(data_parsed$pCutoff)
data_parsed$fcCutoff  <- as.numeric(data_parsed$fcCutoff)
data_parsed$padLength <- as.integer(data_parsed$padLength)
head(data_parsed)

final_data <- cbind(data, data_parsed) %>% subset(pCutoff == 3 & padLength == 50000 & fcCutoff == 2)
final_data$pAdj <- p.adjust(pmin(final_data$pEnrich, final_data$pDeplete))
final_data$minusLogPAdj <- -log10(final_data$pAdj)
ggscatter(final_data, "log2Enrich", "minusLogPAdj", color="Filename2", label="Filename1")

ggscatter(final_data, "padLength", "minusLogPAdj", color="log2Enrich", facet.by="direction")

pdf("degBalloonPlot.pdf", height=3, width=4)
ggballoonplot(final_data, x="Filename2", y="direction", fill="log2Enrich", size="minusLogPAdj") + 
  theme_classic() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0) +
  labs(
    x="",
    y="Human/Great Ape DEG Direction"
  )
dev.off()
