library(ggpubr)


data <- read.csv("testIdr.txt", header=FALSE, sep="\t")
colnames(data) <- c(
  "chrom", "chromStart", "chromEnd", "name", "score", "strand",
  "signalValue", "pValue", "qValue",
  "rep1_chromStart", "rep1_chromEnd", "OverlapLength",
  "rep1_signalValue", "globalIDR"
)
head(data)

data$Rank1 <- rank(-data$Signal1)
data$Rank2 <- rank(-data$Signal2)

