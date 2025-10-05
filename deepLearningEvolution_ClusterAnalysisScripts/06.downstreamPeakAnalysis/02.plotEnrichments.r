library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)

data <- read.csv("enrichments.Summary.txt", header=FALSE, sep="\t")
head(data)
data <- select(data, c("V2", "V3", "V9", "V10", "V11"))
head(data)
colnames(data) <- c("File1", "File2", "Enrich", "pEnrich", "pDeplete")
head(data)

data$rawP <- pmin(1, 2 * pmin(data$pEnrich, data$pDeplete))
head(data)
data$pAdj <- p.adjust(data$rawP, method = "fdr")
head(data)

data$minusLog10pAdj <- -log10(data$pAdj)
data$log2FC <- log2(data$Enrich)
data$Biosample <- sub(".*/(.*)\\..*\\.bed", "\\1", data$File1)
data$Status <- sub(".*/.*\\.(.*)\\.bed", "\\1", data$File1)
data$Roi <- sub(".*/(.*)\\.bed", "\\1", data$File2)

haqerPanel <- ggplot(subset(data, Roi == "haqer.ordered"), aes(x = log2FC, y = minusLog10pAdj, color = Status, label = Biosample)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(
    max.overlaps = 10,
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +  # Custom colors
  labs(
    x = "log2(Overlap Enrichment)",
    y = "Enrichment Significance [-log10(pAdj)]",
    title = "HAQER"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

harPanel <- ggplot(subset(data, Roi == "HAR_Walsh_List_hg38"), aes(x = log2FC, y = minusLog10pAdj, color = Status, label = Biosample)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(
    max.overlaps = 10,
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +  # Custom colors
  labs(
    x = "log2(Overlap Enrichment)",
    y = "Enrichment Significance [-log10(pAdj)]",
    title = "HAR"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

hgePanel <- ggplot(subset(data, Roi == "HGE_merged"), aes(x = log2FC, y = minusLog10pAdj, color = Status, label = Biosample)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(
    max.overlaps = 10,
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +  # Custom colors
  labs(
    x = "log2(Overlap Enrichment)",
    y = "Enrichment Significance [-log10(pAdj)]",
    title = "Human Gained Enhancer/Promoter"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

pdf("roiEnrichmentPlots.pdf", height=10, width=7)
grid.arrange(haqerPanel, harPanel, hgePanel)
dev.off()

haqerBox <- ggboxplot(subset(data, Roi == "haqer.ordered"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) + 
  scale_fill_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "HAQER"
  )

harBox <- ggboxplot(subset(data, Roi == "HAR_Walsh_List_hg38"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) + 
  scale_fill_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "HAR"
  )

hgeBox <- ggboxplot(subset(data, Roi == "HGE_merged"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) + 
  scale_fill_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "Human Gained Enhancer/Promoter"
  )

pdf("roiBoxPlots.pdf", height=8, width=5)
grid.arrange(haqerBox, harBox, hgeBox)
dev.off()

a <- ggboxplot(
  subset(data, Roi == "allDown"), 
  "Status", "minusLog10pAdj", 
  fill = "Status", color = "Status", alpha = 0.5, add = "jitter"
) +
  theme(legend.position = "none") + 
  ylab("Enrichment Significance -log10(pAdj)") + 
  ggtitle("hCONDELs with human loss of MPRA activity") + 
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  scale_fill_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C"))

b <- ggboxplot(
  subset(data, Roi == "allUp"), 
  "Status", "minusLog10pAdj", 
  fill = "Status", color = "Status", alpha = 0.5, add = "jitter"
) +
  theme(legend.position = "none") + 
  xlab("Predicted Accessibility Change (Enformer)") + 
  ylab("-log10(pAdj)") + 
  ggtitle("hCONDELs with human gain of MPRA activity") +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  scale_fill_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C")) +
  scale_color_manual(values = c("downregulated" = "#1E3A8A", "upregulated" = "#DC143C"))

pdf("hCondelValidation.pdf", height=8, width=5)
grid.arrange(a, b)
dev.off()
