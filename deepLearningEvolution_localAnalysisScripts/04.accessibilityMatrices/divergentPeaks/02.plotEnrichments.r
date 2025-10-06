library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggrepel)

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
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +  # Custom colors
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
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +  # Custom colors
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
  geom_point(alpha = 1, size = 3) +
  geom_text_repel(
    max.overlaps = 10,
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +  # Custom colors
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

pdf("pdf/roiEnrichmentPlots.pdf", height=10, width=7)
grid.arrange(haqerPanel, harPanel, hgePanel)
dev.off()

haqerBox <- ggboxplot(subset(data, Roi == "haqer.ordered"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) + 
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "HAQER"
  )

harBox <- ggboxplot(subset(data, Roi == "HAR_Walsh_List_hg38"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) + 
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "HAR"
  )

hgeBox <- ggboxplot(subset(data, Roi == "HGE_merged"), "Status", "log2FC", add="jitter", color="Status", fill="Status", alpha=0.5) + 
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) + 
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(
    x="Predicted Accessibility Change (Enformer)",
    y="log2(Overlap Enrichment)",
    title = "Human Gained Enhancer/Promoter"
  )

hge_data <- subset(data, Roi == "HGE_merged" & is.finite(log2FC))
wilcox.test(log2FC ~ Status, data = hge_data)

pdf("roiBoxPlots.pdf", height=8, width=5)
grid.arrange(haqerBox, harBox, hgeBox)
dev.off()

pdf("pdf/hgePanels.pdf", height=5, width=10)
grid.arrange(hgePanel, hgeBox, nrow=1)
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
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e"))

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
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e"))

pdf("pdf/hCondelValidation.pdf", height=4, width=8)
grid.arrange(a, b, nrow=1)
dev.off()

## labeled hCONDEL

library(ggrepel)

aLabeled <- ggboxplot(
  subset(data, Roi == "allDown"), 
  "Status", "minusLog10pAdj", 
  fill = "Status", color = "Status", alpha = 0.5, add = "jitter"
) +
  theme(legend.position = "none") + 
  ylab("Enrichment Significance -log10(pAdj)") + 
  ggtitle("hCONDELs with human loss of MPRA activity") + 
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  geom_text_repel(
    data = subset(data, Roi == "allDown" & minusLog10pAdj > 3), 
    aes(x = Status, y = minusLog10pAdj, label = Biosample),
    size = 3, box.padding = 0.5, max.overlaps = 10
  )

bLabeled <- ggboxplot(
  subset(data, Roi == "allUp"), 
  "Status", "minusLog10pAdj", 
  fill = "Status", color = "Status", alpha = 0.5, add = "jitter"
) +
  theme(legend.position = "none") + 
  xlab("Predicted Accessibility Change (Enformer)") + 
  ylab("-log10(pAdj)") + 
  ggtitle("hCONDELs with human gain of MPRA activity") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
  scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
  geom_text_repel(
    data = subset(data, Roi == "allUp" & minusLog10pAdj > 3), 
    aes(x = Status, y = minusLog10pAdj, label = Biosample),
    size = 3, box.padding = 0.5, max.overlaps = 10
  )


pdf("pdf/hCondelValidationLabeled.pdf", height=4, width=8)
grid.arrange(aLabeled, bLabeled, nrow=1)
dev.off()



#### HGE Extended Figure ###
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)

makeVolcano <- function(data) {
  plot_title <- deparse(substitute(data))
  p <- ggplot(data, aes(x = log2FC, y = minusLog10pAdj, color = Status, label = Biosample)) +
    geom_point(alpha = 1, size = 3) +
    geom_text_repel(
      max.overlaps = 10,
      size = 3,
      box.padding = 0.3,
      point.padding = 0.3
    ) +
    scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +  # Custom colors
    labs(
      x = NULL,  # Remove x-axis label for alignment
      y = "Enrichment Significance [-log10(pAdj)]",
      title = plot_title
    ) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  return(p)
}

makeBox <- function(data) {
  plot_title <- deparse(substitute(data))
  ggboxplot(data, "Status", "log2FC", add = "jitter", color = "Status", fill = "Status", alpha = 0.5) + 
    scale_color_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) + 
    scale_fill_manual(values = c("downregulated" = "#254857", "upregulated" = "#e96d4e")) +
    scale_x_discrete(labels = c("downregulated" = "Human Loss", "upregulated" = "Human Gain")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),  # Remove y-axis label for alignment
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_blank()     # Remove title to avoid redundancy
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = NULL,
      y = NULL
    ) + 
    coord_flip()
}

Frontal_12Weeks_Enhancer_data <- subset(data, Roi == "GSE63648_12Fpcw_ac_Hu_gain.hg38")
Frontal_8point5Weeks_Enhancer_data <- subset(data, Roi == "GSE63648_8_5pcw_ac_Hu_gain.hg38")
Frontal_7Weeks_Enhancer_data <- subset(data, Roi == "GSE63648_7pcw_ac_Hu_gain.hg38")
Occipital_12Weeks_Enhancer_data <- subset(data, Roi == "GSE63648_12Opcw_ac_Hu_gain.hg38")

Frontal_12Weeks_Promoter_data <- subset(data, Roi == "GSE63648_12Fpcw_me2_Hu_gain.hg38")
Frontal_8point5Weeks_Promoter_data <- subset(data, Roi == "GSE63648_8_5pcw_me2_Hu_gain.hg38")
Frontal_7Weeks_Promoter_data <- subset(data, Roi == "GSE63648_7pcw_me2_Hu_gain.hg38")
Occipital_12Weeks_Promoter_data <- subset(data, Roi == "GSE63648_12Opcw_me2_Hu_gain.hg38")


Frontal_12Weeks_Enhancer_volcano <- makeVolcano(Frontal_12Weeks_Enhancer_data)
Frontal_12Weeks_Enhancer_box <- makeBox(Frontal_12Weeks_Enhancer_data)

# Generate all volcano and box plots
Frontal_12Weeks_Enhancer_volcano <- makeVolcano(Frontal_12Weeks_Enhancer_data)
Frontal_12Weeks_Enhancer_box <- makeBox(Frontal_12Weeks_Enhancer_data)

Frontal_8point5Weeks_Enhancer_volcano <- makeVolcano(Frontal_8point5Weeks_Enhancer_data)
Frontal_8point5Weeks_Enhancer_box <- makeBox(Frontal_8point5Weeks_Enhancer_data)

Frontal_7Weeks_Enhancer_volcano <- makeVolcano(Frontal_7Weeks_Enhancer_data)
Frontal_7Weeks_Enhancer_box <- makeBox(Frontal_7Weeks_Enhancer_data)

Occipital_12Weeks_Enhancer_volcano <- makeVolcano(Occipital_12Weeks_Enhancer_data)
Occipital_12Weeks_Enhancer_box <- makeBox(Occipital_12Weeks_Enhancer_data)

Frontal_12Weeks_Promoter_volcano <- makeVolcano(Frontal_12Weeks_Promoter_data)
Frontal_12Weeks_Promoter_box <- makeBox(Frontal_12Weeks_Promoter_data)

Frontal_8point5Weeks_Promoter_volcano <- makeVolcano(Frontal_8point5Weeks_Promoter_data)
Frontal_8point5Weeks_Promoter_box <- makeBox(Frontal_8point5Weeks_Promoter_data)

Frontal_7Weeks_Promoter_volcano <- makeVolcano(Frontal_7Weeks_Promoter_data)
Frontal_7Weeks_Promoter_box <- makeBox(Frontal_7Weeks_Promoter_data)

Occipital_12Weeks_Promoter_volcano <- makeVolcano(Occipital_12Weeks_Promoter_data)
Occipital_12Weeks_Promoter_box <- makeBox(Occipital_12Weeks_Promoter_data)

Frontal_12Weeks_Enhancer_pair <- plot_grid(
  Frontal_12Weeks_Enhancer_volcano, Frontal_12Weeks_Enhancer_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Frontal_8point5Weeks_Enhancer_pair <- plot_grid(
  Frontal_8point5Weeks_Enhancer_volcano, Frontal_8point5Weeks_Enhancer_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Frontal_7Weeks_Enhancer_pair <- plot_grid(
  Frontal_7Weeks_Enhancer_volcano, Frontal_7Weeks_Enhancer_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Occipital_12Weeks_Enhancer_pair <- plot_grid(
  Occipital_12Weeks_Enhancer_volcano, Occipital_12Weeks_Enhancer_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Frontal_12Weeks_Promoter_pair <- plot_grid(
  Frontal_12Weeks_Promoter_volcano, Frontal_12Weeks_Promoter_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Frontal_8point5Weeks_Promoter_pair <- plot_grid(
  Frontal_8point5Weeks_Promoter_volcano, Frontal_8point5Weeks_Promoter_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Frontal_7Weeks_Promoter_pair <- plot_grid(
  Frontal_7Weeks_Promoter_volcano, Frontal_7Weeks_Promoter_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

Occipital_12Weeks_Promoter_pair <- plot_grid(
  Occipital_12Weeks_Promoter_volcano, Occipital_12Weeks_Promoter_box,
  ncol = 1, align = "v", rel_heights = c(2, 1)
)

final_grid <- plot_grid(
  Frontal_7Weeks_Enhancer_pair, Frontal_8point5Weeks_Enhancer_pair, Frontal_12Weeks_Enhancer_pair, Occipital_12Weeks_Enhancer_pair, 
  Frontal_7Weeks_Promoter_pair, Frontal_8point5Weeks_Promoter_pair, Frontal_12Weeks_Promoter_pair, Occipital_12Weeks_Promoter_pair, 
  ncol = 4,
  align = "hv"
)


pdf("pdf/hgeExtendedGrid.pdf", height=12, width=15)
final_grid
dev.off()

# perform T tests. Will add results to figure manually in illustrator
perform_t_test <- function(data) {
  data_filtered <- subset(data, is.finite(log2FC))  # Exclude -Inf and Inf values
  t_test_result <- t.test(log2FC ~ Status, data = data_filtered)
  return(t_test_result)
}

# Function to perform Wilcoxon rank-sum test (including -Inf values)
perform_wilcoxon_test <- function(data) {
  wilcox_test_result <- wilcox.test(log2FC ~ Status, data = data, exact = FALSE)  # exact=FALSE for large samples
  return(wilcox_test_result)
}

perform_wilcoxon_test(Frontal_12Weeks_Enhancer_data)
perform_wilcoxon_test(Frontal_8point5Weeks_Enhancer_data)
perform_wilcoxon_test(Frontal_7Weeks_Enhancer_data)
perform_wilcoxon_test(Occipital_12Weeks_Enhancer_data)

perform_wilcoxon_test(Frontal_12Weeks_Promoter_data)
perform_wilcoxon_test(Frontal_8point5Weeks_Promoter_data)
perform_wilcoxon_test(Frontal_7Weeks_Promoter_data)
perform_wilcoxon_test(Occipital_12Weeks_Promoter_data)








#### Plot HAQER/HAR/UNICORN Enrichments
library(dplyr)
library(ggpubr)
library(gridExtra)
data <- read.csv("enrichment.HaqerHarUnicorn.summary.txt", header=FALSE, sep="\t")
data <- select(data, c("V2", "V3", "V9", "V10", "V11"))
colnames(data) <- c("File1", "File2", "Enrich", "pEnrich", "pDeplete")
data$rawP <- pmin(1, 2 * pmin(data$pEnrich, data$pDeplete))
data$pAdj <- p.adjust(data$rawP, method = "fdr")
data$minusLog10pAdj <- -log10(data$pAdj)
data$log2FC <- log2(data$Enrich)
data$Biosample <- sub(".*/(.*)\\..*\\.bed", "\\1", data$File1)
data$Status <- sub(".*/.*\\.(.*)\\.bed", "\\1", data$File1)
data$Roi <- sub(".*/(.*)\\.bed", "\\1", data$File2)

ggscatter(data, x="log2FC", y="minusLog10pAdj", color="Roi", facet.by = "Status")
ggboxplot(data, x="Roi", y="minusLog10pAdj", color="Roi", fill="Roi", alpha=0.5, facet.by="Status", add="jitter")


data <- read.csv("unicornEnrichmentSummary.txt", header=FALSE, sep="\t")
data <- select(data, c("V2", "V3", "V9", "V10", "V11"))
colnames(data) <- c("File1", "File2", "Enrich", "pEnrich", "pDeplete")
data$rawP <- pmin(1, 2 * pmin(data$pEnrich, data$pDeplete))
data$pAdj <- p.adjust(data$rawP, method = "fdr")
data$minusLog10pAdj <- -log10(data$pAdj)
data$log2FC <- log2(data$Enrich)
data$Biosample <- sub(".*/(.*)\\..*\\.bed", "\\1", data$File1)
data$Status <- sub(".*/.*\\.(.*)\\.bed", "\\1", data$File1)
data$Roi <- sub(".*/(.*)\\.bed", "\\1", data$File2)

pdf("pdf/unicornRoiEnrichment.pdf", height=5, width=7)
ggplot(data, aes(
  x = Status, 
  y = Roi, 
  size = minusLog10pAdj, 
  fill = log2FC
)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +  # Bubbles with outline
  scale_size_continuous(name = "-log10(FDR-adjusted p)", range = c(2, 10)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, name = "log2(Enrichment)"
  ) +
  labs(
    x = "Accessibility Status", 
    y = "Region of Interest (ROI)", 
    title = "Enrichment for Regions of Evolutionary Interest"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  )
dev.off()

pdf("pdf/unlabeled.unicornRoiEnrichment.pdf", height=1.5, width=1.5)
ggplot(data, aes(
  x = Status, 
  y = Roi, 
  size = minusLog10pAdj, 
  fill = log2FC
)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +  # Bubbles with outline
  scale_size_continuous(name = "-log10(FDR-adjusted p)", range = c(2, 10)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, name = "log2(Enrichment)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
  ) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  )
dev.off()
