# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Directory containing the BED files
bed_dir <- "byBiosample"

# List all BED files in the directory
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

# Initialize an empty list to store data
peak_data <- list()

# Loop through each BED file
for (file in bed_files) {
  # Extract biosample name from the file name (remove directory and .bed extension)
  biosample <- gsub("\\.bed$", "", basename(file))
  
  # Read the BED file
  bed <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract peak lengths (end - start)
  bed$length <- bed$V3 - bed$V2
  
  # Store the data in the list
  peak_data[[biosample]] <- data.frame(
    Biosample = biosample,
    Length = bed$length
  )
}

# Combine all data into a single data frame
all_peaks <- do.call(rbind, peak_data)
median_peak_length <- median(all_peaks$Length, na.rm = TRUE)

hist_data <- ggplot_build(ggplot(all_peaks, aes(x = Length)) +
                            geom_histogram(binwidth = 50))$data[[1]]
max_count <- max(hist_data$count)

hist_peak_length <- ggplot(all_peaks, aes(x = Length)) +
  geom_histogram(binwidth = 50, fill = "#f4a261") +
  geom_vline(xintercept = median_peak_length, color = "black", linetype = "dashed") +
  annotate("text", x = median_peak_length + 200, y = max(max_count * 0.9), label = paste("Median:", median_peak_length), color = "black", size = 5, hjust = 0) +
  theme_classic() +
  xlab("Peak Length (bp)") +
  ylab("Count") +
  ggtitle("Distribution of Peak Lengths Across All Biosamples")
peak_counts <- all_peaks %>%
  group_by(Biosample) %>%
  summarise(PeakCount = n())
hist_peak_counts <- ggplot(peak_counts, aes(x = PeakCount)) +
  geom_histogram(binwidth = 2500, fill = "#2a9d8f") +
  theme_classic() +
  xlab("Number of Peaks per Biosample") +
  ylab("Count") +
  ggtitle("Number of Peaks per Biosample")
pdf("peakLengthAndCountHistograms.pdf", height=3, width=10)
grid.arrange(hist_peak_length, hist_peak_counts, nrow = 1)
dev.off()

##### Jaccard Matrix #####
library(pheatmap)
data <- read.csv("jaccard_result.txt", header=FALSE, sep=" ")
colnames(data) <- c("Biosample1", "Biosample2", "JaccardIndex")
jaccard_matrix <- dcast(data, Biosample1 ~ Biosample2, value.var = "JaccardIndex")
rownames(jaccard_matrix) <- jaccard_matrix$Biosample1
jaccard_matrix <- jaccard_matrix[, -1]
jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]
pheatmap(jaccard_matrix, show_colnames = FALSE, fontsize_row=5)


#### Complex strategy ###
library(ComplexHeatmap)
library(circlize)
library(dplyr)

jaccard_matrix <- as.matrix(jaccard_matrix)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- rownames(jaccard_matrix)
peak_counts <- peak_counts %>% arrange(Biosample)
rownames(peak_counts) <- peak_counts$Biosample
median_lengths <- all_peaks %>%
  group_by(Biosample) %>%
  summarise(MedianLength = median(Length)) %>%
  arrange(Biosample)

row_ha <- rowAnnotation(
  `Peak Count` = anno_barplot(
    peak_counts$PeakCount,
    border = FALSE,
    gp = gpar(fill = "steelblue"),
    width = unit(3, "cm")
  )
)

pdf("testHeatmap.pdf", height=7, width=10)
Heatmap(
  jaccard_matrix,
  name = "Jaccard Index",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Jaccard",
    at = seq(0, 1, by = 0.2),
    labels = seq(0, 1, by = 0.2)
  ),
  right_annotation = row_ha,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5)
)
dev.off()




#### Enrichment ####
data <- read.csv("enrichmentResults.txt", header=FALSE, sep=" ")
colnames(data) <- c("Set1", "Set2", "Enrichment")
data$log2Enrich <- log2(data$Enrichment)
data$log2Enrich[!is.finite(data$log2Enrich)] <- 0 
data <- select(data, -c("Enrichment"))

wide_data <- data %>%
  pivot_wider(names_from="Set1", values_from="log2Enrich")
wide_data[is.na(wide_data)] <- 0

wide_data <- as.data.frame(wide_data)
rownames(wide_data) <- wide_data$Set2
wide_data <- select(wide_data, -c("Set2"))

pheatmap(wide_data)
