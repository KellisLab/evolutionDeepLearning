library(Matrix)
library(pheatmap)

set.seed(42)

###### IDR Peak Overview ######
library(Matrix)
library(pheatmap)
library(RColorBrewer)

data <- readRDS("biosample_binary_matrix_sparse.rds")
set.seed(42)
sampled_rows <- sample(1:nrow(data), size = 50000)
sampled_data <- data[sampled_rows, ]

metadata <- read.csv("subsampledBiosamples.txt", header=TRUE, sep="\t")
annotation_row <- metadata[, c("Group", "Lifestage", "Type")]
rownames(annotation_row) <- metadata$Biosample

# Define custom annotation colors
# Group colors: generate colors for the number of unique groups
group_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(annotation_row$Group)))

# Type colors
type_colors <- c(
  "Cell Line" = "#1c3144",
  "Tissue" = "#d00000",
  "In Vitro Differentiated" = "#ffba08"
)

# Lifestage colors
lifestage_colors <- c(
  "Cell Line" = "#219ebc",
  "Developing" = "#FFD700",
  "Adult" = "#8ecae6",
  "Newborn" = "#fb8500"
)

# Map colors to annotation levels
annotation_colors <- list(
  Group = setNames(group_colors, sort(unique(annotation_row$Group))),
  Lifestage = lifestage_colors,
  Type = type_colors
)
# Define a color scale for the heatmap
color_scale <- colorRampPalette(c("white", "black"))(50)
# Transpose the sparse matrix for plotting
dense_matrix <- as.matrix(sampled_data) # Convert to dense matrix
transposed_matrix <- t(dense_matrix) # Transpose matrix for plotting
png("mainBiosampleMatrix.png", height=3000, width=6000, res=300)  # High resolution
heatmap_obj <- pheatmap(
  mat = transposed_matrix,
  annotation_row = annotation_row, 
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,  # Perform clustering here to get the order
  cluster_cols = TRUE, 
  show_rownames = FALSE,  # Hide row names
  show_colnames = FALSE, 
  treeheight_col = 0,
  color = color_scale
)
dev.off()

ordered_samples <- rownames(transposed_matrix)[heatmap_obj$tree_row$order]
annotation_row <- annotation_row[ordered_samples, , drop = FALSE]  # Reorder to match clustering
set.seed(42)
sampled_columns <- sample(ncol(transposed_matrix), size = 500)

pdf("mainBiosampleMatrix_annotations.pdf", height=10, width=16)
pheatmap(
  mat = transposed_matrix[ordered_samples, sampled_columns, drop = FALSE],  # Keep row order same as heatmap
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # NO clustering to preserve order
  cluster_cols = FALSE,
  show_rownames = TRUE,  # Keep row names
  show_colnames = FALSE,
  treeheight_col = 0,
  color = color_scale,
  fontsize_row = 6,
  main = "Open Chromatin Elements by Biosample"
)
dev.off()

hist(
  rowSums(data),
  breaks = 100,
  xlab = "Accessibility Breadth (Number of Biosamples with Peak Accessibility)",
  ylab = "Number of Peaks",
  col = "darkgreen"
)

Li####### Peak by Haplotypes, per biosample #######
library(pheatmap)
rds_dir <- "rds/"
output_dir <- "rdsFigures/"
dir.create(output_dir, showWarnings = FALSE)

sample_groups <- c(
  # Modern Humans
  "HG00525.A" = "ModernHuman", "HG00525.B" = "ModernHuman",
  "HG02615.A" = "ModernHuman", "HG02615.B" = "ModernHuman",
  "HG03079.A" = "ModernHuman", "HG03079.B" = "ModernHuman",
  "NA18858.A" = "ModernHuman", "NA18858.B" = "ModernHuman",
  "NA20517.A" = "ModernHuman", "NA20517.B" = "ModernHuman",
  
  # Archaic Hominins
  "altai.Neanderthal.A" = "ArchaicHominin", "altai.Neanderthal.B" = "ArchaicHominin",
  "chagyrskaya.neanderthal.A" = "ArchaicHominin", "chagyrskaya.neanderthal.B" = "ArchaicHominin",
  "vindija.neanderthal.A" = "ArchaicHominin", "vindija.neanderthal.B" = "ArchaicHominin",
  "denisova.A" = "ArchaicHominin", "denisova.B" = "ArchaicHominin",
  
  # Great Apes (phylogenetic order: Chimpanzee → Bonobo → Gorilla → Orangutan)
  "clara.panTroTro.A" = "GreatApe", "clara.panTroTro.B" = "GreatApe",
  "jimmie.panTroVerus.A" = "GreatApe", "jimmie.panTroVerus.B" = "GreatApe",
  "nakuu.panTroSchweinfurthii.A" = "GreatApe", "nakuu.panTroSchweinfurthii.B" = "GreatApe",
  "taweh.panTroEllioti.A" = "GreatApe", "taweh.panTroEllioti.B" = "GreatApe",
  
  "natalie.panPaniscus.A" = "GreatApe", "natalie.panPaniscus.B" = "GreatApe",  # Bonobo
  
  "victoria.gorBeringeiGraueri.A" = "GreatApe", "victoria.gorBeringeiGraueri.B" = "GreatApe",  # Gorilla
  "dian.gorGorGor.A" = "GreatApe", "dian.gorGorGor.B" = "GreatApe",
  
  "dunja.ponAbe.A" = "GreatApe", "dunja.ponAbe.B" = "GreatApe"  # Orangutan
)

# Define custom annotation colors
annotation_colors <- list(
  Group = c(
    "ModernHuman" = "#F4A361",      # Orange
    "ArchaicHominin" = "#EAC365",   # Yellow
    "GreatApe" = "#27A092"          # Teal
  )
)

# Define a fixed order for rows (Modern → Archaic → Great Ape (phylogenetic order))
fixed_row_order <- names(sample_groups)

# Define a color scale for the heatmap
color_scale <- colorRampPalette(c("white", "darkgreen"))(50)

# Process each .rds file
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

for (file in rds_files) {
  cat("Processing:", file, "\n")
  
  # Read RDS file
  data <- readRDS(file)
  
  # Filter peaks to remove low-accessibility elements
  filtered_data <- data[, colSums(data) >= 5]
  
  # Sample a subset of columns for visualization
  sampled_cols <- sample(1:ncol(filtered_data), size = min(30000, ncol(filtered_data)))  # Avoid sampling more than available columns
  sampled_data <- filtered_data[, sampled_cols]
  
  # Ensure rows follow the fixed order (subset to existing ones in the data)
  ordered_rows <- intersect(fixed_row_order, rownames(sampled_data))
  sampled_data <- sampled_data[ordered_rows, , drop = FALSE]  # Reorder
  
  # Create annotation data frame
  annotation_row <- data.frame(Group = sample_groups[rownames(sampled_data)])
  rownames(annotation_row) <- rownames(sampled_data)  # Ensure row names match
  
  # Generate a unique output filename
  biosample_name <- gsub("\\.rds$", "", basename(file))
  output_png <- file.path(output_dir, paste0(biosample_name, "_heatmap.png"))
  
  # Generate and save heatmap
  png(output_png, height = 6, width = 15, units="in", res=300)
  pheatmap(
    mat = sampled_data,
    show_colnames = FALSE,
    show_rownames = FALSE,
    cluster_rows = FALSE,  # KEEP ROW ORDER FIXED
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    legend = FALSE,
    color = color_scale,
    treeheight_col = 0
  )
  dev.off()
  
  cat("Saved heatmap:", output_png, "\n")
}


##### REG ELEMENTS ENRICHMENTS ####
library(ggpubr)
data <- read.csv("regElementsEnrichments.txt", header=TRUE, sep="\t")
data$Biosample <- sub("unionSetsHg38/(.*)\\.unionSet\\.bed", "\\1", data$Filename2)
data$PropRegElement <- data$OverlapCount / data$LenElements2
pdf("regOverlapPanel.pdf", height=4, width=6)
ggplot(data, aes(x=PropRegElement)) +
  geom_histogram(binwidth = 0.002, fill="#264653") + 
  geom_vline(aes(xintercept = mean(PropRegElement)), color = "black", linetype="dashed", size = 1) +
  annotate("text", x = mean(data$PropRegElement), y = max(table(cut(data$PropRegElement, breaks=20))) * 0.9, 
           label = paste("Mean:", round(mean(data$PropRegElement), 3)), color = "black", size = 5, hjust = -0.1) + 
  theme_classic() + 
  xlab("Proportion of Open Chromatin Elements Overlapping EpiMap ChromHMM Regulatory States") +
  ylab("Number of Biosamples") +
  ggtitle("Regulatory State Overlap of Enformer-Predicted Open Chromatin Elements")
dev.off()
