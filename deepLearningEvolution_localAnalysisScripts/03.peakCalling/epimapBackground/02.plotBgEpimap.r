library(ggpubr)
library(dplyr)
library(stringr)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(gridExtra)

file_list <- list.files(pattern = "*.txt")  # Adjust pattern if needed

# Initialize list for heatmaps
heatmap_list <- list()

# Read and process data from all files
data_list <- lapply(file_list, function(file) {
  data <- read.csv(file, header=FALSE, sep="\t") %>%
    select(c("V2", "V3", "V4", "V5", "V6", "V9", "V10", "V11")) %>%
    rename(Filename1 = V2, Filename2 = V3, LenElements1 = V4, LenElements2 = V5, 
           OverlapCount = V6, Enrichment = V9, EnrichPValue = V10, DepletePValue = V11) %>%
    mutate(
      log2Enrich = log2(Enrichment),
      pRaw = pmin(EnrichPValue, DepletePValue) * 2,
      pAdj = p.adjust(pRaw, method="fdr"),
      EpimapTissue = str_extract(basename(Filename1), "^[^.]+"),
      EpimapState = str_extract(Filename1, "(?<=\\.)[^.]+(?=\\.bed)"),
      EnformerBiosample = str_extract(basename(Filename2), "^[^.]+")
    )
  return(data)
})

# Combine all data into a single dataframe
all_data <- bind_rows(data_list)

# Extract EnhA1 data for clustering
enha1_data <- all_data %>% filter(EpimapState == "EnhA1")

# Reshape EnhA1 data for clustering
enha1_heatmap <- enha1_data %>%
  select(EpimapTissue, EnformerBiosample, log2Enrich) %>%
  pivot_wider(names_from = EnformerBiosample, values_from = log2Enrich)

# Set row names and remove EpimapTissue column
enha1_heatmap <- as.data.frame(enha1_heatmap)
rownames(enha1_heatmap) <- enha1_heatmap$EpimapTissue
enha1_heatmap <- select(enha1_heatmap, -c("EpimapTissue"))

# Convert to numeric matrix
enha1_matrix <- as.matrix(enha1_heatmap)

# Perform clustering to get row and column order
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)
enha1_clustering <- pheatmap(t(enha1_matrix), cluster_rows = TRUE, cluster_cols = TRUE, silent = TRUE, color = custom_colors, 
                             breaks = seq(-4, 4, length.out = 101))
row_order <- enha1_clustering$tree_row$order
col_order <- enha1_clustering$tree_col$order

# Generate heatmaps for all states using the same order
unique_states <- unique(all_data$EpimapState)

for (state in unique_states) {
  state_data <- all_data %>% filter(EpimapState == state)
  
  # Reshape for heatmap
  heatmap_data <- state_data %>%
    select(EpimapTissue, EnformerBiosample, log2Enrich) %>%
    pivot_wider(names_from = EnformerBiosample, values_from = log2Enrich)
  
  # Set row names and remove EpimapTissue column
  heatmap_data <- as.data.frame(heatmap_data)
  rownames(heatmap_data) <- heatmap_data$EpimapTissue
  heatmap_data <- select(heatmap_data, -c("EpimapTissue"))
  
  # Convert to numeric matrix
  heatmap_matrix <- as.matrix(heatmap_data)
  
  # Reorder rows and columns based on EnhA1 clustering
  heatmap_matrix <- heatmap_matrix[col_order, row_order]
  
  # Define color scale
  custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Generate pheatmap
  heatmap_plot <- pheatmap(t(heatmap_matrix), 
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           main = state, 
                           fontsize_row = 5, 
                           color = custom_colors, 
                           breaks = seq(-4, 4, length.out = 101),
                           silent = TRUE, legend = FALSE)
  
  # Store heatmap
  heatmap_list[[state]] <- heatmap_plot
}

pdf("epiMapBackground.Unlabeled.pdf", height=10, width=13)
grid.arrange(heatmap_list[["EnhA1"]]$gtable, 
heatmap_list[["EnhA2"]]$gtable, 
heatmap_list[["EnhBiv"]]$gtable, 
heatmap_list[["EnhG1"]]$gtable, 
heatmap_list[["EnhG2"]]$gtable, 
heatmap_list[["Het"]]$gtable,
heatmap_list[["TssA"]]$gtable,
heatmap_list[["TssBiv"]]$gtable,
heatmap_list[["Tx"]]$gtable,
heatmap_list[["ZNF_Rpts"]]$gtable,
heatmap_list[["ReprPC"]]$gtable,
heatmap_list[["Quies"]]$gtable,
nrow=2)
dev.off()

pdf("labeled.order.pdf", height=8, width=8)
enha1_clustering
dev.off()
