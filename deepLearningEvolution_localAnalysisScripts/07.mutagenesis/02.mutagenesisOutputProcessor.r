library(ggplot2)
library(ggseqlogo)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(optparse)
library(pheatmap)
library(tidyverse)
library(universalmotif)
library(gridExtra)
#BiocManager::install("TFBSTools")
library(TFBSTools) #https://doi.org/10.1093/bioinformatics/btw024
#BiocManager::install("JASPAR2022")
library(JASPAR2022)
library(Biostrings)
library(ggpubr)

read_mutagenesis_data <- function(inputFile) {
  data_df <- read.delim(inputFile, header = TRUE)
  data_df$Amplitude <- as.numeric(as.character(data_df$Amplitude))
  region_name <- unique(data_df$Name)
  if (length(region_name) != 1) {
    stop("Error: The input file should contain data for only one region.")
  }
  region_data <- data_df[, c("Position", "Mutation", "Amplitude")]
  wide_data <- region_data %>%
    pivot_wider(names_from = Position, values_from = Amplitude) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate_all(~replace(., is.nan(.), 0)) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  mat_ready <- wide_data[, -1]  # Exclude mutation column
  data_matrix <- as.matrix(mat_ready)
  rownames(data_matrix) <- c('A', 'C', 'G', 'T')
  return(data_matrix)
}

plot_mutagenesis_matrix <- function(data_matrix) {
  color_palette <- colorRampPalette(c("blue", "#EEEEEE", "red"))(100)
  max_abs_val <- max(abs(data_matrix))
  breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)
  p <- pheatmap(
    data_matrix, 
    cluster_cols = FALSE, 
    cluster_rows = FALSE, 
    border_color = NA,
    show_colnames = FALSE, 
    color = color_palette,
    breaks = breaks  # Ensures white is at 0
  )
  return(p)
}

compare_motif_scores <- function(pfm_list, motif_ids, ape_seq, human_seq, min_score = "0%") {
  results <- data.frame(
    motif_id = character(),
    motif_name = character(),
    ape_best_score = numeric(),
    human_best_score = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (id in motif_ids) {
    pfm <- pfm_list[[id]]
    pwm <- toPWM(pfm)
    
    ape_hits <- searchSeq(pwm, ape_seq, min.score = min_score, strand = "*")
    human_hits <- searchSeq(pwm, human_seq, min.score = min_score, strand = "*")
    
    ape_best <- if (length(ape_hits) > 0) max(score(ape_hits)) else NA_real_
    human_best <- if (length(human_hits) > 0) max(score(human_hits)) else NA_real_
    
    results <- rbind(results, data.frame(
      motif_id = id,
      motif_name = pfm@name,
      ape_best_score = ape_best,
      human_best_score = human_best,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

find_ape_only_hits <- function(pfm_list, ape_seq, human_seq, min_score = "85%") {
  hits_in_ape_only <- list()
  for (id in names(pfm_list)) {
    motif <- toPWM(pfm_list[[id]])
    ape_hit <- searchSeq(motif, ape_seq, min.score = min_score, strand = "*")
    human_hit <- searchSeq(motif, human_seq, min.score = min_score, strand = "*")
    if (length(ape_hit) > 0 && length(human_hit) == 0) {
      hits_in_ape_only[[id]] <- ape_hit
    }
  }
  return(hits_in_ape_only)
}

find_human_only_hits <- function(pfm_list, ape_seq, human_seq, min_score = "85%") {
  hits_in_human_only <- list()
  for (id in names(pfm_list)) {
    motif <- toPWM(pfm_list[[id]])
    human_hit <- searchSeq(motif, human_seq, min.score = min_score, strand = "*")
    ape_hit <- searchSeq(motif, ape_seq, min.score = min_score, strand = "*")
    if (length(human_hit) > 0 && length(ape_hit) == 0) {
      hits_in_human_only[[id]] <- human_hit
    }
  }
  return(hits_in_human_only)
}

subset_mutagenesis_matrix <- function(data_matrix, center_pos, window_size) {
  positions <- as.numeric(colnames(data_matrix))
  
  if (!(center_pos %in% positions)) {
    stop("Center position not found in matrix column names.")
  }
  
  center_index <- which(positions == center_pos)
  
  # Half-open: start inclusive, end exclusive
  start_index <- center_index - window_size
  end_index <- center_index + window_size
  
  if (start_index < 1 || end_index > ncol(data_matrix) + 1) {
    stop("Requested window is out of bounds.")
  }
  
  subset_matrix <- data_matrix[, start_index:(end_index - 1)]
  return(subset_matrix)
}

subset_mutagenesis_matrix_start_end <- function(data_matrix, start_pos, end_pos) {
  positions <- as.numeric(colnames(data_matrix))
  
  # Validate range
  if (!(start_pos %in% positions)) {
    stop("Start position not found in matrix column names.")
  }
  if (!((end_pos - 1) %in% positions)) {
    stop("End position - 1 not found in matrix column names.")
  }
  
  # Find column indices corresponding to positions
  start_index <- which(positions == start_pos)
  end_index <- which(positions == (end_pos - 1)) + 1  # make end exclusive
  
  if (length(start_index) == 0 || length(end_index) == 0) {
    stop("Start or end index could not be resolved.")
  }
  
  subset_matrix <- data_matrix[, start_index:(end_index - 1)]
  return(subset_matrix)
}

revcomp_pfm <- function(pfm) {
  freq_mat <- pfm@profileMatrix
  # Reverse columns (motif length)
  freq_mat <- freq_mat[, ncol(freq_mat):1]
  
  # Complement rows: swap A<->T, C<->G
  comp_order <- c("T", "G", "C", "A")  # target row order
  rownames(freq_mat) <- c("A", "C", "G", "T")
  freq_mat <- freq_mat[comp_order, ]
  
  return(freq_mat)
}

################################################################################
######################  Analysis Starts Here  ##################################
################################################################################

#inputFile <- "satMutResults/chr20.45095781.45096191.satMut.txt"
inputFile <- "satMutResults/chr1.201381602.201382791.track.347.txt"

data_matrix <- read_mutagenesis_data(inputFile)
p <- plot_mutagenesis_matrix(data_matrix)
opts <- list(species = 9606, collection = "CORE", all_versions = FALSE)
pfm_list <- getMatrixSet(JASPAR2022, opts)
human_seq <- DNAString("ACACCCACAAGGAGTCAG")
ape_seq   <- DNAString("ACACCCACAATGAGTCAG")
ape_only_hits <- find_ape_only_hits(pfm_list,ape_seq, human_seq)
score_df <- compare_motif_scores(pfm_list, names(ape_only_hits), ape_seq, human_seq)

ggscatter(score_df, x="human_best_score", y="ape_best_score", label="motif_name")+
  xlim(0,15)+
  ylim(0,15)+
  geom_abline(slope=1, intercept=0, linetype="dashed")

pdf("pdf/MA1135.1.pdf", height=4, width=8)
ggseqlogo(pfm_list$MA1135.1@profileMatrix)
dev.off()

subset_matrix <- subset_mutagenesis_matrix(data_matrix, 201382001, 50)
p <- plot_mutagenesis_matrix(subset_matrix)

pdf("pdf/hepatocyte.ZoomIn.pdf", height=2, width=8)
p
dev.off()

###### HEART ####
inputFile <- "satMutResults/heart_embryo_80_days.chr1.201381090.201381626.bed.track.570.txt"

data_matrix <- read_mutagenesis_data(inputFile)
p <- plot_mutagenesis_matrix(data_matrix)
pdf("pdf/heart.full.pdf", height=3, width=14)
p
dev.off()

subset_matrix <- subset_mutagenesis_matrix_start_end(data_matrix, 201381282, 201381343)
# subset range : 201381262 201381363
p <- plot_mutagenesis_matrix(subset_matrix)
pdf("pdf/heart.zoomIn.pdf", height=4, width=12)
p
dev.off()

# TFBS
human_seq <- DNAString("CTCCTGGCCATGGCAGCCGGCC")
ape_seq <- DNAString("CTCCTGGCCACGGCAGCCGGCC")
opts <- list(species = 9606, collection = "CORE", all_versions = FALSE)
pfm_list <- getMatrixSet(JASPAR2022, opts)
human_only_hits <- find_human_only_hits(pfm_list,ape_seq, human_seq)

score_df <- compare_motif_scores(pfm_list, names(human_only_hits), ape_seq, human_seq)

ggscatter(score_df, x="human_best_score", y="ape_best_score", label="motif_name")+
  geom_abline(slope=1, intercept=0, linetype="dashed") + 
  labs(x="Motif Score, Human Allele", y="Motif Score, Ape Allele")

currentMotif <- pfm_list$MA0775.1 # MEIS3
currentMotif <- pfm_list$MA0904.2 # HOXB5
currentMotif <- pfm_list$MA1495.1 #HOXA1
rcMat <- revcomp_pfm(currentMotif)
pdf("pdf/MEIS_motif.pdf", height=4, width=8)
ggseqlogo(currentMotif@profileMatrix)
dev.off()
ggseqlogo(rcMat)



##############################################################################
#################### Brain Analysis ##########################################
##############################################################################
opts <- list(species = 9606, collection = "CORE", all_versions = FALSE)
pfm_list <- getMatrixSet(JASPAR2022, opts)

cacna1c_file <- "satMutResults/chr12.1933407.1934070.track.370.txt"
cacna1c_data_matrix <- read_mutagenesis_data(cacna1c_file)
p_cacna1c <- plot_mutagenesis_matrix(cacna1c_data_matrix)

grid1_file <- "satMutResults/chr10.86261934.86262389.track.370.txt"
grid1_data_matrix <- read_mutagenesis_data(grid1_file)
p_grid1 <- plot_mutagenesis_matrix(grid1_data_matrix)

############ GRIN3B #############
grin3b_file <- "satMutResults/chr19.1795882.1796497.track.0.txt"
grin3b_data_matrix <- read_mutagenesis_data(grin3b_file)
p_grin3b <- plot_mutagenesis_matrix(grin3b_data_matrix)
grin3b_subset_data_matrix <- subset_mutagenesis_matrix_start_end(grin3b_data_matrix, 1795939, 1795999)
p_grin3b_subset <- plot_mutagenesis_matrix(grin3b_subset_data_matrix)
pdf("pdf/grin3b_data_matrix.1795939.1795999.pdf", height=3, width=12)
p_grin3b_subset
dev.off()

grin_hum_seq <- DNAString("GTACAGGCAGATGGCGCCCTCTGGCTTT")
grin_ape_seq <- DNAString("GTACAGGCAGATGGTGCCCTCTGGCTTT")
grin_human_only_hits <- find_human_only_hits(pfm_list, grin_ape_seq, grin_hum_seq)
grin_ape_only_hits <- find_ape_only_hits(pfm_list,grin_ape_seq, grin_hum_seq)
score_df <- compare_motif_scores(pfm_list, names(grin_human_only_hits), grin_ape_seq, grin_hum_seq)
ggscatter(score_df, x="human_best_score", y="ape_best_score", label="motif_name")+
  geom_abline(slope=1, intercept=0, linetype="dashed") + 
  labs(x="Motif Score, Human Allele", y="Motif Score, Ape Allele") + xlim(0, 12) + ylim(0, 12)
e2f1_motif <- pfm_list$MA0024.3 # e2f1

pdf("pdf/e2f1.motif.pdf", height=4, width=8)
ggseqlogo(e2f1_motif@profileMatrix)
dev.off()

ggseqlogo(e2f1_motif@profileMatrix)
knl1_file <- "satMutResults/chr15.40266842.40267256.track.77.txt"
knl1_data_matrix <- read_mutagenesis_data(knl1_file)
p_knl1 <- plot_mutagenesis_matrix(knl1_data_matrix)
knl1_subset_data_matrix <- subset_mutagenesis_matrix_start_end(knl1_data_matrix, 40266952, 40267052)
p_knl1_subset <- plot_mutagenesis_matrix(knl1_subset_data_matrix)
pdf("pdf/knl1_data_matrix.40266952.40267052.pdf", height=3, width=12)
p_knl1_subset
dev.off()
knl1_hum_seq <- DNAString("AGAGGAGGAGACTGATAGAGAGAGGGTGAGTCACCCCCAGGCATGATGGAGAGTGAGTCA")
knl1_ape_seq <- DNAString("AGAGGAGGAGACTGATAGAGAGAGGGTGAGACACCCCCAGGCATGATGGAGAGTGAGTCA")
knl_human_only_hits <- find_human_only_hits(pfm_list, knl1_ape_seq, knl1_hum_seq)
knl1_score_df <- compare_motif_scores(pfm_list, names(knl_human_only_hits), knl1_ape_seq, knl1_hum_seq)
ggscatter(knl1_score_df, x="human_best_score", y="ape_best_score", label="motif_name")+
  geom_abline(slope=1, intercept=0, linetype="dashed") + 
  labs(x="Motif Score, Human Allele", y="Motif Score, Ape Allele")
fos_motif <- pfm_list$MA0476.1
pdf("pdf/fos.motif.pdf", height=4, width=8)
ggseqlogo(fos_motif@profileMatrix)
dev.off()


nedd9_file <- "satMutResults/chr6.163756292.163756452.track.179.txt"
nedd9_data_matrix <- read_mutagenesis_data(nedd9_file)
p_nedd9 <- plot_mutagenesis_matrix(nedd9_data_matrix)



