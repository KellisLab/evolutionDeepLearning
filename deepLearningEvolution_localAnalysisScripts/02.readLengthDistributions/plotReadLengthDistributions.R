library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggridges)
library(patchwork)
myTheme <- theme(legend.position = "none", axis.title.x = element_blank())

directory <- "readLengthDistributions"
myTheme <- theme(legend.position = "none", axis.title.x = element_blank())
filePaths <- list.files(directory, full.names = TRUE)
dataList <- list()
for (filepath in filePaths) {
  name <- tools::file_path_sans_ext(basename(filepath))
  tryCatch({
    dataList[[name]] <- read.table(filepath, header = TRUE, sep = "\t")
  }, error = function(e) {
    message(sprintf("Error reading file: %s. Skipping...", filepath))
  })
}
names(dataList) <- gsub("\\.readLengthDistribution$", "", names(dataList))

#### Riley's attempt ###
composite_list <- lapply(names(dataList), function(sample_name) {
  # Add a Sample column to the current data frame
  df <- dataList[[sample_name]]
  df$Sample <- sample_name
  return(df)
})

# Combine all data frames into one composite data frame
composite_df <- do.call(rbind, composite_list)

# add group identifiers
composite_df$Group <- ifelse(
  grepl("Neanderthal|denisova", composite_df$Sample, ignore.case = TRUE), "ArchaicHominin",
  ifelse(
    grepl("panTro|panPaniscus", composite_df$Sample, ignore.case = TRUE), "PanGenus",
    ifelse(
      grepl("gor|pon", composite_df$Sample, ignore.case = TRUE), "OutGroup",
      ifelse(
        grepl("^HG|^NA", composite_df$Sample, ignore.case = TRUE), "ModernHuman",
        NA # Catch-all in case of unexpected sample names
      )
    )
  )
)

breaks <- seq(0, 200, by = 10)  # Adjust the range as needed
labels <- c(paste0(breaks[-length(breaks)], "-", breaks[-1]), ">200")  # Create labels for the buckets
breaks <- c(breaks, Inf)  # Append a final 'Inf' bucket

# Create a bucketed data frame based on ReadLength
bucketed_df <- composite_df %>%
  mutate(readLengthBucket = cut(ReadLength, breaks = breaks, labels = labels, right = FALSE)) %>%
  group_by(readLengthBucket, Group) %>%
  summarize(Count = sum(Count), .groups = 'drop')

# Define a custom order for the groups
custom_order <- c("PanGenus", "ModernHuman", "ArchaicHominin", "OutGroup")
bucketed_df$Group <- factor(bucketed_df$Group, levels = custom_order)
bucketed_df <- bucketed_df[order(bucketed_df$Group), ]

# Calculate total counts for each group
total_counts <- bucketed_df %>%
  group_by(Group) %>%
  summarize(Total = sum(Count), .groups = 'drop')

# Join total counts back to the bucketed data frame and calculate percentages
proportional_df <- bucketed_df %>%
  left_join(total_counts, by = "Group") %>%
  mutate(Percentage = (Count / Total) * 100) %>%
  select(readLengthBucket, Group, Percentage)

# Define group colors
group_colors <- c(
  "PanGenus" = "#264653",
  "ModernHuman" = "#E9C46A",
  "ArchaicHominin" = "#2A9D8F",
  "OutGroup" = "#E76F51"
)

ggplot(proportional_df, aes(x = readLengthBucket, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "ReadLength Distributions Across Sample Groups",
       x = "Read Length",
       y = "Percentage") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Group, ncol = 1)


### Riley sample level ###
aggregated_df <- composite_df %>%
  group_by(ReadLength, Sample, Group) %>%
  summarize(Count = sum(Count), .groups = 'drop')

# Define breaks and labels
breaks <- seq(0, 200, by = 10)  # Adjust the range as needed
labels <- c(paste0(breaks[-length(breaks)], "-", breaks[-1]), ">200")
breaks <- c(breaks, Inf)  # Append a >200 bucket

# Create a bucketed data frame
bucketed_df <- aggregated_df %>%
  mutate(readLengthBucket = cut(ReadLength, breaks = breaks, labels = labels, right = FALSE)) %>%
  group_by(readLengthBucket, Sample, Group) %>%
  summarize(Count = sum(Count), .groups = 'drop')

# Calculate total counts for each sample
total_counts <- bucketed_df %>%
  group_by(Sample) %>%
  summarize(Total = sum(Count), .groups = 'drop')

# Create proportional data
proportional_df <- bucketed_df %>%
  left_join(total_counts, by = "Sample") %>%
  mutate(Percentage = (Count / Total) * 100) %>%
  select(readLengthBucket, Sample, Group, Percentage)

# Sample-specific colors
sample_colors <- c(
  "clara.panTroTro" = "#264653",
  "jimmie.panTroVerus" = "#264653",
  "taweh.panTroEllioti" = "#264653",
  "nakuu.panTroSchweinfurthii" = "#264653",
  "natalie.panPaniscus" = "#264653",
  "HG02615" = "#E9C46A",
  "NA18858" = "#E9C46A",
  "HG00525" = "#E9C46A",
  "HG03079" = "#E9C46A",
  "NA20517" = "#E9C46A",
  "altai.Neanderthal" = "#2A9D8F",
  "chagyrskaya.neanderthal" = "#2A9D8F",
  "denisova" = "#2A9D8F",
  "vindija.neanderthal" = "#2A9D8F",
  "dian.gorGorGor" = "#E76F51",
  "victoria.gorBeringeiGraueri" = "#E76F51",
  "dunja.ponAbe" = "#E76F51"
)

# Create a plot for each group
plot_modernhuman <- ggplot(proportional_df %>% filter(Group == "ModernHuman"), aes(x = readLengthBucket, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(x = "Read Length", y = "Percent of Reads") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.5)), legend.position = "none") +
  facet_wrap(~Sample, ncol = 5) +
  ggtitle("ModernHuman")

plot_pangenus <- ggplot(proportional_df %>% filter(Group == "PanGenus"), aes(x = readLengthBucket, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(x = "Read Length", y = "Percent of Reads") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.5)), legend.position = "none") +
  facet_wrap(~Sample, ncol = 5) +
  ggtitle("PanGenus")

plot_archaichominin <- ggplot(proportional_df %>% filter(Group == "ArchaicHominin"), aes(x = readLengthBucket, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(x = "Read Length", y = "Percent of Reads") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.5)), legend.position = "none") +
  facet_wrap(~Sample, ncol = 4) +
  ggtitle("ArchaicHominin")

plot_outgroup <- ggplot(proportional_df %>% filter(Group == "OutGroup"), aes(x = readLengthBucket, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(x = "Read Length", y = "Percent of Reads") +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.5)), legend.position = "none") +
  facet_wrap(~Sample, ncol = 3) +
  ggtitle("OutGroup")

pdf("readLengthDistributions.pdf", height=5, width=12)
grid.arrange(plot_modernhuman, plot_archaichominin, plot_pangenus, plot_outgroup)
dev.off()



