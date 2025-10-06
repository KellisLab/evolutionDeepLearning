library(ggpubr)
library(tidyr)
library(dplyr)
library(pheatmap)

data <- read.csv("brain_male_embryo_105_days.chr20.45095781.45096191.txt", header=FALSE, sep="\t")
#data <- read.csv("heart_embryo_80_days.chr1.201381090.201381626.txt", header=FALSE, sep="\t")
#data <- read.csv("hepatocyte_from_H9.chr1.201381602.201382791.txt", header=FALSE, sep="\t")
head(data)
haplotypes <- data$V1
data <- select(data, -c(V1))
mat <- as.matrix(data)
rownames(mat) <- haplotypes

manual_order <- c(
  "NA18858.A", "NA18858.B", "NA20517.A", "NA20517.B", "HG00525.A", "HG00525.B", 
  "HG02615.A", "HG02615.B", "HG03079.A", "HG03079.B",
  "vindija.neanderthal.A", "vindija.neanderthal.B", 
  "chagyrskaya.neanderthal.A", "chagyrskaya.neanderthal.B",
  "altai.Neanderthal.A", "altai.Neanderthal.B", "denisova.A", "denisova.B",
  "clara.panTroTro.A", "clara.panTroTro.B", "jimmie.panTroVerus.A", "jimmie.panTroVerus.B",
  "taweh.panTroEllioti.A", "taweh.panTroEllioti.B", 
  "nakuu.panTroSchweinfurthii.A", "nakuu.panTroSchweinfurthii.B",
  "natalie.panPaniscus.A", "natalie.panPaniscus.B",
  "dian.gorGorGor.A", "dian.gorGorGor.B", "victoria.gorBeringeiGraueri.A", "victoria.gorBeringeiGraueri.B",
  "dunja.ponAbe.A", "dunja.ponAbe.B"
)

# Reorder the matrix
mat <- mat[manual_order, ]

central_width <- 5000
num_columns <- ncol(mat)
start_index <- max(1, floor((num_columns - central_width) / 2) + 1)
end_index <- min(num_columns, start_index + central_width - 1)
mat_central <- mat[, start_index:end_index]

pheatmap(mat_central, cluster_cols=FALSE, show_colnames = FALSE, cluster_rows = FALSE, color = colorRampPalette(c("white", "red"))(100))


group_assignments <- data.frame(
  haplotype = c(
    "NA18858.A", "NA18858.B", "NA20517.A", "NA20517.B", "HG00525.A", "HG00525.B",
    "HG02615.A", "HG02615.B", "HG03079.A", "HG03079.B", # Modern Human
    "vindija.neanderthal.A", "vindija.neanderthal.B", "chagyrskaya.neanderthal.A",
    "chagyrskaya.neanderthal.B", "altai.Neanderthal.A", "altai.Neanderthal.B",
    "denisova.A", "denisova.B", # Archaic Hominin
    "clara.panTroTro.A", "clara.panTroTro.B", "jimmie.panTroVerus.A", "jimmie.panTroVerus.B",
    "taweh.panTroEllioti.A", "taweh.panTroEllioti.B", "nakuu.panTroSchweinfurthii.A",
    "nakuu.panTroSchweinfurthii.B", "natalie.panPaniscus.A", "natalie.panPaniscus.B",
    "dian.gorGorGor.A", "dian.gorGorGor.B", "victoria.gorBeringeiGraueri.A",
    "victoria.gorBeringeiGraueri.B", "dunja.ponAbe.A", "dunja.ponAbe.B" # Great Ape
  ),
  group = c(
    rep("Modern Human", 10),
    rep("Archaic Hominin", 8),
    rep("Great Ape", 16)
  )
)

data_long <- as.data.frame(mat) 
data_long$haplotype <- rownames(mat)  
data_long <- pivot_longer(
  data_long,
  cols = -haplotype,
  names_to = "position",
  values_to = "accessibility"
) %>%
  left_join(group_assignments, by = "haplotype") %>%  # Add group info
  mutate(
    position = as.numeric(gsub("V", "", position)),  # Convert position to numeric
    center = (max(position) + min(position)) / 2,   # Calculate center
    position = position - center                    # Make position relative to center
  )

# Ensure haplotypes are ordered as desired
data_long$haplotype <- factor(data_long$haplotype, levels = group_assignments$haplotype)

#pdf("brain_male_embryo_105_days.chr20.45095781.45096191.pdf", height=16, width=24)
png("brain_male_embryo_105_days.chr20.45095781.45096191.png", height=1600, width=2400, units="px")
ggplot(data_long, aes(x = position, y = accessibility, group = haplotype, fill = group)) +
  geom_area(alpha = 1) +  # Fill under the curve
  facet_wrap(~ haplotype, ncol = 1, strip.position = "right") +
  scale_fill_manual(
    values = c("Modern Human" = "#F4A361", "Archaic Hominin" = "#EAC466", "Great Ape" = "#27A092")
  ) +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
  ) +
  labs(
    x = "Distance from Peak Center",
    y = "Predicted Accessibility",
    title = "Chromatin Accessibility Across Haplotypes",
    fill = "Group"
  )
dev.off()
