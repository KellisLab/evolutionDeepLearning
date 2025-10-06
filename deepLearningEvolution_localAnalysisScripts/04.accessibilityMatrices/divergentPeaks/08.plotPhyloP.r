library(dplyr)
library(ggpubr)
library(stringr)
library(gridExtra)


files <- list.files("phyloP", pattern = "\\.txt$", full.names = TRUE)

mergedData <- list()
for (file in files) {
  rawData <- read.delim(file, header = FALSE)
  if (ncol(rawData) >= 5) {
    formattedData <- rawData %>%
      select(V1, V5) %>%
      rename(Chr = V1, PhyloP = V5) %>%
      mutate(Filename = basename(file))  # Extract just the filename
    mergedData[[file]] <- formattedData
  }
}

finalData <- bind_rows(mergedData)
finalData <- finalData %>%
  mutate(
    Condition = str_extract(Filename, "^[^.]+"),  # Everything before first '.'
    PhyloP_Type = str_extract(Filename, "(mammal|primate)PhyloP"),  # Extract "MammalPhyloP" or "PrimatePhyloP"
    State = sapply(str_split(Filename, "\\."), function(x) x[length(x) - 2])
  )

collapsedData <- subset(finalData, 
                        State %in% c("TssA", "OtherRegulatoryStates") & 
                        PhyloP_Type %in% c("mammalPhyloP", "primatePhyloP"))

# Clean up labels for facets
collapsedData$State <- factor(collapsedData$State,
                              levels = c("TssA", "OtherRegulatoryStates"),
                              labels = c("TssA Elements", "Other Regulatory Elements"))

collapsedData$PhyloP_Type <- factor(collapsedData$PhyloP_Type,
                                    levels = c("mammalPhyloP", "primatePhyloP"),
                                    labels = c("Mammal Conservation", "Primate Conservation"))

custom_colors <- c("#E96D4E", "#254857", "#6D6E71", "#000000")
# Create a single facetted plot
collapsedPlot <- ggboxplot(collapsedData, x="Condition", y="PhyloP", fill="Filename", color="Filename", 
                           alpha=0.5, outlier.shape=NA, palette=custom_colors) + 
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() +
  ylim(-1, 1) +
  facet_grid(PhyloP_Type ~ State) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x="", y="PhyloP", title="Conservation Across Element Classes and Lineages")


tssMammalData <- subset(collapsedData, PhyloP_Type == "Mammal Conservation" & State == "TssA Elements")
pairwise.t.test(
  x = tssMammalData$PhyloP,
  g = tssMammalData$Filename,
  p.adjust.method = "bonferroni"
)

otherRegMammalData <- subset(collapsedData, PhyloP_Type == "Mammal Conservation" & State == "Other Regulatory Elements")
pairwise.t.test(
  x = otherRegMammalData$PhyloP,
  g = otherRegMammalData$Filename,
  p.adjust.method = "bonferroni"
)

tssPrimateData <- subset(collapsedData, PhyloP_Type == "Primate Conservation" & State == "TssA Elements")
pairwise.t.test(
  x = tssPrimateData$PhyloP,
  g = tssPrimateData$Filename,
  p.adjust.method = "bonferroni"
)

otherRegPrimateData <- subset(collapsedData, PhyloP_Type == "Primate Conservation" & State == "Other Regulatory Elements")
pairwise.t.test(
  x = otherRegPrimateData$PhyloP,
  g = otherRegPrimateData$Filename,
  p.adjust.method = "bonferroni"
)




######################## PhyloP No TssA/Enh Separation ##################
