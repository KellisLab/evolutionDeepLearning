library(ggplot2)
library(ggpubr)
library(ggridges)
myTheme <- theme(legend.position = "none", axis.title.x = element_blank())
fileTable <- read.table("~/kellisUrop/samCoverage/coverageFileTable.txt", header = TRUE, sep = "\t")

for (i in 1:nrow(fileTable)) {
  name <- fileTable$Name[i]
  filepath <- as.character(fileTable$Filename[i])
  data <- read.table(filepath, header=TRUE, sep = "\t")
  assign(name, subset(data, Group == "Empirical"))
}

plotList <- list(Altai, Chagyrskaya, Denisova, panTroEllioti, panTroSchweinfurthii, 
                   panTroTro, panTroVerus, Vindija, Gambian, Yoruba, hanChinese, Mende, Toscani,
                 panPaniscus, gorGorGor, gorBeringeiGraueri, pongoAbelii)
plotNames <- c("altai.Neanderthal", "chagyrskaya.neanderthal", "denisova", "taweh.panTroEllioti", "nakuu.panTroSchweinfurthii", 
                  "clara.panTroTro", "jimmie.panTroVerus", "vindija.neanderthal", "HG02615 (Gambian)", "NA18858 (Yoruba)", 
               "HG00525 (hanChinese)", "HG03079 (Mende)", "NA20517 (Toscani)", "natalie.panPaniscus", "dian.gorGorGor", 
               "victoria.gorBeringeiGraueri", "dunja.pongoAbelii")
for (j in 1:length(plotList)) {
    plotList[[j]]$Sample <- plotNames[j]
}

composite <- do.call(rbind, plotList)
composite <- subset(composite, select = -c(Group))
composite <- subset(composite, Coverage < 100)

mean_coverage <- by(composite, composite$Sample, FUN = function(sub_df) {
  weighted.mean(sub_df$Coverage, w = sub_df$Pileups)
})
means_df <- as.data.frame(mean_coverage)
means_df$Sample <- rownames(means_df)

custom_order <- c("victoria.gorBeringeiGraueri", "dian.gorGorGor", "dunja.pongoAbelii", "vindija.neanderthal", "denisova", 
                  "altai.Neanderthal", "chagyrskaya.neanderthal", "HG03079 (Mende)", "HG02615 (Gambian)", "NA18858 (Yoruba)", 
                  "NA20517 (Toscani)", "HG00525 (hanChinese)", "nakuu.panTroSchweinfurthii", "clara.panTroTro", "jimmie.panTroVerus", 
                  "taweh.panTroEllioti", "natalie.panPaniscus")

composite$Sample <- factor(composite$Sample, levels = custom_order)
composite <- composite[order(composite$Sample), ]

ggplot(composite, aes(x = Coverage, y = Sample, height = Pileups, fill = Pileups, group = Sample)) +
  geom_density_ridges_gradient(stat = "identity", scale = 0.9, rel_min_height = 0.01) +
  theme_classic() + theme(legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 20, l = 10)) +
  scale_fill_gradient(low = "#AAC5FF", high = "#0439AA") +
  labs(xlab = "Coverage")+
  geom_segment(data= subset(means_df, Sample == "victoria.gorBeringeiGraueri"), aes(x=x, xend = x, y = 1, yend = 1.35, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "dian.gorGorGor"), aes(x=x, xend = x, y = 2, yend = 2.55, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "dunja.pongoAbelii"), aes(x=x, xend = x, y = 3, yend = 3.42, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "vindija.neanderthal"), aes(x=x, xend = x, y = 4, yend = 4.4, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "denisova"), aes(x=x, xend = x, y = 5, yend = 5.68, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "altai.Neanderthal"), aes(x=x, xend = x, y = 6, yend = 6.38, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "chagyrskaya.neanderthal"), aes(x=x, xend = x, y = 7, yend = 7.72, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "HG03079 (Mende)"), aes(x=x, xend = x, y = 8, yend = 8.77, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "HG02615 (Gambian)"), aes(x=x, xend = x, y = 9, yend = 9.7, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "NA18858 (Yoruba)"), aes(x=x, xend = x, y = 10, yend = 10.7, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "NA20517 (Toscani)"), aes(x=x, xend = x, y = 11, yend = 11.75, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "HG00525 (hanChinese)"), aes(x=x, xend = x, y = 12, yend = 12.77, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "nakuu.panTroSchweinfurthii"), aes(x=x, xend = x, y = 13, yend = 13.5, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "clara.panTroTro"), aes(x=x, xend = x, y = 14, yend = 14.55, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "jimmie.panTroVerus"), aes(x=x, xend = x, y = 15, yend = 15.52, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "taweh.panTroEllioti"), aes(x=x, xend = x, y = 16, yend = 16.9, color="red", group = "Sample"), inherit.aes = FALSE)+
  geom_segment(data= subset(means_df, Sample == "natalie.panPaniscus"), aes(x=x, xend = x, y = 17, yend = 17.6, color="red", group = "Sample"), inherit.aes = FALSE)


    
  