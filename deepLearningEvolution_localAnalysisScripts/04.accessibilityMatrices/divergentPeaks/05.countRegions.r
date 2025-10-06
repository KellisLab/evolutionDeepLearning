library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(ggrepel)

# Define directories
dirs <- list(
  BackgroundRegions = "backgroundRegions",
  UpInHominins = "upInHominins",
  UpInGreatApes = "upInGreatApes"
)

count_lines <- function(file) {
  count_str <- system(paste("wc -l", shQuote(file)), intern = TRUE)
  count <- as.integer(str_trim(str_extract(count_str, "^\\s*\\d+")))  # Trim and extract the number
  return(count)
}
results <- list()
for (category in names(dirs)) {
  files <- list.files(dirs[[category]], pattern = "\\.bed$", full.names = TRUE)
  for (file in files) {
    biosample <- str_extract(basename(file), "^[^.]+")
    line_count <- count_lines(file)
    results <- append(results, list(data.frame(Biosample = biosample, Category = category, Count = line_count)))
  }
}

df <- bind_rows(results) %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)  # Fill missing values with 0

data <- as.data.frame(df)
data$Unchanged <- data$BackgroundRegions - data$UpInHominins - data$UpInGreatApes
data$PropUnchanged <- data$Unchanged / data$BackgroundRegions
data$PropUpInHominins <- data$UpInHominins / data$BackgroundRegions
data$PropUpInGreatApes <- data$UpInGreatApes / data$BackgroundRegions
data$TotalPropChanged <- data$PropUpInHominins + data$PropUpInGreatApes
data$Biosample <- factor(data$Biosample, levels = data$Biosample[order(data$TotalPropChanged, decreasing = TRUE)])


tallData <- data %>%
  pivot_longer(cols = c(PropUpInHominins, PropUpInGreatApes), names_to = "Category", values_to = "Prop")
ggplot(tallData, aes(x=Biosample, y=Prop, fill=Category)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("PropUpInHominins" = "#DC143C",
                               "PropUpInGreatApes" = "#1E3A8A")) +
  theme_classic() + 
  coord_flip() + 
  labs(title="Hominin-Great Ape Divergent Regulatory Elements across Biosamples", y="Proportion of Open Chromatin Elements") + 
  theme(axis.text.y=element_text(size=6))




ggscatter(data, "PropUpInHominins", "PropUpInGreatApes", label="Biosample") + geom_abline(intercept=0, slope=1, linetype="dashed")


plot_data <- data %>%
  mutate(UpInGreatApes = -UpInGreatApes)
plot_data <- plot_data %>%
  arrange(desc(abs(UpInHominins) + abs(UpInGreatApes)))
pdf("pdf/homininGreatApeDACounts.pdf", height=6, width=8)
ggplot(plot_data, aes(y = reorder(Biosample, UpInHominins + abs(UpInGreatApes)))) +
  geom_bar(aes(x = UpInHominins), stat = "identity", fill = "#E96D4E") +  # Hominin gains
  geom_bar(aes(x = UpInGreatApes), stat = "identity", fill = "#244857") + # Great Ape gains
  theme_classic() +
  labs(x = "Number of Differentially Accessible Regions", y = "Biosample",
       title = "Directional Differential Accessibility in Hominins vs. Great Apes")
dev.off()

plot_data <- data %>%
  select(Biosample, PropUpInHominins, PropUpInGreatApes) %>%
  pivot_longer(cols = c(PropUpInHominins, PropUpInGreatApes),
               names_to = "Comparison",
               values_to = "Proportion") %>%
  mutate(Proportion = ifelse(Comparison == "PropUpInGreatApes", -Proportion, Proportion))
plot_data <- plot_data %>%
  arrange(desc(abs(PropUpInHominins) + abs(PropUpInGreatApes)))

# Plot proportions
ggplot(plot_data, aes(y = reorder(Biosample, Proportion), x = Proportion, fill = Comparison)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("PropUpInHominins" = "#E96D4E", "PropUpInGreatApes" = "#244857")) +
  theme_minimal() +
  labs(x = "Proportion of Differentially Accessible Regions", y = "Biosample",
       title = "Proportion of Differentially Accessible Regions in Hominins vs. Great Apes") +
  theme(legend.position = "none", axis.text.y = element_text(size = 7))

























################# MODERN ARCHAIC ###################
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

dirs <- list(
  BackgroundRegions = "backgroundRegions",
  UpInModern = "UpInModern",
  UpInArchaic = "UpInArchaic"
)

count_lines <- function(file) {
  count_str <- system(paste("wc -l", shQuote(file)), intern = TRUE)
  count <- as.integer(str_trim(str_extract(count_str, "^\\s*\\d+")))  # Trim and extract the number
  return(count)
}

results <- list()
for (category in names(dirs)) {
  files <- list.files(dirs[[category]], pattern = "\\.bed$", full.names = TRUE)
  for (file in files) {
    biosample <- str_extract(basename(file), "^[^.]+")
    line_count <- count_lines(file)
    results <- append(results, list(data.frame(Biosample = biosample, Category = category, Count = line_count)))
  }
}

df <- bind_rows(results) %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)  # Fill missing values with 0

data <- as.data.frame(df)
data$Unchanged <- data$BackgroundRegions - data$UpInModern - data$UpInArchaic
data$PropUnchanged <- data$Unchanged / data$BackgroundRegions
data$PropUpInModern <- data$UpInModern / data$BackgroundRegions
data$PropUpInArchaic <- data$UpInArchaic / data$BackgroundRegions
data$TotalPropChanged <- data$PropUpInModern + data$PropUpInArchaic
data$Biosample <- factor(data$Biosample, levels = data$Biosample[order(data$TotalPropChanged, decreasing = TRUE)])

plot_data <- data %>%
  mutate(UpInArchaic = -UpInArchaic)
plot_data <- plot_data %>%
  arrange(desc(abs(UpInModern) + abs(UpInArchaic)))
pdf("pdf/modernArchaicDACounts.pdf", height=6, width=8)
ggplot(subset(plot_data, abs(UpInModern) + abs(UpInArchaic) > 0), aes(y = reorder(Biosample, UpInModern + abs(UpInArchaic)))) +
  geom_bar(aes(x = UpInModern), stat = "identity", fill = "#F4A361") + 
  geom_bar(aes(x = UpInArchaic), stat = "identity", fill = "#EAC365") +
  theme_classic() +
  labs(x = "", y = "",
       title = "Number of Differentially Accessible Regions between Modern Humans and Archaic Hominins across biosamples")
dev.off()














