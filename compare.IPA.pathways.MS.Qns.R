# this script is for comparing IPA results of MS and Qns
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/IPA/padj<0.05")
# Step 1: Load your datasets
MS <- read.xlsx("cannonical_pathway_MS_padj.<0.05.xlsx", sheet = 1)
QNS <- read.xlsx("cannonical_pathway_QNS_padj<0.05.xlsx", sheet = 1)

# Step 2: Create protein count column for both datasets (MS and QNS)
MS <- MS %>%
  mutate(protein_count = sapply(strsplit(Molecules, ","), length))  # Count genes in MS dataset

QNS <- QNS %>%
  mutate(protein_count = sapply(strsplit(Molecules, ","), length))  # Count genes in QNS dataset

# Step 3: Remove rows where -log(padj) > 1.3 in the MS dataset
MS_filtered <- MS %>%
  filter(`-log(p-value)` >= 1.3)

# Step 4: Merge the MS and QNS datasets by 'Ingenuity.Canonical.Pathways'
combined_data <- inner_join(MS_filtered, QNS, by = "Ingenuity.Canonical.Pathways")

# Step 5: Function to calculate similarity between two gene lists
similarity_percentage <- function(list1, list2) {
  genes1 <- unlist(strsplit(list1, ","))
  genes2 <- unlist(strsplit(list2, ","))
  common_genes <- intersect(genes1, genes2)
  total_genes <- union(genes1, genes2)
  
  similarity <- length(common_genes) / length(total_genes)
  return(similarity)
}

# Step 6: Merge rows with more than 70% similarity in the 'Molecules.x' column
i <- 1
while (i < nrow(combined_data)) {
  j <- i + 1
  while (j <= nrow(combined_data)) {
    list1 <- combined_data$Molecules.x[i]
    list2 <- combined_data$Molecules.x[j]
    
    similarity <- similarity_percentage(list1, list2)
    
    if (similarity > 0.7) {
      # Merge gene lists
      combined_genes <- union(unlist(strsplit(list1, ",")), unlist(strsplit(list2, ",")))
      combined_data$Molecules.x[i] <- paste(combined_genes, collapse = ",")
      
      # Remove the second row
      combined_data <- combined_data[-j, ]
    } else {
      j <- j + 1
    }
  }
  i <- i + 1
}

# Step 7: Sort combined data by protein count and absolute z-score
combined_data <- combined_data %>%
  mutate(abs_z_score_MS = abs(`z-score.x`),  # Absolute value of z-score for MS
         abs_z_score_QNS = abs(`z-score.y`)) %>%
  arrange(desc(protein_count.x), desc(abs_z_score_MS))

#write the combined data file 
write.xlsx(combined_data, "MS.QNS.compare.merged70%similarities.protein.xlsx")
# Step 8: Reshape the combined dataset into a long format for plotting
comparison_long <- combined_data %>%
  select(Ingenuity.Canonical.Pathways, `-log(p-value).x`, `z-score.x`, protein_count.x,
         `-log(p-value).y`, `z-score.y`, protein_count.y) %>%
  pivot_longer(cols = c(`-log(p-value).x`, `z-score.x`, protein_count.x,
                        `-log(p-value).y`, `z-score.y`, protein_count.y),
               names_to = c(".value", "dataset"),
               names_pattern = "(.*)\\.(x|y)") %>%
  mutate(dataset = ifelse(dataset == "x", "MS", "QNS"))

# Step 9: Get the min and max values of z-score to set axis range and gradient color scale
z_min <- min(comparison_long$`z-score`, na.rm = TRUE)
z_max <- max(comparison_long$`z-score`, na.rm = TRUE)

# Step 10: Create the bubble plot
ggplot(comparison_long[1:40,], aes(x = `z-score`, y = Ingenuity.Canonical.Pathways)) +
  # Add bubbles where size is based on protein count and color is based on -log(p-value)
  geom_point(aes(size = protein_count, color = `-log(p-value)`)) +
  
  # Customize bubble sizes for protein counts (adjust the range as needed)
  scale_size_continuous(range = c(3, 15), name = "Protein Count") +
  
  # Use a gradient color scale for p-values (adjust the colors as needed)
  scale_color_gradient(low = "blue", high = "red", name = "-log(p-value)") +
  
  # Facet the plot based on the dataset (MS vs QNS)
  facet_wrap(~dataset) +
  
  # Customize the axis labels and remove plot title or y-axis title
  labs(x = "z-score", y = NULL) +
  
  # Customize the theme
  theme_classic(base_size = 10) +
  
  # Customize the font for Ingenuity.Canonical
  # Customize the font for Ingenuity.Canonical.Pathways (term descriptions)
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black"),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10), size = rel(1.1)),
        legend.position = "right",  # Position the legend on the right
        legend.direction = "vertical") +
  
  # Customize the x-axis text angle (optional if terms are too long)
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Customize the x-axis breaks for z-score (e.g., 0, 0.25, 0.5, 0.75, 1)
  scale_x_continuous(breaks = seq(z_min, z_max, by = 2))
##step 8 to 10 is for data visualization side by side
## below given optional steps for getting the pathways are changed in same direction
# New Step: Identify pathways with the same z-score direction (both positive or both negative)
combined_data <- combined_data %>%
  mutate(same_direction = ifelse(sign(`z-score.x`) == sign(`z-score.y`), "Same Direction", "Opposite Direction"))

# Filter pathways that are in the same direction
same_direction_pathways <- combined_data %>%
  filter(same_direction == "Same Direction")

# Optionally, view the pathways in the same direction
head(same_direction_pathways)

#write the pathways in same direction
write.xlsx(same_direction_pathways, "same.direction.IPA.pathways.MS.Qns.xlsx")
