#this script is for IPA analysis, for incoporating protein counts for pathways and upstream regulators
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)


setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/IPA/MSQNScombined.padj<0.05")      
# Step 1: Load your datasets
data = read.xlsx("final_logFC>1pdj.0.01/human_pathways.logFC+-1padj.0.01.xlsx",sheet = 1)


# Step 2: Create protein count column 
data <- data %>%
  mutate(protein_count = sapply(strsplit(Molecules, ","), length))


##save the data only with protein count
write.csv(data, "Human.DEP.bothsex.q<0.05logFC+-1.proteincount.csv")
# Step 3: Remove rows where -log(padj) > 1.3 in the to keep only significat pathways
data <- data %>%
  filter(`-log(p-value)` > 1.3)




# Step 4: Function to calculate similarity between two gene lists
similarity_percentage <- function(list1, list2) {
  genes1 <- unlist(strsplit(list1, ","))
  genes2 <- unlist(strsplit(list2, ","))
  common_genes <- intersect(genes1, genes2)
  total_genes <- union(genes1, genes2)
  
  similarity <- length(common_genes) / length(total_genes)
  return(similarity)
}


# Step 5: Merge rows with more than 60% similarity
i <- 1
while (i < nrow(data)) {
  j <- i + 1
  while (j <= nrow(data)) {
    list1 <- data$Molecules[i]
    list2 <- data$Molecules[j]
    
    similarity <- similarity_percentage(list1, list2)
    
    if (similarity > 0.6) {
      combined_genes <- union(unlist(strsplit(list1, ",")), unlist(strsplit(list2, ",")))
      data$Molecules[i] <- paste(combined_genes, collapse = ",")
      data <- data[-j, ]
    } else {
      j <- j + 1
    }
  }
  i <- i + 1
}
# Sort by protein count and absolute Z-score
data <- data %>% arrange(desc(protein_count), desc(`z-score`))
#save the table
write.csv(data, "Human.DEP.bothsex.q<0.05logFC+-1.simililarity.merged.csv")


#####for upstream regulators
# Step 1: Load your datasets
data = read.xlsx("final_logFC>1pdj.0.01/upstream.regulator.logFC+-1padj.0.01.xlsx",sheet = 1)


# Step 2: Create protein count column 
data <- data %>%
  mutate(protein_count = sapply(strsplit(Target.Molecules.in.Dataset, ","), length))


##save the data only with protein count
write.csv(data, "IPA.Upstream.regulators.for.Human.DEP.bothsex.q<0.05logFC+-1.proteincount.csv")
# Step 3: Remove all non-significant regulators/rows 
data <- data %>%
  filter(`p-value.of.overlap` < 0.05)




# Step 4: Function to calculate similarity between two gene lists
similarity_percentage <- function(list1, list2) {
  genes1 <- unlist(strsplit(list1, ","))
  genes2 <- unlist(strsplit(list2, ","))
  common_genes <- intersect(genes1, genes2)
  total_genes <- union(genes1, genes2)
  
  similarity <- length(common_genes) / length(total_genes)
  return(similarity)
}


# Step 5: Merge rows with more than 60% similarity in protein count
i <- 1
while (i < nrow(data)) {
  j <- i + 1
  while (j <= nrow(data)) {
    list1 <- data$Target.Molecules.in.Dataset[i]
    list2 <- data$Target.Molecules.in.Dataset[j]
    
    similarity <- similarity_percentage(list1, list2)
    
    if (similarity > 0.6) {
      combined_genes <- union(unlist(strsplit(list1, ",")), unlist(strsplit(list2, ",")))
      data$Target.Molecules.in.Dataset[i] <- paste(combined_genes, collapse = ",")
      data <- data[-j, ]
    } else {
      j <- j + 1
    }
  }
  i <- i + 1
}
# Sort by protein count and absolute Z-score
data <- data %>% arrange(desc(protein_count), desc(`z-score`))
#save the table
write.csv(data, "IPA.upstream.regulators.similarity.merged.Human.DEP.bothsex.q<0.05logFC+-1.simililarity.merged.csv")




#write this data for any further analysis
write.xlsx(combined_data, "MS_QNS_compare.xlsx")
# Step 6: Reshape the combined dataset into a long format for plotting
comparison.long <- combined_data %>%
  select(Ingenuity.Canonical.Pathways, `-log(p-value).x`, `z-score.x`, protein_count.x,
         `-log(p-value).y`, `z-score.y`, protein_count.y) %>%
  pivot_longer(cols = c(`-log(p-value).x`, `z-score.x`, protein_count.x,
                        `-log(p-value).y`, `z-score.y`, protein_count.y),
               names_to = c(".value", "dataset"),
               names_pattern = "(.*)\\.(x|y)") %>%
  mutate(dataset = ifelse(dataset == "x", "MS", "QNS"))


# Step 7: Get the min and max values of z-score for axis range and gradient color scale
z_min <- min(comparison.long$`z-score`, na.rm = TRUE)
z_max <- max(comparison.long$`z-score`, na.rm = TRUE)


# Step 8: Create the bubble plot
ggplot(comparison.long[1:40,], aes(x = `z-score`, y = Ingenuity.Canonical.Pathways)) +
  # Add bubbles where size is based on protein count and color is based on -log(p-value)
  geom_point(aes(size = protein_count, color = `-log(p-value)`)) +
  
  # Customize bubble sizes for protein counts (adjust the range as needed)
  scale_size_continuous(range = c(3, 15), name = "Protein Count") +
  
  # Use a gradient color scale for p-values (adjust the colors as needed)
  scale_color_gradient(low = "blue", high = "red", name = "-log(p-value)") +
  
  # Facet the plot based on the dataset (MS vs QNS)
  facet_wrap(~dataset) +
  
  # Customize the axis labels and remove plot title or y-axis title
  labs(x = "Z-Score", y = NULL) +
  
  # Customize the theme
  theme_classic(base_size = 10) +
  
  # Customize the font for Ingenuity.Canonical.Pathways (term descriptions)
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black"),
        axis.title.x = element_text(face = "bold", margin = margin(t = 10), size = rel(1.1)),
        legend.position = "right",  # Position the legend on the right
        legend.direction = "vertical") +
  
  # Customize the x-axis text angle (optional if terms are too long)
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Customize the x-axis breaks for z-score (e.g., 0, 0.25, 0.5, 0.75, 1)
  scale_x_continuous(breaks = seq(z_min, z_max, by = 0.25))
