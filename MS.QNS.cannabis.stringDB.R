###thgis script is made to compare the string result of MS and QNS , sort them for 
# 4> gene count in each category, FDR< 0.05 , and merge rows as per 70% similarities
#load necessary libraries and data
library("data.table")
library("openxlsx")
library("readxl")
library("tidyverse")
library("ggplot2")
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/stringdb/")
MS.controls = read.xlsx("MS_unique_control/string_results_unique_control_MS.xlsx", sheet = 1)
MS.cannabis = read.xlsx("MS_unique_cannabis/stringDB_results_MS_cannabis.xlsx", sheet = 1)
QNS.controls = read.xlsx("QNS_unique_control/string_results_unique_control_Qns.xlsx", sheet = 1)
QNS.cannabis = read.xlsx("QNS_unique_cannabis/string_results_unique_cannabis_Qns.xlsx", sheet = 1)

library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Merge datasets and remove duplicates
cannabis.common <- inner_join(MS.cannabis, QNS.cannabis, by = "term.ID") %>%
  distinct(term.ID, .keep_all = TRUE) %>%
  filter(!(`#category.x` %in% c("GO Component", "COMPARTMENTS", "DISEASES")))

# Step 2: Function to calculate similarity between two gene lists
similarity_percentage <- function(list1, list2) {
  genes1 <- unlist(strsplit(list1, ","))
  genes2 <- unlist(strsplit(list2, ","))
  common_genes <- intersect(genes1, genes2)
  total_genes <- union(genes1, genes2)
  
  similarity <- length(common_genes) / length(total_genes)
  return(similarity)
}

# Step 3: Merge rows with more than 70% similarity
i <- 1
while (i < nrow(cannabis.common)) {
  j <- i + 1
  while (j <= nrow(cannabis.common)) {
    list1 <- cannabis.common$`matching.proteins.in.your.network.(labels).x`[i]
    list2 <- cannabis.common$`matching.proteins.in.your.network.(labels).x`[j]
    
    similarity <- similarity_percentage(list1, list2)
    
    if (similarity > 0.7) {
      combined_genes <- union(unlist(strsplit(list1, ",")), unlist(strsplit(list2, ",")))
      cannabis.common$`matching.proteins.in.your.network.(labels).x`[i] <- paste(combined_genes, collapse = ",")
      cannabis.common <- cannabis.common[-j, ]
    } else {
      j <- j + 1
    }
  }
  i <- i + 1
}

# Step 4: Reshape the data into a long format
cannabis.long <- cannabis.common %>%
  select(term.ID, term.description.x, signal.x, observed.gene.count.x, false.discovery.rate.x,
         term.description.y, signal.y, observed.gene.count.y, false.discovery.rate.y) %>%
  pivot_longer(cols = c(term.description.x, signal.x, observed.gene.count.x, false.discovery.rate.x,
                        term.description.y, signal.y, observed.gene.count.y, false.discovery.rate.y),
               names_to = c(".value", "dataset"),
               names_pattern = "(.*)\\.(x|y)") %>%
  mutate(dataset = ifelse(dataset == "x", "MS.cannabis", "QNS.cannabis"))

# Step 5: Create the bubble plot with facets
p <- ggplot(cannabis.long, aes(x = signal, y = term.description)) +
  geom_point(aes(size = observed.gene.count, color = false.discovery.rate)) +
  scale_size_continuous(range = c(3, 15), 
                        limits = c(min(cannabis.long$observed.gene.count, na.rm = TRUE), 
                                   max(cannabis.long$observed.gene.count, na.rm = TRUE))) +  
  scale_color_gradient(low = "blue", high = "red", 
                       limits = c(min(cannabis.long$false.discovery.rate, na.rm = TRUE), 
                                  max(cannabis.long$false.discovery.rate, na.rm = TRUE)),
                       name = "FDR") +
  facet_wrap(~dataset) +
  labs(x = "Signal", 
       size = "Gene Count",  
       color = paste("FDR (", 
                     round(min(cannabis.long$false.discovery.rate, na.rm = TRUE), 4), " to ", 
                     round(max(cannabis.long$false.discovery.rate, na.rm = TRUE), 4), ")", sep = "")) +  
  theme_classic(base_size = 20) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1.1), color = 'black'),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),  # Customize font for y-axis labels
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14, face = "bold"))

# Print the plot
print(p)
