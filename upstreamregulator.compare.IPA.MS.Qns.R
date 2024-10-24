## this script is for comparing cannonical pathway files from MS and Qns cohort
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/IPA/padj<0.05")
# Step 1: Load your datasets
MS <- read.xlsx("upstream.regulator.MS.padj<0.05.xlsx", sheet = 2)
QNS <- read.xlsx("upstream.regulator_QNS_padj<0.05.xlsx", sheet = 2)

# Step 2: Create protein count column for both datasets (MS and QNS)
MS <- MS %>%
  mutate(protein_count = sapply(strsplit(Target.Molecules.in.Dataset, ","), length))  # Count genes in MS dataset

QNS <- QNS %>%
  mutate(protein_count = sapply(strsplit(Target.Molecules.in.Dataset, ","), length))  # Count genes in QNS dataset

# Step 3: Filter the MS dataset based on:
# - p-value.of.overlap <= 0.05
# - Target.Molecules.in.Dataset >= 3
# - Predicted.Activation.State is not blank, NA, or whitespace
MS_filtered <- MS %>%
  filter(`p-value.of.overlap` <= 0.05,                          # Remove rows where p-value.of.overlap > 0.05
         Target.Molecules.in.Dataset >= 3,                     # Remove rows where Target.Molecules.in.Dataset < 3
         !is.na(Predicted.Activation.State),                    # Remove rows where Predicted.Activation.State is NA
         trimws(Predicted.Activation.State) != "")              # Remove rows where Predicted.Activation.State is blank or whitespace

# Step 4: Merge the MS and QNS datasets by 'Upstream.Regulator'
combined_data <- inner_join(MS_filtered, QNS, by = "Upstream.Regulator")

# Step to ensure the columns are numeric
combined_data <- combined_data %>%
  mutate(`Activation.z-score.x` = as.numeric(`Activation.z-score.x`),
         `Activation.z-score.y` = as.numeric(`Activation.z-score.y`))
# If there are warnings about NAs introduced, you can replace them with a value or handle them as needed

# Step to create absolute z-score columns for both MS and QNS
combined_data <- combined_data %>%
  mutate(abs_z_score_MS = abs(`Activation.z-score.x`),   # Absolute value of z-score for MS
         abs_z_score_QNS = abs(`Activation.z-score.y`))  # Absolute value of z-score for QNS
# If there are warnings about NAs introduced, you can replace them with a value or handle them as needed

# Step 6: Sort the combined data by protein count (descending) and absolute z-score (descending) of MS
sorted_combined_data <- combined_data %>%
  arrange(desc(protein_count.x), desc(abs_z_score_MS))
#save this data file for further analysis

write.xlsx(sorted_combined_data,"sorted_consolidated_comparive.upstream.regulator.MSQns.xlsx")

##after saving this result I also wanted to see how many of them are in same direction
# Step 7: Filter to keep only upstream regulators with the same z-score direction (both positive or both negative)
same_direction_data <- combined_data %>%
  filter((`Activation.z-score.x` > 0 & `Activation.z-score.y` > 0) | 
           (`Activation.z-score.x` < 0 & `Activation.z-score.y` < 0))
# Step 8: Sort the filtered data by protein count (descending) and absolute z-score (descending) for MS only
sorted_same_direction_data <- same_direction_data %>%
  arrange(desc(protein_count.x), desc(abs_z_score_MS))
# Step 8: Sort the filtered data by protein count (descending) and absolute z-score (descending) for MS and Qns both
sorted_same_direction_data <- same_direction_data %>%
  arrange(desc(protein_count.x), desc(abs_z_score_MS), desc(abs_z_score_QNS))


#save this new file too
write.xlsx(sorted_same_direction_data, "sorted_same_direction_data.xlsx")
