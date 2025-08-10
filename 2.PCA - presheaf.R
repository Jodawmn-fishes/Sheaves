# Load necessary libraries
library(dplyr)
library(tidyr)

# Read the original PCA data
pca_data <- read.csv("pca_data.csv")

# Transform to presheaf format by selecting first principal components
presheaf_data <- pca_data %>%
  select(
    species,
    Genomic_PC1,
    Morphological_PC1,
    Geographic_PC1,
  )

# Write the output to a new CSV file
write.csv(presheaf_data, "pca_presheaf.csv", row.names = FALSE)