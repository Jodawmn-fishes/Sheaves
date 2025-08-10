library(readxl)
library(dplyr)
library(tidyr)

# Read all data files
genomics <- read_excel("ligia_genomics.xlsx")
geography <- read_excel("ligia_geography.xlsx")
morphology <- read_excel("ligia_morphology.xlsx")

# Clean and prepare each dataset
prepare_data <- function(df, prefix) {
  df %>%
    select(-starts_with("species")) %>%
    scale() %>%  # Standardize the data
    as.data.frame() %>%
    setNames(paste0(prefix, "_", names(.)))
}

# Perform PCA for each data type using prcomp()
perform_pca <- function(data, prefix) {
  pca_result <- prcomp(data, scale. = TRUE)
  pca_scores <- as.data.frame(pca_result$x)
  colnames(pca_scores) <- paste0(prefix, "_PC", 1:ncol(pca_scores))
  pca_scores
}

# Prepare each dataset
genomic_data <- prepare_data(genomics, "Genomic")
geo_data <- prepare_data(geography %>% select(-latitude, -longitude), "Geographic")
morph_data <- prepare_data(morphology, "Morphological")

# Perform PCA on each dataset
genomic_pca <- perform_pca(genomic_data, "Genomic")
geo_pca <- perform_pca(geo_data, "Geographic")
morph_pca <- perform_pca(morph_data, "Morphological")

# Combine all data into one dataframe
pca_data <- bind_cols(
  species = genomics$species,
  genomic_pca,
  geo_pca,
  morph_pca,
  latitude = geography$latitude,
  longitude = geography$longitude
)

# Save the combined PCA data
write.csv(pca_data, "pca_data.csv", row.names = FALSE)

# Show the structure of the resulting data
str(pca_data)