# Load required libraries
library(dplyr)
library(geometry)
library(Rvcg)
library(rgl)

# === STEP 1: Load PCA Data ===
pca_df <- read.csv("pca_presheaf_synthetic.csv", stringsAsFactors = FALSE)

# === STEP 2: Build Convex Hulls by Species ===
species_meshes <- list()
species_groups <- split(pca_df, pca_df$species)

for (species in names(species_groups)) {
  group <- species_groups[[species]]
  coords <- as.matrix(group[, c("Genomic_PC1", "Morphological_PC1", "Geographic_PC1")])
  
  if (nrow(coords) >= 4) {
    # Compute convex hull
    hull <- convhulln(coords, options = "QJ")
    
    # Create proper mesh3d object
    mesh <- list(
      vb = t(cbind(coords, 1)),  # Homogeneous coordinates
      it = t(hull)               # Triangle indices
    )
    class(mesh) <- "mesh3d"
    
    # Add normals for proper mesh structure
    mesh <- vcgUpdateNormals(mesh)
    
    species_meshes[[species]] <- list(
      mesh = mesh,
      points = group,
      hull = hull
    )
  }
}

# === STEP 3: Detect 3D Polygon Intersections ===
obstructed_indices <- integer(0)
species_list <- names(species_meshes)

# Custom function to safely check intersections
check_intersection <- function(mesh1, mesh2) {
  # Get vertices from both meshes
  verts1 <- t(mesh1$vb[1:3,])
  verts2 <- t(mesh2$vb[1:3,])
  
  # Determine safe number of samples (not exceeding vertex count)
  sample_num <- min(100, nrow(verts1), nrow(verts2))
  
  if (sample_num < 1) return(FALSE)
  
  # Sample points from mesh1 surface
  pts1 <- verts1[sample(nrow(verts1), sample_num), , drop = FALSE]
  
  # Check if points are inside mesh2
  inside1 <- sapply(1:nrow(pts1), function(i) {
    vcgClostKD(matrix(pts1[i,], ncol = 3), mesh2, sign = TRUE)$quality < 0
  })
  
  # Sample points from mesh2 surface
  pts2 <- verts2[sample(nrow(verts2), sample_num), , drop = FALSE]
  
  # Check if points are inside mesh1
  inside2 <- sapply(1:nrow(pts2), function(i) {
    vcgClostKD(matrix(pts2[i,], ncol = 3), mesh1, sign = TRUE)$quality < 0
  })
  
  any(inside1) || any(inside2)
}

for (i in seq_along(species_list)) {
  for (j in seq_along(species_list)) {
    if (i >= j) next  # Skip duplicate pairs and self-checks
    
    sp1 <- species_list[i]
    sp2 <- species_list[j]
    
    mesh1 <- species_meshes[[sp1]]$mesh
    mesh2 <- species_meshes[[sp2]]$mesh
    
    # Check bounding box overlap first for efficiency
    bb1 <- apply(mesh1$vb[1:3,], 1, range)
    bb2 <- apply(mesh2$vb[1:3,], 1, range)
    
    if (all(bb1[1,] <= bb2[2,]) && all(bb2[1,] <= bb1[2,])) {
      # Check for intersection using our safe function
      if (check_intersection(mesh1, mesh2)) {
        obstructed_indices <- union(obstructed_indices, as.numeric(rownames(species_meshes[[sp1]]$points)))
        obstructed_indices <- union(obstructed_indices, as.numeric(rownames(species_meshes[[sp2]]$points)))
      }
    }
  }
}

# === STEP 4: Tag and Export Obstructions ===
pca_df$refined_obstructed <- ifelse(1:nrow(pca_df) %in% obstructed_indices, "Yes", "No")

# Save to CSV
output_file <- "obstructions_synthetic.csv"
write.csv(pca_df, output_file, row.names = FALSE)