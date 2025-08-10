# === Complete Intersection Analysis with Relative Overlap ===
library(Rvcg)
library(geometry)
library(plotly)
library(dplyr)

# 1. Calculate volumes with exact intersection (if available) or sampling fallback
calculate_overlaps <- function(sp1_points, sp2_points, sample_size = 10000) {
  # Ensure minimum points
  if (nrow(sp1_points) < 4 || nrow(sp2_points) < 4) return(list(vol1 = 0, vol2 = 0, intersection = 0))
  
  # Function to create watertight mesh
  create_mesh <- function(points) {
    hull <- convhulln(points, options = "Qt")
    mesh <- list(
      vb = t(cbind(points, 1)),
      it = t(hull),
      primitivetype = "triangle"
    )
    class(mesh) <- "mesh3d"
    Rvcg::vcgClean(mesh, silent = TRUE)
  }
  
  # Calculate individual volumes
  vol1 <- suppressWarnings(Rvcg::vcgVolume(create_mesh(sp1_points)))
  vol2 <- suppressWarnings(Rvcg::vcgVolume(create_mesh(sp2_points)))
  
  # Try exact boolean intersection first
  intersection_vol <- tryCatch({
    mesh1 <- create_mesh(sp1_points)
    mesh2 <- create_mesh(sp2_points)
    intersection <- Rvcg::vcgBoolean(mesh1, mesh2, operation = "intersection")
    if (!is.null(intersection$vb)) Rvcg::vcgVolume(intersection) else 0
  }, error = function(e) {
    # Fallback to sampling method
    set.seed(123) # For reproducibility
    samples <- rbind(
      sp1_points[sample(nrow(sp1_points), min(sample_size, nrow(sp1_points))), ],
      sp2_points[sample(nrow(sp2_points), min(sample_size, nrow(sp2_points))), ]
    )
    
    in_sp1 <- geometry::inhulln(convhulln(sp1_points), samples)
    in_sp2 <- geometry::inhulln(convhulln(sp2_points), samples)
    overlap_points <- samples[in_sp1 & in_sp2, ]
    
    if (nrow(overlap_points) < 4) return(0)
    hull <- convhulln(overlap_points, options = "Qt")
    suppressWarnings(Rvcg::vcgVolume(list(
      vb = t(cbind(overlap_points, 1)),
      it = t(hull)
    )))
  })
  
  # Calculate relative overlaps
  rel_overlap_sp1 <- ifelse(vol1 > 0, intersection_vol/vol1 * 100, 0)
  rel_overlap_sp2 <- ifelse(vol2 > 0, intersection_vol/vol2 * 100, 0)
  
  list(
    sp1_volume = vol1,
    sp2_volume = vol2,
    absolute_intersection = intersection_vol,
    sp1_relative = rel_overlap_sp1,
    sp2_relative = rel_overlap_sp2
  )
}

# 2. Enhanced visualization function
visualize_intersection <- function(pca_data, sp1, sp2) {
  sp1_points <- pca_data %>% filter(species == sp1) %>% 
    select(Genomic_PC1, Morphological_PC1, Geographic_PC1) %>% as.matrix()
  sp2_points <- pca_data %>% filter(species == sp2) %>% 
    select(Genomic_PC1, Morphological_PC1, Geographic_PC1) %>% as.matrix()
  
  hull1 <- convhulln(sp1_points, "Qt")
  hull2 <- convhulln(sp2_points, "Qt")
  
  # Find exact overlapping points
  in_sp1 <- geometry::inhulln(hull1, sp2_points)
  in_sp2 <- geometry::inhulln(hull2, sp1_points)
  overlap_sp2_in_sp1 <- sp2_points[in_sp1, ]
  overlap_sp1_in_sp2 <- sp1_points[in_sp2, ]
  
  # Calculate volumes and overlaps
  overlaps <- calculate_overlaps(sp1_points, sp2_points)
  
  plot_ly() %>%
    add_trace(
      type = 'mesh3d',
      x = sp1_points[,1], y = sp1_points[,2], z = sp1_points[,3],
      i = hull1[,1]-1, j = hull1[,2]-1, k = hull1[,3]-1,
      color = I('blue'), opacity = 0.2, name = paste(sp1, "volume:", round(overlaps$sp1_volume, 3))
    ) %>%
    add_trace(
      type = 'mesh3d',
      x = sp2_points[,1], y = sp2_points[,2], z = sp2_points[,3],
      i = hull2[,1]-1, j = hull2[,2]-1, k = hull2[,3]-1,
      color = I('red'), opacity = 0.2, name = paste(sp2, "volume:", round(overlaps$sp2_volume, 3))
    ) %>%
    add_trace(
      type = 'scatter3d',
      x = overlap_sp2_in_sp1[,1], y = overlap_sp2_in_sp1[,2], z = overlap_sp2_in_sp1[,3],
      mode = 'markers',
      marker = list(size = 5, color = 'green'),
      name = paste(round(overlaps$sp1_relative, 1), "% of", sp1)
    ) %>%
    add_trace(
      type = 'scatter3d',
      x = overlap_sp1_in_sp2[,1], y = overlap_sp1_in_sp2[,2], z = overlap_sp1_in_sp2[,3],
      mode = 'markers',
      marker = list(size = 5, color = 'yellow'),
      name = paste(round(overlaps$sp2_relative, 1), "% of", sp2)
    ) %>%
    layout(
      title = paste("Intersection: ", sp1, "vs", sp2, 
                   "<br>Absolute Overlap Volume:", round(overlaps$absolute_intersection, 3)),
      scene = list(
        xaxis = list(title = "Genomic PC1"),
        yaxis = list(title = "Morphological PC1"),
        zaxis = list(title = "Geographic PC1")
      )
    )
}

# 3. Main execution
species_pairs <- combn(unique(pca_data$species), 2, simplify = FALSE)

if (length(species_pairs) > 0) {
  results <- lapply(species_pairs, function(pair) {
    sp1 <- pair[1]
    sp2 <- pair[2]
    
    # Generate visualization
    p <- visualize_intersection(pca_data, sp1, sp2)
    
    # Calculate numerical results
    sp1_points <- pca_data %>% filter(species == sp1) %>% 
      select(Genomic_PC1, Morphological_PC1, Geographic_PC1) %>% as.matrix()
    sp2_points <- pca_data %>% filter(species == sp2) %>% 
      select(Genomic_PC1, Morphological_PC1, Geographic_PC1) %>% as.matrix()
    
    overlaps <- calculate_overlaps(sp1_points, sp2_points)
    
    # Print results
    cat(sprintf(
      "\n=== %s vs %s ===
      Absolute Volume: %.4f
      %s Relative Overlap: %.1f%%
      %s Relative Overlap: %.1f%%\n",
      sp1, sp2,
      overlaps$absolute_intersection,
      sp1, overlaps$sp1_relative,
      sp2, overlaps$sp2_relative
    ))
    
    # Show plot
    print(p)
    
    overlaps
  })
}