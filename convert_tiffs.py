#!/usr/bin/env Rscript
# NOTE: Despite the .py extension from the original plan, this script
# was implemented in R since R + magick is available on this system.
# Run with: Rscript convert_tiffs.py
#
# Converts large TIFF microscopy images to web-optimized JPEGs.
# Also crops slider comparison images to remove ggplot legends.

library(magick)

root <- getwd()
input_dir <- file.path(root, "data", "1-ImageJ_cell_ROI")
output_dir <- file.path(root, "website", "converted")
protein_map_dir <- file.path(root, "output", "RD5-final_protein_maps", "final_protein_maps")
cell_type_map <- file.path(root, "output", "RD1-ROI_mapping_and_cell_type_assignment", "cell_type_assignment.png")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- TIFF to JPEG Conversion ---
cat("\n=== TIFF to JPEG Conversion ===\n\n")

tiff_conversions <- list(
  list(src = "Slide61_Max_Projection_DAPI.tif", dst = "DAPI_fluorescence.jpg"),
  list(src = "Result of Slide61_Max_Projection_INS_GCG_DAPI.tif", dst = "INS_GCG_DAPI_composite.jpg"),
  list(src = "cropped.tif", dst = "tissue_cropped.jpg"),
  list(src = "Slide61-islet23_with_cell_boundaries_with_edges.tif", dst = "cell_boundaries_overlay.jpg")
)

for (conv in tiff_conversions) {
  src_path <- file.path(input_dir, conv$src)
  dst_path <- file.path(output_dir, conv$dst)

  if (!file.exists(src_path)) {
    cat("  SKIP:", conv$src, "not found\n")
    next
  }

  cat("  Converting:", conv$src, "\n")
  img <- image_read(src_path)
  info <- image_info(img)

  # Resize if longest edge > 2000px
  max_dim <- max(info$width, info$height)
  if (max_dim > 2000) {
    scale_pct <- round(2000 / max_dim * 100)
    img <- image_resize(img, paste0(scale_pct, "%"))
  }

  # Normalize to 8-bit if needed (magick handles this automatically)
  img <- image_normalize(img)

  image_write(img, path = dst_path, format = "jpeg", quality = 85)
  size_kb <- round(file.size(dst_path) / 1024)
  new_info <- image_info(img)
  cat("    ->", conv$dst, paste0("(", size_kb, " KB, ", new_info$width, "x", new_info$height, ")\n"))
}

# --- Slider Image Crops ---
cat("\n=== Slider Image Crops ===\n\n")

slider_crops <- list(
  list(src = cell_type_map, dst = "cell_type_cropped.png"),
  list(src = file.path(protein_map_dir, "INS_P01308.png"), dst = "INS_cropped.png"),
  list(src = file.path(protein_map_dir, "GCG_P01275.png"), dst = "GCG_cropped.png")
)

for (crop in slider_crops) {
  if (!file.exists(crop$src)) {
    cat("  SKIP:", basename(crop$src), "not found\n")
    next
  }

  cat("  Cropping:", basename(crop$src), "\n")
  img <- image_read(crop$src)
  info <- image_info(img)

  # Crop to left ~72% of width (remove right-side ggplot legend)
  crop_width <- round(info$width * 0.72)
  crop_geom <- paste0(crop_width, "x", info$height, "+0+0")
  img <- image_crop(img, crop_geom)

  dst_path <- file.path(output_dir, crop$dst)
  image_write(img, path = dst_path)
  size_kb <- round(file.size(dst_path) / 1024)
  new_info <- image_info(img)
  cat("    ->", crop$dst, paste0("(", size_kb, " KB, ", new_info$width, "x", new_info$height, ")\n"))
}

cat("\n=== Done! ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:", length(list.files(output_dir)), "\n")
