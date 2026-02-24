#!/usr/bin/env Rscript
# generate_slider_images.R
# Generates legend-free slider comparison images for the website.
# All three images use the SAME sf data source (pix_to_ms) so the
# polygon geometries and coordinate extents are identical.

library(ggplot2)
library(sf)

out_dir <- "website/converted"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Shared dimensions and theme (no legend, no axes, no margins)
W <- 1100
H <- 1100
DPI <- 300

shared_theme <- theme(
  legend.position = "none",
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = margin(0, 0, 0, 0),
  panel.border = element_blank()
)

# Load the protein data (same source for ALL three images)
cat("Loading protein data...\n")
pix_to_ms <- readRDS("output/RD3-ROI_and_pixel_to_MS/RD3-Slide61_cell_roi_ms_absolute.rds")

# --- 1. Cell Type Map (from pix_to_ms, same geometry as protein maps) ---
cat("Generating cell type map (no legend, from pix_to_ms)...\n")
p_ct <- ggplot(pix_to_ms) +
  geom_sf(aes(fill = cell_type)) +
  scale_fill_brewer(palette = "Dark2") +
  shared_theme

png(file.path(out_dir, "cell_type_cropped.png"), width = W, height = H, res = DPI)
plot(p_ct)
dev.off()
cat("  -> cell_type_cropped.png\n")
rm(p_ct)

# --- 2. INS Protein Map ---
cat("Generating INS map (no legend)...\n")
p_ins <- ggplot(pix_to_ms[, c("P01308", "geometry")]) +
  geom_sf(aes(fill = .data[["P01308"]])) +
  scale_fill_gradient(low = "midnightblue", high = "yellow") +
  shared_theme

png(file.path(out_dir, "INS_cropped.png"), width = W, height = H, res = DPI)
plot(p_ins)
dev.off()
cat("  -> INS_cropped.png\n")
rm(p_ins)

# --- 3. GCG Protein Map ---
cat("Generating GCG map (no legend)...\n")
p_gcg <- ggplot(pix_to_ms[, c("P01275", "geometry")]) +
  geom_sf(aes(fill = .data[["P01275"]])) +
  scale_fill_gradient(low = "midnightblue", high = "yellow") +
  shared_theme

png(file.path(out_dir, "GCG_cropped.png"), width = W, height = H, res = DPI)
plot(p_gcg)
dev.off()
cat("  -> GCG_cropped.png\n")

cat("\nDone! All 3 slider images are", W, "x", H, "px, no legends, same sf source.\n")
