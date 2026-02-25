library(magick)

out_dir <- "website/converted"

# --- Process cell type assignment ---
cat("Processing cell type map...\n")
ct <- image_read("output/RD1-ROI_mapping_and_cell_type_assignment/cell_type_assignment.png")
ct_info <- image_info(ct)
cat("  Original:", ct_info$width, "x", ct_info$height, "\n")

# Crop to remove the ggplot legend (right ~30%) and trim margins
# First remove right side legend
ct_crop <- image_crop(ct, paste0(round(ct_info$width * 0.70), "x", ct_info$height, "+0+0"))
# Then trim whitespace/light gray borders
ct_trim <- image_trim(ct_crop, fuzz = 10)
ct_info2 <- image_info(ct_trim)
cat("  After crop+trim:", ct_info2$width, "x", ct_info2$height, "\n")
image_write(ct_trim, file.path(out_dir, "cell_type_cropped.png"))
rm(ct, ct_crop, ct_trim)
gc()

# --- Process INS protein map ---
cat("Processing INS map...\n")
ins <- image_read("output/RD5-final_protein_maps/final_protein_maps/INS_P01308.png")
ins_info <- image_info(ins)
cat("  Original:", ins_info$width, "x", ins_info$height, "\n")

# Same crop: remove legend, then trim
ins_crop <- image_crop(ins, paste0(round(ins_info$width * 0.70), "x", ins_info$height, "+0+0"))
ins_trim <- image_trim(ins_crop, fuzz = 10)
ins_info2 <- image_info(ins_trim)
cat("  After crop+trim:", ins_info2$width, "x", ins_info2$height, "\n")
image_write(ins_trim, file.path(out_dir, "INS_cropped.png"))
rm(ins, ins_crop, ins_trim)
gc()

# --- Process GCG protein map ---
cat("Processing GCG map...\n")
gcg <- image_read("output/RD5-final_protein_maps/final_protein_maps/GCG_P01275.png")
gcg_info <- image_info(gcg)
cat("  Original:", gcg_info$width, "x", gcg_info$height, "\n")

gcg_crop <- image_crop(gcg, paste0(round(gcg_info$width * 0.70), "x", gcg_info$height, "+0+0"))
gcg_trim <- image_trim(gcg_crop, fuzz = 10)
gcg_info2 <- image_info(gcg_trim)
cat("  After crop+trim:", gcg_info2$width, "x", gcg_info2$height, "\n")
image_write(gcg_trim, file.path(out_dir, "GCG_cropped.png"))
rm(gcg, gcg_crop, gcg_trim)
gc()

cat("\nDone! Check the trimmed dimensions above.\n")
cat("If they differ, we need to resize them to match.\n")
