library(magick)

# Use the FULL-FIELD CellPose mask (4095x4096) that matches the original DAPI TIFF
cat("Reading full-field CellPose mask...\n")
mask <- image_read("data/1-ImageJ_cell_ROI/Cellpose 2.0/Slide61_Max_Projection_DAPI_cellpose_mask.png")
info <- image_info(mask)
cat("  Original:", info$width, "x", info$height, "\n")

# Normalize to make cell labels visible (raw mask is near-black integer labels)
cat("Normalizing...\n")
mask <- image_normalize(mask)

# Boost contrast so background is truly black, cells are bright
mask <- image_level(mask, black_point = 5, white_point = 80)

# Apply cyan tint: convert to grayscale then tint
# First make it a clean grayscale
mask <- image_convert(mask, colorspace = "gray")
# Then tint cyan by compositing with a cyan image
cyan_layer <- image_blank(image_info(mask)$width, image_info(mask)$height, color = "#00d4ff")
mask <- image_composite(cyan_layer, mask, operator = "Multiply")
rm(cyan_layer)

# Resize to match DAPI_fluorescence.jpg (2007x2007)
cat("Resizing to 2007x2007...\n")
mask_resized <- image_resize(mask, "2007x2007!")
rm(mask)
gc()

image_write(mask_resized, "website/converted/DAPI_mask_resized.png")
info2 <- image_info(mask_resized)
cat("  Resized:", info2$width, "x", info2$height, "\n")
cat("Done\n")
