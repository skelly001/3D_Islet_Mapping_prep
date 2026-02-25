library(magick)

cat("Converting DAPI - straight conversion, no normalization, no boost...\n")
img <- image_read("data/1-ImageJ_cell_ROI/Slide61_Max_Projection_DAPI.tif")

# Resize to web-friendly dimensions
info <- image_info(img)
max_dim <- max(info$width, info$height)
if (max_dim > 2000) {
  scale_pct <- round(2000 / max_dim * 100)
  img <- image_resize(img, paste0(scale_pct, "%"))
}

# Straight conversion to sRGB JPEG - no normalize, no level adjust
img <- image_convert(img, colorspace = "sRGB")

image_write(img, "website/converted/DAPI_fluorescence.jpg", format = "jpeg", quality = 90)
size_kb <- round(file.size("website/converted/DAPI_fluorescence.jpg") / 1024)
cat("  -> DAPI_fluorescence.jpg (", size_kb, "KB)\n")
