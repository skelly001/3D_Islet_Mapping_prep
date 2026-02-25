library(magick)

cat("Reading Slide61_mask_8bit_color.png...\n")
img <- image_read("data/2-ROI_to_pixel/Slide61_mask_8bit_color.png")
info <- image_info(img)
cat("  Original:", info$width, "x", info$height, "\n")

# Resize only — no normalization or color changes
cat("Resizing to 2007x2007...\n")
img_resized <- image_resize(img, "2007x2007!")
rm(img)
gc()

image_write(img_resized, "website/converted/Slide61_mask_8bit_color.jpg", format = "jpeg", quality = 95)
fsize <- file.size("website/converted/Slide61_mask_8bit_color.jpg")
cat("  File size:", round(fsize / 1024), "KB\n")
cat("Done\n")
