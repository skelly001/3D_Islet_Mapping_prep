library(readxl)
library(tidyverse)
library(sf)


# Import data -------------------------------------------------------------

# ROI data
(load("output/RD2-ROI_to_pixel/RD2-Slide61_roi_w_edges_and_pixel_names.RData"))

# Protein abundance
prot_data <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_protein_abundance.xlsx",
                   sheet = "noinputation_combined", 
                   na = "NaN")

# Cell counts
alpha_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "alpha")

beta_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "beta")

acinar_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "acinar")



# Prepare MS Data for Link to Cell ROI by Pixel ---------------------------

dir.create("output/RD3-ROI_and_pixel_to_MS")

# Get accession num
prot_data <- prot_data %>% 
    mutate(PROTID = sub(".*?\\|(.*?)\\|.*", "\\1", PROTID)) %>% 
    as.data.frame()
    

prot_data_2 <- prot_data %>% 
    `rownames<-`(., .$PROTID) %>% 
    t() %>% 
    as.data.frame() %>% 
    .[-c(1:2), ] %>% 
    mutate(pixel = rownames(.),
           .before = colnames(.[1])) %>% 
    `rownames<-`(., NULL)

# These pixels have poor protein coverage.
bad_pixels <- c("x60E4", "x68A4", "x72A11")

transform_counts <- function(df, col_name) {
    df %>% 
        pivot_longer(cols = "3":"11", names_to = "pixel_location", values_to = col_name) %>% 
        mutate(pixel = paste0(pixel, pixel_location)) %>%
        select(-c(Slide, pixel_location)) %>% 
        filter(!pixel %in% bad_pixels)
}

cell_counts <- list(
    transform_counts(alpha_cell_count, "alpha_cell_count"),
    transform_counts(beta_cell_count, "beta_cell_count"),
    transform_counts(acinar_cell_count, "acinar_cell_count")
)

prot_data_2 <- reduce(cell_counts, inner_join, by = "pixel", .init = prot_data_2)



# Link MS Data to Cell ROI by Pixel ---------------------------------------
cell_to_pix <- cell_to_pix %>% 
    rename(pixel = pix.MS.name)

# Annotation columns
anno_cols <- c("pixel", "cell_type", "cell.roi", "geometry", 
"acinar_cell_count", "alpha_cell_count", "beta_cell_count") 


pix_to_ms <- full_join(cell_to_pix, prot_data_2, by = "pixel") %>% 
    mutate_at(.vars = vars(!all_of(anno_cols)), .funs = as.double)

# Only slide 60
pix_to_ms <- pix_to_ms %>%
    filter(!is.na(cell.roi))

# Check missingness
missingness <- as.data.frame(pix_to_ms) %>%
    select(!all_of(anno_cols)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x))}))

hist(as.numeric(missingness[1,]), 100)

# remove proteins w/ no slide 60 observations
pix_to_ms <- pix_to_ms %>% 
    select(all_of(anno_cols), where(\(x){mean(is.na(x)) < 1}))


# Undo log2 transform
prot_cols <- setdiff(colnames(pix_to_ms), anno_cols)
pix_to_ms[prot_cols] <- as.data.frame(pix_to_ms) %>% 
    select(all_of(prot_cols)) %>% 
    lapply(\(x){2^x})

# Save
saveRDS(pix_to_ms, file = "output/RD3-ROI_and_pixel_to_MS/RD3-Slide61_roi_w_edges_pixels_ms_absolute.rds")
