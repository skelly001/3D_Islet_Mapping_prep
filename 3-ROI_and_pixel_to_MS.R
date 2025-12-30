library(tidyverse)
# library(RImageJROI)
library(sf)
# library(maptools)
library(readxl)


# Import data -------------------------------------------------------------

(load("output/RD2-ROI_to_pixel/RD2-Slide61_roi_w_edges_and_pixel_names.RData"))

p_data <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_no_imputation.xlsx",
                   sheet = "noinputation_combined_r", 
                   na = "NaN")

alpha_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "alpha") #James count

beta_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "beta") #James count

acinar_cell_count <- read_excel(path = "data/3-ROI_and_pixel_to_MS/3D_mapping_IF_summary.xlsx",
                               sheet = "acinar") #James count



# Prepare MS Data for Link to Cell ROI by Pixel ---------------------------

dir.create("output/RD3-ROI_and_pixel_to_MS")

# Get accession num
p_data <- p_data %>% 
    mutate(PROTID = sub(".*?\\|(.*?)\\|.*", "\\1", PROTID)) %>% 
    as.data.frame()
    

p_data_2 <- p_data %>% 
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
    df |> 
        pivot_longer(cols = "3":"11", names_to = "pixel_location", values_to = col_name) |> 
        mutate(pixel = paste0(pixel, pixel_location)) |>
        select(-c(Slide, pixel_location)) |> 
        filter(!pixel %in% bad_pixels)
}

alpha_cell_count_long <- transform_counts(alpha_cell_count, "alpha_cell_count")
beta_cell_count_long <- transform_counts(beta_cell_count, "beta_cell_count")
acinar_cell_count_long <- transform_counts(acinar_cell_count, "acinar_cell_count")


p_data_2 <- inner_join(alpha_cell_count_long, p_data_2, by = "pixel")
p_data_2 <- inner_join(beta_cell_count_long, p_data_2, by = "pixel")
p_data_2 <- inner_join(acinar_cell_count_long, p_data_2, by = "pixel")



# Link MS Data to Cell ROI by Pixel ---------------------------------------
cell_to_pix <- cell_to_pix %>% 
    rename(pixel = pix.MS.name)

# WHY P00167 IS MISSING FROM UNIMPUTED DATA??
# pix_to_ms <- full_join(cell_to_pix, p_data_2, by = "pixel") %>% 
#     mutate_at(.vars = vars(P00167:ncol(.)), .funs = as.double)

pix_to_ms <- full_join(cell_to_pix, p_data_2, by = "pixel") %>% 
    mutate_at(.vars = vars(A0A024RBG1:ncol(.)), .funs = as.double)

old_log2 <- pix_to_ms

# Only slide 60
pix_to_ms <- pix_to_ms %>%
    filter(!is.na(cell.roi))

# Check missingness
missingness <- as.data.frame(pix_to_ms) %>%
    select(8:ncol(.)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x))}))

hist(as.numeric(missingness[1,]), 100)

# remove proteins w/ no slide 60 observations
pix_to_ms <- pix_to_ms %>% 
    select(colnames(pix_to_ms)[1:8], where(\(x){mean(is.na(x)) < 1}))


# Undo log2 transform
pix_to_ms[8:ncol(pix_to_ms)] <- lapply(as.data.frame(pix_to_ms)[8:ncol(pix_to_ms)], \(x){2^x})

old <- pix_to_ms



save(old, file = "output/RD3-ROI_and_pixel_to_MS/RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData")





# Calculate relative abundance
pix_to_ms <- mutate(pix_to_ms, slide = str_extract(pixel, "^x\\d{2}")) %>% 
    relocate(slide, .before = colnames(.)[1])

na <- pix_to_ms %>% 
    filter(is.na(slide))

not_na <- pix_to_ms %>% 
    filter(!is.na(slide))

not_na_relative <- not_na %>%
    group_by(slide) %>%
    mutate(across(!matches(colnames(pix_to_ms)[1:8]), \(x) {x / max(x, na.rm = T)})) %>%
    ungroup()

# not_na_relative <- not_na %>% 
#     group_by(slide) %>% 
#     mutate(across(!matches(colnames(pix_to_ms)[1:8]), \(x) {(x - min(x, na.rm = T)) / (range(x, na.rm=T)[2] - range(x, na.rm=T)[1])})) %>%
#     ungroup()

pix_to_ms <- bind_rows(not_na_relative, na)


# Check transformation by comparing to transformation performed on a single column
old_m <- mutate(old, slide = str_extract(pixel, "^x\\d{2}")) %>% 
    relocate(slide, .before = colnames(.)[1]) %>% 
    select("slide", "A0A0B4J2D5")

old_m_na <- old_m %>% 
    filter(is.na(slide))

old_m_not_na <- old_m %>% 
    filter(!is.na(slide))

old_m_not_na_relative <- old_m_not_na %>%
    group_by(slide) %>%
    mutate(across(matches("A0A0B4J2D5"), \(x) {x / max(x, na.rm = T)})) %>%
    ungroup()
# old_m_not_na_relative <- old_m_not_na %>% 
#     group_by(slide) %>% 
#     mutate(across(matches("A0A0B4J2D5"), \(x) {(x - min(x, na.rm = T)) / (range(x, na.rm=T)[2] - range(x, na.rm=T)[1])})) %>%
#     ungroup()

old_m_pix_to_ms <- bind_rows(old_m_not_na_relative, old_m_na)

new_pix_to_ms <- pix_to_ms[c("slide", "A0A0B4J2D5")]

identical(new_pix_to_ms, old_m_pix_to_ms)

# save
# save(pix_to_ms, file = "RD3_2-Slide61_roi_w_edges_pixels_ms_relative.RData")




plot(old_log2)
plot(old)
plot(pix_to_ms)


# low <- "grey90"
# high <- "midnightblue"
low <- "midnightblue"
high <- "yellow"

# source <-  old_log2
# source <- old
source <-  pix_to_ms

## Islet Markers -----
# Glucagon (P01275)
ggplot(source)+
    geom_sf(aes(fill = P01275))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Glucagon")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Glucagon_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")


# Insulin (P01308)
ggplot(source)+
    geom_sf(aes(fill = P01308))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Insulin")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Insulin_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")


# Secretogranin-2 (P13521)
ggplot(source)+
    geom_sf(aes(fill = P13521))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Secretogranin-2")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Secretogranin-2_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")


## Acinar Markers -----
# Lithostathine-1-alpha (P05451)
ggplot(source)+
    geom_sf(aes(fill = P05451))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Lithostathine-1-alpha")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Lithostathine-1-alpha_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")


# Lithostathine-1-beta (P48304)
ggplot(source)+
    geom_sf(aes(fill = P48304))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Lithostathine-1-beta")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Lithostathine-1-beta_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")


# Carboxypeptidase A1 (P15085)
ggplot(source)+
    geom_sf(aes(fill = P15085))+
    scale_fill_gradient(
        low = low,
        high = high)+
    # ggtitle("Carboxypeptidase A1")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD3-ROI_and_pixel_to_MS/Carboxypeptidase_A1_pixel.png", 
dpi = 300, width = 1570, height = 1100, units = "px")





# # ACAA2 (P42765)
# t <- ggplot(source)+
#     geom_sf(aes(fill = P42765))+
#     scale_fill_gradient(
#         low = low,
#         high = high)+
#     # ggtitle("ACAA2")+
#     # annotate("text", label = "test", x = Inf, y = 1)+
#     theme(plot.title = element_text(hjust = 0.5))+
#     labs(fill = "Relative Intensity")+
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           axis.text = element_blank())

# g = ggplotGrob(t)
# g$layout$clip[g$layout$name == "panel"] = "off"
# grid.draw(g)


# library(gridExtra)




# walk(colnames(pix_to_ms[, 8:10]), function(x){
#     ggplot(pix_to_ms[, x])+
#         geom_sf(aes(fill = x))+
#         scale_fill_gradient(
#             low = "#F5F5F5",
#             high = "#131F56")
# })



# Future work:
# 0. filter MSF pixel data by slide (need 4 total, each with filtered MS data 
#    & new shape files)
# 1. make df for cells/pixel (James C) #DONE
# 2. assign cell type (e.g., 0001-0040 = beta)
# 3. assign cells to pixels
# 4. distribute pixel protein intensity by RNAseq and # cells/pixel (from James C?)
#    mult intensity by RNAseq cell ratio and divide by # cells per
#    cell type in each pixel
# 5. plot map over islet image(?) 
# 6. undo log transform?

# Note: current test image data from slide 61



