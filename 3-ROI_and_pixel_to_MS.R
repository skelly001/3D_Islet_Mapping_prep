library(tidyverse)
library(RImageJROI)
library(sf)
library(maptools)
library(readxl)


# Import data -------------------------------------------------------------

load("RD2-Slide61_roi_w_edges_and_pixel_names.RData")

p_data <- read_excel(path = "../../../Protein_Data/MS_output/processed/3D_mapping_no_imputation.xlsx",
                   sheet = "noinputation_combined_r", 
                   na = "NaN")

alpha_cell_count <- read_excel(path = "../../../Protein_Data/3D_mapping_IF_summary.xlsx",
                               sheet = "alpha") #James count

beta_cell_count <- read_excel(path = "../../../Protein_Data/3D_mapping_IF_summary.xlsx",
                               sheet = "beta") #James count

acinar_cell_count <- read_excel(path = "../../../Protein_Data/3D_mapping_IF_summary.xlsx",
                               sheet = "acinar") #James count



# Prepare MS Data for Link to Cell ROI by Pixel ---------------------------

# Get accession num
p_data <- p_data %>% 
    mutate(PROTID = sub(".*?\\|(.*?)\\|.*", "\\1", PROTID))
    

p_data_2 <- p_data %>% 
    `rownames<-`(., .$PROTID) %>% 
    t() %>% 
    as.data.frame() %>% 
    .[-c(1:2), ] %>% 
    mutate(pixel = rownames(.),
           .before = colnames(.[1])) %>% 
    `rownames<-`(., NULL)



pixels_slide_61_beta <- vector(mode = "list", length = nrow(beta_cell_count))
    
pixel_list <- beta_cell_count %>% 
    pivot_longer(cols = "3":"11",
                 names_to = "pixel_location",
                 values_to = "cell_count") %>% 
    mutate(pixel = paste0(pixel, pixel_location)) %>% 
    # filter(!is.na(cell_count)) %>% 
    select(Slide, pixel) %>% 
    filter(pixel != c("x60E4", "x68A4", "x72A11")) %>% #These pixels have poor protein coverage. Yumi labeled as missing samples
    filter(pixel != "x60E4") #dont know why this isnt filtered out in above line.


alpha_cell_count_long <- alpha_cell_count %>% 
    pivot_longer(cols = "3":"11",
                 names_to = "pixel_location",
                 values_to = "alpha_cell_count") %>% 
    mutate(pixel = paste0(pixel, pixel_location)) %>%
    select(-c(Slide, pixel_location)) %>% 
    filter(pixel != c("x60E4", "x68A4", "x72A11")) %>% #These pixels have poor protein coverage. Yumi labeled as missing samples
    filter(pixel != "x60E4") #dont know why this isnt filtered out in above line.

beta_cell_count_long <- beta_cell_count %>% 
    pivot_longer(cols = "3":"11",
                 names_to = "pixel_location",
                 values_to = "beta_cell_count") %>% 
    mutate(pixel = paste0(pixel, pixel_location)) %>%
    select(-c(Slide, pixel_location)) %>% 
    filter(pixel != c("x60E4", "x68A4", "x72A11")) %>% #These pixels have poor protein coverage. Yumi labeled as missing samples
    filter(pixel != "x60E4") #dont know why this isnt filtered out in above line.

acinar_cell_count_long <- acinar_cell_count %>% 
    pivot_longer(cols = "3":"11",
                 names_to = "pixel_location",
                 values_to = "acinar_cell_count") %>% 
    mutate(pixel = paste0(pixel, pixel_location)) %>%
    select(-c(Slide, pixel_location)) %>% 
    filter(pixel != c("x60E4", "x68A4", "x72A11")) %>% #These pixels have poor protein coverage. Yumi labeled as missing samples
    filter(pixel != "x60E4") #dont know why this isnt filtered out in above line.

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


if(!file.exists("RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData")){
    save(old, file = "RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData")
}

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
if(!file.exists("RD3_2-Slide61_roi_w_edges_pixels_ms_relative.RData")){
    save(pix_to_ms, file = "RD3_2-Slide61_roi_w_edges_pixels_ms_relative.RData")
}

plot(old_log2)
plot(old)
plot(pix_to_ms)


# low <- "grey90"
# high <- "midnightblue"
low <- "midnightblue"
high <- "yellow"

source <-  old_log2
source <- old
source <-  pix_to_ms

## Islet Markers
# Pro-glucagon (P01275)
ggplot(source)+
    geom_sf(aes(fill = P01275))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Pro-glucagon")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Pro-glucagon_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


# Insulin (P01308)
ggplot(source)+
    geom_sf(aes(fill = P01308))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Insulin")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Insulin_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


# Secretogranin-2 (P13521)
ggplot(source)+
    geom_sf(aes(fill = P13521))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Secretogranin-2")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Secretogranin-2_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


## Acinar Markers
# Lithostathine-1-alpha (P05451)
ggplot(source)+
    geom_sf(aes(fill = P05451))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Lithostathine-1-alpha")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Lithostathine-1-alpha_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


# Lithostathine-1-beta (P48304)
ggplot(source)+
    geom_sf(aes(fill = P48304))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Lithostathine-1-beta")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Lithostathine-1-beta_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


# Carboxypeptidase A1 (P15085)
ggplot(source)+
    geom_sf(aes(fill = P15085))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("Carboxypeptidase A1")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./Pixel_Level_Sample_images/Img_Carboxypeptidase_A1_pixel.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")





# ACAA2 (P42765)
t <- ggplot(source)+
    geom_sf(aes(fill = P42765))+
    scale_fill_gradient(
        low = low,
        high = high)+
    ggtitle("ACAA2")+
    # annotate("text", label = "test", x = Inf, y = 1)+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text = element_blank())

g = ggplotGrob(t)
g$layout$clip[g$layout$name == "panel"] = "off"
grid.draw(g)


library(gridExtra)




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



