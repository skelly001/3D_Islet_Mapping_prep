library(readxl)
library(RColorBrewer)
library(tidyverse)
library(sf)

# Import Data -------------------------------------------------------------

# Read segmented cell (sc) ROI (8-bit color) histogram measurements
file_list_sc <- list.files(path = r"(data\2-ROI_to_pixel\1-ImageJ_cell_ROI_to_pixel_histograms)",
                        full.names = T)
sc_hist <- map(file_list_sc, data.table::fread, data.table = F)


# Read pixel (pix) ROI (8-bit color) histogram measurements
file_list_pix <- list.files(path = r"(data\2-ROI_to_pixel\2-ImageJ_pixel_ROI_histograms)",
                        full.names = T)
pix_hist <- map(file_list_pix, data.table::fread, data.table = F)


# Read pixel ROI name to Mass Spec pixel name key
pix_ms_names <- read_xlsx(path = r"(data\2-ROI_to_pixel\Slide61_Pixel_Flu_RoiSet_IJ names_to_pixel_names.xlsx)")

# Load segmented cell ROIs
(load("output/RD1-ROI_mapping_and_cell_type_assignment/ROI_polygons_with_cell_types.RData"))


# Process Segmented Cell Histogram Data -----------------------------------

dir.create("output/RD2-ROI_to_pixel")

# Extract file order from file name
file_list_sc_2 <- file_list_sc %>% 
    as.data.frame() %>% 
    mutate(index = str_extract(`.`, "[[:digit:]]{1,4}.txt")) %>% 
    mutate(index = sub(".txt", "", index)) %>% 
    select(index) %>% 
    map_df(as.integer)

# Add file order column to each file
for (i in 1:nrow(file_list_sc_2)) {
    sc_hist[[i]] <- sc_hist[[i]] %>% 
        mutate(index = rep(file_list_sc_2$index[[i]], 256))
}

# Single dataframe and arrange files in correct order
sc_hist <- bind_rows(sc_hist) %>% 
    group_by(index) %>% 
    arrange(index) %>% 
    ungroup()

# Assign ROI to dominant color bin
sc_hist <- sc_hist %>% 
    group_by(index) %>% 
    filter(Count > 0) %>% 
    filter(Count == max(Count)) %>% 
    ungroup() %>% 
    select(-c(Count, V1)) %>% 
    rename(index.sc = index)
# Note: index.sc will match to segmented cell ROI names


# Process Pixel Histogram Data --------------------------------------------


# Extract file order from file name
file_list_pix_2 <- file_list_pix %>% 
    as.data.frame() %>% 
    mutate(index = str_extract(`.`, "[[:digit:]]{1,2}.txt")) %>% 
    mutate(index = sub(".txt", "", index)) %>% 
    select(index) %>% 
    map_df(as.integer)

# Add file order column to each file
for (i in 1:nrow(file_list_pix_2)) {
    pix_hist[[i]] <- pix_hist[[i]] %>% 
        mutate(index = rep(file_list_pix_2$index[[i]], 256))
}

# Single dataframe and arrange files in correct order (low to high)
pix_hist <- bind_rows(pix_hist) %>% 
    group_by(index) %>% 
    arrange(index) %>% 
    ungroup()

# Assign ROI to dominant color bin
pix_hist <- pix_hist %>% 
    group_by(index) %>% 
    filter(Count > 0) %>% 
    ungroup() %>% 
    select(-c(Count, V1)) %>% 
    rename(index.pix = index)
# Note: all pixels are unique by "Value". index.pix will match to pixel ROI names and 
# pixel MS names.  



# Process Pixel MS Names --------------------------------------------------


pix_ms_names <- pix_ms_names %>% 
    select(-c(Name, Row, Col)) %>% 
    rename(index.pix = Index,
           pix.MS.name = MS.pixel.name)


# Link Cell ROIs to Pixels ------------------------------------------------


# Correlate cell ROIs to pixels by "Value" (aka color bin)
cell_to_pix <- left_join(sc_hist, pix_hist, by = "Value")
# Note: Rows with NA in count.pix & index.pix columns are acinar cells


# Add MS pixel names
cell_to_pix <- left_join(cell_to_pix, pix_ms_names, by = "index.pix") %>% 
    select(-index.pix)

# Make key for segmented cell ROIs
roi_all_sf_polygon <- roi_all_sf_polygon %>% 
    mutate(index.sc = 0:(nrow(roi_all_sf_polygon)-1)) %>% 
    rename(cell.roi = roi)

# Correlate segmented cell ROIs with pixels
cell_to_pix  <- full_join(roi_all_sf_polygon, cell_to_pix, by = "index.sc")%>% 
        select(-c(index.sc, Value)) %>% 
    relocate(pix.MS.name, .before = cell.roi) %>% 
    relocate(cell_type, .before = cell.roi)


save(cell_to_pix, file = "output/RD2-ROI_to_pixel/RD2-Slide61_roi_w_edges_and_pixel_names.RData")


# ROI level information
plot(cell_to_pix)


## Visualize cell ROI pixel assignments -----
pal <- brewer.pal(n = 12, name = "Paired")

ggplot(cell_to_pix)+
    geom_sf(aes(fill = `pix.MS.name`))+
    scale_fill_manual(values = rep(pal, 6))+
    theme_dark()+
    theme(legend.title = element_text("Cell Type"),
          legend.position = "none")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD2-ROI_to_pixel/Pixel Assignment.png", 
dpi = 300, width = 1200, height = 1100, units = "px")
