# ==============================================================================
# ROI Mapping and Cell Type Assignment from Immunofluorescence
# ==============================================================================
# Script: 1b_1-ROI_mapping_and_cell_type_assignment.R
# Description: Processes ImageJ ROI data to create cell boundary polygons and
#              assigns cell types (alpha, beta, acinar) based on RGB
#              fluorescence intensity measurements. Converts ROI coordinates from 
#              ImageJ (Fiji) cell segmentation to sf polygon objects and visualizes 
#              cell segmentation results.
#
# Input: - data/1-ImageJ_cell_ROI/ROI_RGB_values/*.csv
#        - data/1-ImageJ_cell_ROI/ROI.zip
# Output: - output/RD1-ROI_mapping_and_cell_type_assignment/ROI_polygons_with_cell_types.RData
#         - output/RD1-ROI_mapping_and_cell_type_assignment/segmented_cells.png
#         - output/RD1-ROI_mapping_and_cell_type_assignment/cell_type_assignment.png
# ==============================================================================

library(RColorBrewer)
library(RImageJROI)
library(tidyverse)
library(sf)

# Import ------------------------------------------------------------------


# Read ROI RGB color measurements
file.list <- list.files(path = r"(data\1-ImageJ_cell_ROI\ROI_RGB_values)",
                        full.names = T)
data <- map(file.list, data.table::fread, data.table = F)


# Read ROI from mask
roi_all <- read.ijzip(file = r"(data\1-ImageJ_cell_ROI\ROI.zip)", 
                      verbose = F)



# Process ----------------------------------------------------------------

dir.create("output/RD1-ROI_mapping_and_cell_type_assignment")



# Create ROI polygons ----------------------------------------------------


# Extract coordinates to long format
roi_all_coords <- vector("list", length(roi_all))

for (i in 1:length(roi_all)) {
    roi_all_coords[[i]] <- as.data.frame(roi_all[[i]][["coords"]]) %>%
        mutate(roi = names(roi_all[i]))
}

roi_all_coords <- bind_rows(roi_all_coords)

# Flip y-axis to match image orientation
roi_all_coords[, "y"] <- roi_all_coords[, "y"] * -1

# Convert to sf and create convex hull polygons
roi_all_sf <- st_as_sf(roi_all_coords, coords = c("x","y"))

# Points to polygon
roi_all_sf_polygon <- roi_all_sf %>% 
    group_by(roi) %>%
    summarise() %>% 
    st_cast("POLYGON") %>% 
    st_convex_hull()

# Visualize ROIs
ggplot(roi_all_sf_polygon)+
    geom_sf(aes(fill = roi))+
    scale_fill_manual(values = rep(brewer.pal(12, "Set3"), 
		length.out = length(unique(roi_all_sf_polygon$roi))))+
    theme(legend.position = "none")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD1-ROI_mapping_and_cell_type_assignment/segmented_cells.png",
dpi = 300, width = 1200, height = 1100, units = "px")



# Associate ROI with color ------------------------------------------------

# Parse file indices from filenames
file.list.2 <- file.list %>% 
    as.data.frame() %>% 
    mutate(index = str_extract(`.`, "[[:digit:]]{1,4}.csv")) %>% 
    mutate(index = sub(".csv", "", index)) %>% 
    select(index) %>% 
    map_df(as.integer)

# Add index column to each file
for (i in 1:nrow(file.list.2)) {
    data[[i]] <- data[[i]] %>%
        mutate(index = rep(file.list.2$index[[i]], 3))
}

# Combine and sort by index
data <- bind_rows(data) %>% 
    group_by(index) %>% 
    arrange(index) %>% 
    ungroup()

# Assign cell types based on RGB channel intensities
data_2 <- data %>% 
    group_by(index) %>% 
    mutate(grouped_mean = mean(mean),
           grouped_max = max(mean)) %>% 
    arrange(mean, .by_group = TRUE) %>% 
    mutate(cell_type = case_when(near(grouped_mean[which.max(mean)], grouped_max[which.max(mean)], 2) ~ "acinar",
                                 channel[which.max(mean)] == "green" ~ "alpha",
                                 channel[which.max(mean[-which.max(mean)])] == "green" & mean[which.max(mean[-which.max(mean)])] > 0.5*mean[which.max(mean)] ~ "alpha",
                                 channel[which.min(mean)] == "green" & mean[which.max(mean[-which.max(mean)])] > 0.9*mean[which.max(mean)] ~ "alpha",
                                 channel[which.max(mean)] == "red" ~ "beta",
                                 channel[which.max(mean)] == "blue" ~ "beta"
                                 )
    ) %>%
    distinct(index, .keep_all = T) %>%
    select(c(index, cell_type))

# Verify cell type proportions
data_2 %>% 
    ungroup() %>% 
    group_by(cell_type) %>% 
    tally()


# Add cell types to ROI polygons
roi_all_sf_polygon$cell_type <- data_2$cell_type

# Save segmented cell ROIs
save(roi_all_sf_polygon, file = "output/RD1-ROI_mapping_and_cell_type_assignment/ROI_polygons_with_cell_types.RData")

## Visualize cell type assignments ----
roi_all_sf_polygon_plot <-  rename(roi_all_sf_polygon, "Cell Type" = cell_type)

ggplot(roi_all_sf_polygon_plot)+
    geom_sf(aes(fill = `Cell Type`))+
    scale_fill_brewer(palette = "Dark2")+
    theme(legend.title = element_text("Cell Type"))+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("output/RD1-ROI_mapping_and_cell_type_assignment/cell_type_assignment.png",
dpi = 300, width = 1570, height = 1100, units = "px")
