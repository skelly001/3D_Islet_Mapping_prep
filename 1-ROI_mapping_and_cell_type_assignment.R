# library(plyr)
library(tidyverse)
library(RImageJROI)
# library(spatstat)
# library(sp)
library(sf)
# library(maptools)
library(RColorBrewer)

# Import ------------------------------------------------------------------


# Read ROI RGB color measurements *WITH EDGE ROIs*
file.list <- list.files(path = r"(data\1-ImageJ_cell_ROI\ROI_RGB_values)",
                        full.names = T)
# file.list <- list.files(path = r"(Slide61_cell_seg\Fiji_cell_seg\Slide61_fiji_seg_w_edges\Slide61-islet23_wcellbnd_edge_crop\Imagej_py_Slide61_cropped_ROI_RGB_values)",
#                         full.names = T)

data <- map(file.list, data.table::fread, data.table = F)


# Read ROI from mask *WITH EDGE ROIs*
## Note: using develop branch of skelly001/RImageJROI (to allow for dtrType == "traced")
roi_all <- read.ijzip(file = r"(data\1-ImageJ_cell_ROI\ROI.zip)", 
                      verbose = F)
# roi_all <- read.ijzip(file = r"(Slide61_cell_seg\Fiji_cell_seg\Slide61_fiji_seg_w_edges\Slide61-islet23_wcellbnd_edge_crop\ROI.zip)", 
#                       verbose = F)


# Process ----------------------------------------------------------------

dir.create("output/RD1-ROI_mapping_and_cell_type_assignment")


# Create ROI polygons *WITH EDGE ROIs* ------------------------------------------------------------
## Code copied from 0-import_data.R

# Long format x,y data
roi_all_coords <- vector("list", length(roi_all))

for (i in 1:length(roi_all)) {
    roi_all_coords[[i]] <- as.data.frame(roi_all[[i]][["coords"]]) %>%
        mutate(roi = names(roi_all[i]))
}

roi_all_coords <- bind_rows(roi_all_coords)

#Flip image 180 degrees across the x-axis (R plots the image upside down..)
roi_all_coords[, "y"] <- roi_all_coords[, "y"]*-1

# Create sf obj
roi_all_sf <- st_as_sf(roi_all_coords, coords = c("x","y"))

# Points to polygon
roi_all_sf_polygon <- roi_all_sf %>% 
    group_by(roi) %>%
    summarise() %>% 
    st_cast("POLYGON") %>% 
    st_convex_hull()

# plot(roi_all_sf_polygon)


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



# Correlate ROI with color ------------------------------------------------

# Extract file order from file name
file.list.2 <- file.list %>% 
    as.data.frame() %>% 
    mutate(index = str_extract(`.`, "[[:digit:]]{1,4}.csv")) %>% 
    mutate(index = sub(".csv", "", index)) %>% 
    select(index) %>% 
    map_df(as.integer)

# Add file order column to each file
for (i in 1:nrow(file.list.2)) {
    data[[i]] <- data[[i]] %>% 
        mutate(index = rep(file.list.2$index[[i]], 3))
}

# Single dataframe and arrange files in correct order (low to high)
data <- bind_rows(data) %>% 
    group_by(index) %>% 
    arrange(index) %>% 
    ungroup()



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
    
# Check proportion of cell types
data_2 %>% 
    ungroup() %>% 
    group_by(cell_type) %>% 
    tally()


# Combine cell_type with ROI polygons
roi_all_sf_polygon$cell_type <- data_2$cell_type

# Save segmented cell ROIs in correct order (i.e., order in which they were created)
save(roi_all_sf_polygon, file = "output/RD1-ROI_mapping_and_cell_type_assignment/ROI_polygons_with_cell_types.RData")


## plot ----
# library(RColorBrewer)
# library(extrafont)
# font_import()
# loadfonts(device = "win")
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


# scale_fill_discrete(type = c("pink", "yellow", "midnightblue"))

# plot(select(roi_all_sf_polygon, "cell_type"))



