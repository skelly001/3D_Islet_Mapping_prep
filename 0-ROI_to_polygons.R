library(tidyverse)
library(RImageJROI)
# library(spatstat)
# library(sp)
library(sf)
library(maptools)

# Import ROI --------------------------------------------------------------


# roi1 <- read.ijroi(file = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Data\Image\Immunostaining\Fluorescence\Custom_composites\Fiji_mask\Slide61_mask_cell_boundaries\0001-0040-1.roi)", 
#                    verbose = F)  # old file path
# roi1[["strType"]] <- "polyline" #this tricks RImageJROI::ij2spatstat() to produce psp object

# Note: using develop branch of skelly001/RImageJROI (to allow for dtrType == "traced")
# *NO EDGES*
roi_all <- read.ijzip(file = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Map\Prototype_Slide_61\61_Unimputed_Maps\Slide61_cell_seg\Fiji_cell_seg\Slide61_fiji_segmentation_no_edges\segmented_cells_final\Slide61-islet23_with_cell_boundaries_ROI.zip)", 
                      verbose = F)


# Strategy 1 ------------------------------------------------------------
## This strategy works ##

# Long format x,y data
roi_all_coords <- vector("list", length(roi_all))

for (i in 1:length(roi_all)) {
    roi_all_coords[[i]] <- as.data.frame(roi_all[[i]][["coords"]]) %>%
        mutate(roi = names(roi_all[i]))
}

roi_all_coords <- bind_rows(roi_all_coords)

# Create sf obj
roi_all_sf <- st_as_sf(roi_all_coords, coords = c("x","y"))
    
# Points to polygon
roi_all_sf_polygon <- roi_all_sf %>% 
    group_by(roi) %>% 
    summarise() %>% 
    st_cast("POLYGON") %>% 
    st_convex_hull()

plot(roi_all_sf_polygon)

if (!file.exists("roi_polygons.RData")) {
    save(roi_all_sf_polygon, file = "roi_polygons.RData")
    }



# Strategy 2 --------------------------------------------------------------


# Extract x,y coords slot
roi_all_2 <- vector("list", length(roi_all))

for (i in 1:length(roi_all)) {
    names(roi_all_2)[i] <- names(roi_all[i])
    roi_all_2[[i]] <- as.data.frame(roi_all[[i]][["coords"]])
}

# create sp objects
roi_all_2 <- map(roi_all_2, SpatialPoints)

# create sf objects (point geom)
roi_all_2_sf <- map(roi_all_2, st_as_sf)

# convert to "polygon" class?
roi_all_2_sf <- st_cast(roi_all_2_sf[[1]], to = 'LINESTRING', warn = TRUE, do_split = TRUE)
roi_all_2_sf2 <- st_cast(roi_all_2_sf[[1]], to = 'POLYGON', warn = TRUE, do_split = TRUE)




# Strategy 3 --------------------------------------------------------------


#ijroi to psp
roi_all_psp <- ij2spatstat(roi_all,
            window = NULL,
            unitname = 1,
            scale = 1,
            return.type = F,
            convert.only = NULL)


#psp to spdf

# roi1.1 <- roi1 %>% as.psp.SpatialLinesDataFrame(from = "psp")

roi1.2 <- st_as_sf(roi1.1, 
                   sf_column_name = "roi1",
                   coords = )


plot(roi_all)


