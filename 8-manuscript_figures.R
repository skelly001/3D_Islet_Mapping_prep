library(tidyverse)
library(sf)
library(glue)
library(gridExtra)



# Import ------------------------------------------------------------------

# pixel-level proteomics maps
(load("RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData"))


# pseudo-single cell proteomics maps
(load("RD5_2-RNA-Seq_distributed_protein_relative_intensities.RData"))




# Plots for manuscript fig 5 ----------------------------------------------

dir.create("manuscript_images")



## pseudo-single cell maps ----

low <- "midnightblue"
high <- "yellow"

source <-  pix_to_ms

width <- 8.5/4
height <- 8.5/4


## Islet Markers
# Pro-glucagon (P01275)
GCG <- ggplot(source)+
    geom_sf(aes(fill = P01275), , lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Pro-glucagon_RNA-Seq.tiff", GCG, device = "tiff", dpi = 600, width = width, height = height, units = "in")

# for the scale bar
# Pro-glucagon (P01275)
GCG_scale_bar <- ggplot(source)+
    geom_sf(aes(fill = P01275), , lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./manuscript_images/Img_Pro-glucagon_RNA-Seq_scale_bar.tiff", GCG_scale_bar, device = "tiff", dpi = 600, width = width, height = height, units = "in")

# Insulin (P01308)
INS <- ggplot(source)+
    geom_sf(aes(fill = P01308), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Insulin_RNA-Seq.tiff", plot = INS, device = "tiff", dpi = 600, width = width, height = height, units = "in")



# Secretogranin-2 (P13521)
SCG2 <- ggplot(source)+
    geom_sf(aes(fill = P13521), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Secretogranin_2_RNA-Seq.tiff", plot = SCG2, device = "tiff", dpi = 600, width = width, height = height, units = "in")


# Chromogranin-A (P10645)
CHGA <- ggplot(source)+
    geom_sf(aes(fill = P10645), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Chromogranin_A_RNA-Seq.tiff", plot = CHGA, device = "tiff", dpi = 600, width = width, height = height, units = "in")


# Carboxypeptidase E (P16870)
CPE <- ggplot(source)+
    geom_sf(aes(fill = P16870), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Carboxypeptidase_E_RNA-Seq.tiff", plot = CPE, device = "tiff", dpi = 600, width = width, height = height, units = "in")






## Acinar Markers
# Lithostathine-1-alpha (P05451)
REG1A <- ggplot(source)+
    geom_sf(aes(fill = P05451), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Lithostathine-1-alpha_RNA-Seq.tiff", plot = REG1A, device = "tiff", dpi = 600, width = width, height = height, units = "in")




# Carboxypeptidase A1 (P15085)
CPA1 <- ggplot(source)+
    geom_sf(aes(fill = P15085), lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Carboxypeptidase_A1_RNA-Seq.tiff", plot = CPA1, device = "tiff", dpi = 600, width = width, height = height, units = "in")


# cell type
ggplot(source)+
    geom_sf(aes(fill = cell_type))+
    # scale_fill_gradient(
    #     low = low,
    #     high = high)+
    ggtitle("Cell type")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Cell Type")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    paletteer::scale_fill_paletteer_d("ggthemes::wsj_rgby")

# cell type
cell_type_plot <- ggplot(source)+
    geom_sf(aes(fill = cell_type), lwd = 0.165)+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Cell Type")+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"))
    # scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02"))

ggsave("./manuscript_images/Img_Cell_type_RNA-Seq.tiff", plot = cell_type_plot, 
       device = "tiff", dpi = 600, width = width+1, height = height, units = "in")






## pixel level maps ----

low <- "midnightblue"
high <- "yellow"

source <-  old

width <- 8.5/4
height <- 8.5/4


# Insulin (P01308)
INS_pixel <- ggplot(source)+
    geom_sf(aes(fill = P01308), , lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")

ggsave("./manuscript_images/Img_Insulin_pixel_level_RNA-Seq.tiff", INS_pixel, device = "tiff", dpi = 600, width = width, height = height, units = "in")

# Glucagon (P01275)
GCG_pixel <- ggplot(source)+
    geom_sf(aes(fill = P01275), , lwd=0.1)+
    scale_fill_gradient(
        low = low,
        high = high)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggsave("./manuscript_images/Img_Glucagon_pixel_level_RNA-Seq.tiff", GCG_pixel, device = "tiff", dpi = 600, width = width, height = height, units = "in")









# cellpose mask -----------------------------------------------------------


library(RImageJROI)
library(sf)
library(tidyverse)

roi_all <- read.ijzip(file = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Map\Prototype_Slide_61\61_Unimputed_Maps\manuscript_images\cellpose_mask\RoiSet.zip)", 
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

plot(roi_all_sf_polygon)


pal <- paletteer::scale_fill_paletteer_d("fishualize::Halichoeres_brasiliensis")
pal <- paletteer::scale_fill_paletteer_d("khroma::smoothrainbow")
pal <- pal$palette()

ggplot(roi_all_sf_polygon)+
    geom_sf(aes(fill = roi), lwd=0, alpha=0.7)+
    theme_void()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.background = element_blank(),
          # plot.background = element_rect(fill = 'black', colour = 'black'),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.border = element_blank())+
    scale_fill_manual(values = rep(pal, nrow(roi_all_sf_polygon)))


ggsave("./manuscript_images/cellpose_mask/DAPI_mask_for_fig_5.png", width = 2.125, height = 2.125, dpi = 600)























