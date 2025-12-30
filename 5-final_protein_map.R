library(tidyverse)
library(sf)
library(glue)



# Import Data -------------------------------------------------------------


load("output/RD3-ROI_and_pixel_to_MS/RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData")
pix_to_ms <- old


(load("RD4_1-Proteins_w_no_RNASeq.RData"))

(load("RD4_2-RNA-Seq_cell_protein_ratios.RData"))


gene_name <- read_tsv(file = "../../../../Uniprot Downloads/Human_sp_tr_primary_gene_name_20220531.tab")

# Process -----------------------------------------------------------------

# Distribute protein intensity using RNA-Seq data
pix_to_ms_data_save <- filter(pix_to_ms, !is.na(cell.roi))

pix_to_ms_data <- filter(pix_to_ms, !is.na(cell.roi))

pix_to_ms_na <- filter(pix_to_ms, is.na(cell.roi))

# Observed proteins without RNASeq data
c1 <- colnames(pix_to_ms)
c2 <- colnames(protein_ratios)
missing_rna <- anti_join(data.frame(uniprot=c1), data.frame(uniprot=c2)) %>% 
    slice(8:nrow(.)) %>% 
    left_join(gene_name, by = c("uniprot" = "Entry")) %>% 
    select(-Status)
# writexl::write_xlsx(x = missing_rna, path = "Missing_RNASeq_Proteins.xlsx")

# This takes >2hrs
if(!file.exists("RD5_1-RNA-Seq_distributed_protein_absolute_intensities.RData")){

    for(j in colnames(select(protein_ratios, 2:last_col()))){
        for(i in 1:nrow(pix_to_ms_data)){
            pix_to_ms_data[i, j] <- pix_to_ms_data[[i, j]] * protein_ratios[protein_ratios$celltype == pix_to_ms_data[[i, "cell_type"]], j]
        }
    }
    
    pix_to_ms <- bind_rows(pix_to_ms_data, pix_to_ms_na)
    
    save(pix_to_ms, file = "RD5_1-RNA-Seq_distributed_protein_absolute_intensities.RData")
} else {
    load("RD5_1-RNA-Seq_distributed_protein_absolute_intensities.RData")
}

plot(pix_to_ms)


low <- "midnightblue"
high <- "yellow"

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

# cell type
ggplot(source)+
    geom_sf(aes(fill = cell_type))+
    # scale_fill_gradient(
    #     low = low,
    #     high = high)+
    ggtitle("Cell type")+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

# Pixel
ggplot(pix_to_ms_data)+
    geom_sf(aes(fill = pixel))+
    # scale_fill_gradient(
    #     low = low,
    #     high = high)+
    ggtitle("Pixel")+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")+
    labs(fill = "Relative Intensity")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

# Secretogranin-2 (P13521)
ggplot(old)+
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

pix_to_ms <- bind_rows(not_na_relative, na)

# Check missingness
as.data.frame(pix_to_ms) %>%
    select(8:ncol(.)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x))})) %>% 
    as.numeric(.[1,]) %>% 
    hist(100)

# get proteins with 100% missingness created by relative abundance calc
# these are created by 0/0 = NaN after x / max(x)
new_missing <- as.data.frame(pix_to_ms) %>% 
    select(8:ncol(.)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x)) == 1})) %>% 
    select(which(as.logical(as.vector(.[1,])) == TRUE))

# Remove proteins w/ no data
pix_to_ms <- pix_to_ms %>% 
    select(!colnames(new_missing))


if(!file.exists("RD5_2-RNA-Seq_distributed_protein_relative_intensities.RData")){
    save(pix_to_ms, file = "RD5_2-RNA-Seq_distributed_protein_relative_intensities.RData")
} else {
    load("RD5_2-RNA-Seq_distributed_protein_relative_intensities.RData")
}


low <- "midnightblue"
high <- "yellow"

# low <- "#0072B2"
# high <- "#F0E442"
# high <- "#D55E00"

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

# ggsave("./Final_Sample_Images/Img_Pro-glucagon_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")

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

# ggsave("./Final_Sample_Images/Img_Insulin_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")

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

# ggsave("./Final_Sample_Images/Img_Secretogranin-2_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


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

# ggsave("./Final_Sample_Images/Img_Lithostathine-1-alpha_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


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

# ggsave("./Final_Sample_Images/Img_Lithostathine-1-beta_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


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

# ggsave("./Final_Sample_Images/Img_Carboxypeptidase_A1_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")


# Save all ~2000 protein maps ----

# get gene names
uniprotID <- pix_to_ms %>% 
    as.data.frame() %>% 
    filter(slide == "x60") %>% 
    select(-c(slide:alpha_cell_count))

uniprotID <- uniprotID[, !is.na(colMeans(uniprotID, na.rm=T))]

uniprotID <- uniprotID %>%
    colnames() %>% 
    as.data.frame() %>% 
    rename(uniprot = colnames(.)[[1]])

gene_name <- gene_name %>% 
    filter(grepl("HUMAN", `Entry name`))

# join and spot fix missing gene name (missing uniprot accession too)
protein_names <- left_join(uniprotID, gene_name, by = c("uniprot" = "Entry")) %>% 
    mutate(`Gene names  (primary )` = if_else(is.na(`Gene names  (primary )`),
                                            uniprot,
                                            `Gene names  (primary )`)) %>% 
    mutate(`Gene names  (primary )` = sub("(^.*?);.*", "\\1", `Gene names  (primary )`))


# list of proteins w/ no RNASeq data (to skip)
skip <- sct_counts2_remain


calc_missing <- \(x) {
    select_prot <- pix_to_ms[, c(x, "pixel")] %>% 
    filter(!is.na(pixel))
    
    missingness <- mean(is.na(select_prot[[x]]))
    return(missingness)
}

calc_missing_v <- Vectorize(calc_missing)

# Remove proteins with no RNASeq data or too much missingness (skip plotting)
protein_names <- protein_names %>% 
    select(uniprot, `Gene names  (primary )`) %>% 
    filter(!protein_names$uniprot %in% skip$Entry) %>% 
    mutate(missingness = calc_missing_v(uniprot)) 

# save(protein_names, file = "RD5_4-protein_names_with missingness.RData")


plot_maps <- \(x, y) {
    
    
    pix_to_ms_floop <- pix_to_ms[, x]
    
        p <- ggplot(pix_to_ms_floop)+
                geom_sf(aes(fill = eval(parse(text = x))))+
                scale_fill_gradient(
                    low = low,
                    high = high)+
                ggtitle(glue("Protein: {y}, Uniprot ID: {x}"))+
                theme(plot.title = element_text(hjust = 0.5))+
                labs(fill = "Relative Intensity")+
                theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      plot.title = element_text(size = 20))
        


        return(p)

}


save_maps <- FALSE

if (save_maps) {
    
    # write maps to RData (takes >30 min)
    if (!file.exists("RD5_3-All_protein_maps.RData")) {
        maps <- map2(.x = protein_names$uniprot, .y = protein_names$`Gene names  (primary )`, 
                                        .f = plot_maps,
                                        .progress = TRUE)
        
        names(maps) <- protein_names$`Gene names  (primary )`
    
        save(maps, file = "RD5_3-All_protein_maps.RData")
    } else {
        load("RD5_3-All_protein_maps.RData")
    }

    
    # write maps to png (takes >30 min) 
    missing_cut <- 0.5
    
    plots_to_save <- protein_names %>% 
        filter(missingness < missing_cut)
    
    save_plots <- \(x, y) {
        ggsave(filename = glue("{y}_{x}.png"),
               plot = maps[[y]],
               path = "Slide_61_Final_Protein_Maps_missingness_50pct_Unimputed",
               device = "png",
               dpi = 300,
               width = 4,
               height = 3,
               units = "in")
    }
    
    walk2(.x = plots_to_save$uniprot, 
          .y = plots_to_save$`Gene names  (primary )`,
          .f = save_plots)
}

options(error=recover)

# Plots for HuBMAP Update Presentation ------------------------------------


library(gridExtra)
load("RD5_2-RNA-Seq_distributed_protein_relative_intensities.RData")


low <- "midnightblue"
high <- "yellow"

source <-  pix_to_ms

## Islet Markers
# Pro-glucagon (P01275)
GCG <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")

# ggsave("./Final_Sample_Images/Img_Pro-glucagon_RNA-Seq.tiff", device = "tiff", dpi = 300, width = 1570, height = 1180, units = "px")

# Insulin (P01308)
INS <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")


# Secretogranin-2 (P13521)
SCG2 <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")


grob1 <- grid.arrange(SCG2, INS, GCG, nrow = 1)

ggsave(filename = "../HuBMAP_Updates/Islet Protein Markers.tiff", 
       plot = grob1, 
       device = "tiff",
       width = 2315,
       height = 947,
       units = "px")



## Acinar Markers
# Lithostathine-1-alpha (P05451)
REG1A <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")



# Lithostathine-1-beta (P48304)
REG1B <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")



# Carboxypeptidase A1 (P15085)
CPA1 <- ggplot(source)+
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
          axis.ticks.y=element_blank(),
          legend.position = "none")


grob2 <- grid.arrange(REG1A, REG1B, CPA1, nrow = 1)

# ggsave(filename = "../HuBMAP_Updates/Acinar Protein Markers.tiff", 
#        plot = grob2, 
#        device = "tiff",
#        width = 2315,
#        height = 947,
#        units = "px")







