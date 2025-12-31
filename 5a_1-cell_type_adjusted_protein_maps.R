# ==============================================================================
# Cell Type-Adjusted Spatial Protein Map Generation
# ==============================================================================
# Script: 5a_1-cell_type_adjusted_protein_maps.R
# Description: Generates final spatial protein abundance maps by applying
#              RNA-Seq-derived cell type expression ratios to MS protein data.
#              Calculates slide-normalized relative abundances and creates
#              QC visualization plots for islet (Glucagon, Insulin, Secretogranin-2)
#              and acinar (REG1A, REG1B, CPA1) marker proteins. Exports individual
#              protein maps for all proteins meeting quality thresholds.
#
# Input: - output/RD3-ROI_and_pixel_to_MS/RD3-Slide61_cell_roi_ms_absolute.rds
#        - output/RD4-RNA-Seq/RD4-RNA-Seq_cell_gene_ratios.rds
#        - output/RD4-RNA-Seq/RD4-proteins_no_RNA-Seq.rds
#        - data/Human_sp_tr_primary_gene_name_20220531.tab
# Output: - output/RD5-final_protein_maps/islet_protein_markers.png
#         - output/RD5-final_protein_maps/acinar_protein_markers.png
#         - output/RD5-final_protein_maps/final_protein_maps/*.png
# ==============================================================================

library(gridExtra)
library(tidyverse)
library(sf)
library(furrr)
plan(multisession, workers = 2)



# Import Data -------------------------------------------------------------

# ROI with Protein abundance
pix_to_ms <- readRDS("output/RD3-ROI_and_pixel_to_MS/RD3-Slide61_cell_roi_ms_absolute.rds")

# Proteins with no RNA-seq data
no_rnaseq <- readRDS("output/RD4-RNA-Seq/RD4-proteins_no_RNA-Seq.rds")

# Cell type gene expression ratios
protein_ratios <- readRDS("output/RD4-RNA-Seq/RD4-RNA-Seq_cell_gene_ratios.rds")

# Uniprot accession to gene name key
gene_name <- read_tsv(file = "data/Human_sp_tr_primary_gene_name_20220531.tab")


# Process -----------------------------------------------------------------

dir.create("output/RD5-final_protein_maps")


## Apply gene expression ratios to protein abundance -----

# Annotation columns
anno_cols <- c("pixel", "cell_type", "cell.roi", "geometry", 
"acinar_cell_count", "alpha_cell_count", "beta_cell_count") 

# Observed proteins without RNASeq data
c1 <- colnames(pix_to_ms)
c2 <- colnames(protein_ratios)
missing_rnax <- anti_join(data.frame(uniprot=c1), data.frame(uniprot=c2)) %>% 
    filter(!uniprot %in% anno_cols) %>% 
    left_join(gene_name, by = c("uniprot" = "Entry")) %>% 
    select(-Status)

# Create a lookup matrix from protein_ratios indexed by celltype
ratio_matrix <- protein_ratios %>%
	column_to_rownames("celltype") %>%
	as.matrix()

protein_cols <- colnames(ratio_matrix)

# Get the ratio for each row based on its cell_type
cell_types <- pix_to_ms$cell_type
ratio_lookup <- ratio_matrix[cell_types, protein_cols, drop = FALSE]

# Extract protein data as matrix, multiply element-wise, put back
protein_data <- as.matrix(st_drop_geometry(pix_to_ms)[, protein_cols])
pix_to_ms[, protein_cols] <- protein_data * ratio_lookup



## Calculate relative abundance -----
pix_to_ms <- mutate(pix_to_ms, slide = str_extract(pixel, "^x\\d{2}")) %>% 
    relocate(slide, .before = colnames(.)[1])

na <- pix_to_ms %>% 
    filter(is.na(slide))

not_na <- pix_to_ms %>% 
    filter(!is.na(slide))

anno_cols_ext <- c("slide", anno_cols)

not_na_relative <- not_na %>% 
    group_by(slide) %>% 
    mutate(across(!matches(anno_cols_ext), \(x) {x / max(x, na.rm = T)})) %>%
    ungroup()

pix_to_ms <- bind_rows(not_na_relative, na)

# Check missingness
as.data.frame(pix_to_ms) %>%
    select(!all_of(anno_cols_ext)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x))})) %>% 
    as.numeric(.[1,]) %>% 
    hist(100)

# get proteins with 100% missingness created by relative abundance transformation
to_remove <- as.data.frame(pix_to_ms) %>% 
    select(!all_of(anno_cols_ext)) %>% 
    summarise(across(everything(), \(x) {mean(is.na(x)) == 1})) %>% 
    select(which(as.logical(as.vector(.[1,])) == TRUE))

# Remove proteins with no data
pix_to_ms <- pix_to_ms %>% 
    select(!colnames(to_remove))



## Visualize islet/acinar protein marker abundance -----
low <- "midnightblue"
high <- "yellow"

plot_protein_map <- \(data, protein_col, title) {
    ggplot(data) +
        geom_sf(aes(fill = .data[[protein_col]])) +
        scale_fill_gradient(low = low, high = high) +
        ggtitle(title) +
        labs(fill = "Relative Intensity") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "none")
}

# Islet Markers
GCG <- plot_protein_map(pix_to_ms, "P01275", "Glucagon")
INS <- plot_protein_map(pix_to_ms, "P01308", "Insulin")
SCG2 <- plot_protein_map(pix_to_ms, "P13521", "Secretogranin-2")

grob_islet <- grid.arrange(SCG2, INS, GCG, nrow = 1)
ggsave(filename = "output/RD5-final_protein_maps/islet_protein_markers.png", 
       plot = grob_islet, width = 2315, height = 947, units = "px")

# Acinar Markers
REG1A <- plot_protein_map(pix_to_ms, "P05451", "Lithostathine-1-alpha")
REG1B <- plot_protein_map(pix_to_ms, "P48304", "Lithostathine-1-beta")
CPA1 <- plot_protein_map(pix_to_ms, "P15085", "Carboxypeptidase A1")

grob_acinar <- grid.arrange(REG1A, REG1B, CPA1, nrow = 1)
ggsave(filename = "output/RD5-final_protein_maps/acinar_protein_markers.png", 
       plot = grob_acinar, width = 2315, height = 947, units = "px")



# Save all protein maps ----

# get gene names
uniprotID <- pix_to_ms %>% 
    as.data.frame() %>% 
    filter(slide == "x60") %>% 
    select(!all_of(anno_cols_ext))

uniprotID <- uniprotID[, !is.na(colMeans(uniprotID, na.rm=T))]

uniprotID <- uniprotID %>%
    colnames() %>% 
    as.data.frame() %>% 
    rename(uniprot = colnames(.)[[1]])

gene_name <- gene_name %>% 
    filter(grepl("HUMAN", `Entry name`))

# join and spot fix missing gene name/uniprot accession
protein_names <- left_join(uniprotID, gene_name, by = c("uniprot" = "Entry")) %>% 
    mutate(`Gene names  (primary )` = if_else(is.na(`Gene names  (primary )`),
                                            uniprot,
                                            `Gene names  (primary )`)) %>% 
    mutate(`Gene names  (primary )` = sub("(^.*?);.*", "\\1", `Gene names  (primary )`))


# Compute protein missingness
calc_missing <- \(x) {
    select_prot <- pix_to_ms[, c(x, "pixel")] %>% 
    filter(!is.na(pixel))
    
    missingness <- mean(is.na(select_prot[[x]]))
    return(missingness)
}

calc_missing_v <- Vectorize(calc_missing)

# Remove proteins with no RNASeq data
protein_names <- protein_names %>% 
    select(uniprot, `Gene names  (primary )`) %>% 
    filter(!.$uniprot %in% no_rnaseq$Entry) %>% 
	filter(!duplicated(.$`Gene names  (primary )`)) %>% 
    mutate(missingness = calc_missing_v(uniprot))


## write maps to png -----

# Apply missingness filter (<50%)
missing_threshold <- 0.5
plots_to_save <- protein_names %>% 
	filter(missingness < missing_threshold)  

dir.create("output/RD5-final_protein_maps/final_protein_maps")

# Helper function for saving protein maps
save_plots <- \(x, y) {
	
	pix_to_ms_floop <- pix_to_ms[, c(x, "geometry")]

	p <- ggplot(pix_to_ms_floop)+
		geom_sf(aes(fill = .data[[x]]))+
		scale_fill_gradient(
			low = low,
			high = high)+
		labs(fill = "Relative Intensity")+
		theme(axis.text.x=element_blank(),
				axis.ticks.x=element_blank(),
				axis.text.y=element_blank(),
				axis.ticks.y=element_blank(),
				plot.title = element_text(size = 20))
	
	png(filename = str_glue("output/RD5-final_protein_maps/final_protein_maps/{y}_{x}.png"),
		width = 1570,
		height = 1100,
		res = 300)
	plot(p)  
	dev.off()

}

# Save maps
future_walk2(.x = plots_to_save$uniprot, 
		.y = plots_to_save$`Gene names  (primary )`,
		.f = save_plots, 
		.progress = TRUE,
		.options = furrr_options(seed = NULL))
