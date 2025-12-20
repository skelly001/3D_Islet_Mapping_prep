library(clusterProfiler)
library(org.Hs.eg.db)
library(UniProt.ws)
library(tidyverse)
library(hpar)
library(sf)


# Import ------------------------------------------------------------------

# scRNA-Seq counts filtered for observed proteins
rna_counts <- readRDS("RD4_3-RNA-Seq_counts_all_cell_types.rds")

# protein names with missingness
(load("RD5_4-protein_names_with missingness.RData"))

# scRNA-Seq Cell Type Signatures ------------------------------------------

# observed proteins in slide 61 w/ < 50% missingness
protein_names <- protein_names %>% 
    filter(missingness < 0.50)

# check scRNA-seq counts for proteins unique to cell types
signatures_all <- rna_counts %>% 
    filter(!is.na(celltype), !celltype %in% c("activated_stellate", 
                                              "quiescent_stellate", 
                                              "schwann", 
                                              "mast",
                                              "epsilon")) %>% 
    column_to_rownames("celltype") %>% 
    select(all_of(protein_names$uniprot))

# get cumulative counts (per protein)
signatures_all2 <- sweep(x = signatures_all, 
                         MARGIN = 2, 
                         STATS = colSums(signatures_all), 
                         FUN = "/")

# calc max-to-cumulative count ratio and filter
signatures_all3 <- signatures_all2 %>% 
    mutate(across(everything(), \(x){if_else(x >= 0.75, TRUE, FALSE)})) %>%
    colSums()

# get prot with max-to-cumulative > 0.75
signatures_all4 <- names(signatures_all3[signatures_all3 == 1 & !is.na(signatures_all3)])

# get signatures
signatures_all5 <- signatures_all %>% 
    dplyr::select(all_of(signatures_all4)) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(celltype_marker = colnames(.)[max.col(.)]) %>% 
    rownames_to_column("UniProtAcc") %>% 
    rowwise() %>%
    mutate(max_to_cumulative = max(c_across(acinar:macrophage))/sum(c_across(acinar:macrophage))) %>% 
    ungroup() %>% 
    arrange(celltype_marker, desc(max_to_cumulative))

# annotate uniprot IDs   ##  NOTE: will use uniprot DB instead b/c it is more complete
# key <- bitr(geneID = signatures_all5$UniProtAcc, 
#             fromType = "UNIPROT", 
#             toType = c("SYMBOL", "GENENAME"),
#             OrgDb = "org.Hs.eg.db")

# key1 <- AnnotationDbi::select(org.Hs.eg.db, 
#                       keys = signatures_all5$UniProtAcc,
#                       columns = c("SYMBOL", "GENENAME"), 
#                       keytype = "UNIPROT")
# 
# signatures_all6 <- left_join(signatures_all5, key1, by = c("UniProtAcc" = "UNIPROT"))


# get uniprot IDs 
up <- UniProt.ws(taxId=9606)
columns(up)
keytypes(up)

key2 <- AnnotationDbi::select(up, 
                             keys = signatures_all5$UniProtAcc,
                             columns = c("gene_primary", "id", "protein_name"),
                             keytype = "UniProtKB") %>% 
    select(-From)

signatures_all6 <- left_join(signatures_all5, key2, by = c("UniProtAcc" = "Entry"))

# writexl::write_xlsx(signatures_all6, path = "RD7_2-cell_type_signature_candidates.xlsx")



# Plot pixel level maps ---------------------------------------------------

(load("RD3_2-Slide61_roi_w_edges_pixels_ms_relative.RData"))

pix_to_ms <- filter(pix_to_ms, !is.na(cell.roi))


plot_maps <- \(x, y) {
    
    pix_to_ms_floop <- pix_to_ms[, x]
    
    p <- ggplot(pix_to_ms_floop)+
        geom_sf(aes(fill = eval(parse(text = x))))+
        scale_fill_gradient(
            low = low,
            high = high)+
        ggtitle(glue::glue("Protein: {y}, Uniprot ID: {x}"))+
        theme(plot.title = element_text(hjust = 0.5))+
        labs(fill = "Relative Intensity")+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(size = 20))
    
    return(p)
    
}
 

save_maps <- F

if (save_maps) {
    
    # write maps to RData (takes >60 min)
    if (!file.exists("RD7_1-pixel_level_protein_maps_missingness_50pct.RData")) {
        maps <- map2(.x = protein_names$uniprot, .y = protein_names$`Gene names  (primary )`, 
                     .f = plot_maps,
                     .progress = TRUE)
        
        names(maps) <- protein_names$`Gene names  (primary )`
        
        save(maps, file = "RD7_1-pixel_level_protein_maps_missingness_50pct.RData")
    } else {
        load("RD7_1-pixel_level_protein_maps_missingness_50pct.RData")
    }
    
    
    # write maps to jpeg (takes >60 min) 
    missing_cut <- 0.5
    
    plots_to_save <- protein_names %>% 
        filter(missingness < missing_cut)
    
    save_plots <- \(x, y) {
        ggsave(filename = glue::glue("Protein_{y}, Uniprot ID_{x}, No Impute.jpeg"),
               plot = maps[[y]],
               path = "Slide_61_pixel_level_protein_maps_missingness_50pct",
               device = "jpeg",
               dpi = 100,
               width = 8,
               height = 7,
               units = "in")
    }
    
    walk2(.x = plots_to_save$uniprot, 
          .y = plots_to_save$`Gene names  (primary )`,
          .f = save_plots)
}



# Copy cell type signature maps to subfolders -----------------------------

# copy cell type signature proteins to new folders
to_match <- paste0(signatures_all6$UniProtAcc, collapse = "|")

file_list <- list.files("Slide_61_pixel_level_protein_maps_missingness_50pct", full.names = T) %>% 
    as_tibble() %>%
    filter(grepl(pattern = to_match, x = value))

file.copy(from = file_list$value, to = "Slide_61_pixel_level_protein_maps_missingness_50pct/cell_type_signature_candidates/")

dir <- paste0("Slide_61_pixel_level_protein_maps_missingness_50pct/cell_type_signature_candidates/", 
              unique(signatures_all6$celltype_marker))

fs::dir_create(path = dir, recurse = T)

for (i in unique(signatures_all6$celltype_marker)) {
    
    uni <- signatures_all6[signatures_all6$celltype_marker == i, "UniProtAcc"][[1]]
    
    to_match <- paste0(uni, collapse = "|")
    
    file_list_floop <- file_list %>% 
        filter(grepl(pattern = to_match, x = value))
    
    dir_floop <- paste0("Slide_61_pixel_level_protein_maps_missingness_50pct/cell_type_signature_candidates/", i)
    
    file.copy(from = file_list_floop$value, to = dir_floop)
    
}





# Overlap with MatrisomeDB ------------------------------------------------

# Import DB
matrisome_db <- readxl::read_xlsx("Hs_Matrisome_DB_Masterlist_Naba et al_2012.xlsx") %>% 
    separate_longer_delim(cols = UniProt_IDs, delim = ":")

# cell type signatures with entries in MatrisomeDB
hits1 <- inner_join(signatures_all6, matrisome_db, by = c("UniProtAcc" = "UniProt_IDs"))




# Overlap with MatrixDB ECM -----------------------------------------------

# Import DB (had to change file ext from csv to tsv)
matrix_db <- read_tsv("MatrixDB_extracellular_matrix_proteins.tsv")

hits2 <- inner_join(signatures_all6, matrix_db, by = c("UniProtAcc" = "Uniprot primary AC"))


# all hits2 are found in hits1
dual_hits <- inner_join(hits1, hits2, by = "UniProtAcc")

# writexl::write_xlsx(hits2, "RD7_3-cell_type_signature_matrix_proteins.xlsx")



#









# cross-ref with Human Protein Atlas
nt <- hpaNormalTissue()

nt2 <- nt %>% 
    filter(Reliability == "Approved",
           Level %in% c("High", "Medium"),
           Tissue == "pancreas") 
    
counts_all3 <- left_join(counts_all2, nt2, by = c("ENSEMBL" = "Gene"))

beta <- counts_all3 %>% 
    arrange(desc(beta)) %>% 
    head()

alpha <- counts_all3 %>% 
    arrange(desc(alpha)) %>% 
    head()

acinar <- counts_all3 %>% 
    arrange(desc(acinar)) %>% 
    head()









