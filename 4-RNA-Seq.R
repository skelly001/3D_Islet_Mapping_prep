library(tidyverse)
library(sf)
library(Seurat)

# Note: counts[, 1:50] is used only for quick prototyping. Not used for final product (lines 50/91)
# Import Data -------------------------------------------------------------


# MS data
load("RD3_1-Slide61_roi_w_edges_pixels_ms_absolute.RData")
ms_abs <- old
rm(old)

# RNA-Seq reference data
fullref <- readRDS(r"(..\..\..\azimuth-references\human_pancreas_snakemake\seurat_objects\fullref.Rds)")

# Uniprot accession to gene correlation
# uniprot <- read.delim("../../../../Uniprot Downloads/Human_sp_tr_primary_gene_name_20220531.tab",
#                       sep = "\t")

uniprot <- data.table::fread(file = "../../../../Uniprot Downloads/uniprotkb_Human_2023_10_25.tsv.gz", 
                             data.table = F,
                             stringsAsFactors = F,
                             showProgress = T,
                             na.strings = c("", "NA"), 
                             check.names = T)

# Process -----------------------------------------------------------------

# RNA-Seq meta data
meta <- fullref@meta.data %>% 
    mutate(cell = rownames(.)) %>% 
    select(c(cell, celltype)) %>% 
    `rownames<-`(NULL)

# filter RNA-Seq (counts) ref data for entries in MS data
ms_proteins <- as.data.frame(ms_abs) %>% 
    select(-c(pixel:alpha_cell_count)) %>% 
    colnames()

## Prep uniprot gene names (some RNASeq entries using alternate gene names)
uniprot_sp_pirimary <- uniprot %>% 
    filter(Reviewed == "reviewed") %>% 
    select(Entry, Gene.Names..primary.) %>% 
    filter(Entry %in% ms_proteins) %>% 
    separate_longer_delim(cols = Gene.Names..primary.,
                          delim = "; ")

## first match, using primary gene name
sct_counts1 <- as.matrix(fullref@assays[["SCT"]]@counts) %>% #  fullref@assays[["SCT"]]@counts[, 1:50]
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    `rownames<-`(NULL) %>% 
    right_join(uniprot_sp_pirimary, by = c("gene" = "Gene.Names..primary.")) %>% 
    relocate(c(Entry, gene), .before = colnames(.)[1]) %>% 
    mutate(missing = if_else(is.na(rowSums(.[-c(1:2)])), TRUE, FALSE), .after = gene) %>% 
    arrange(missing) %>% 
    distinct(Entry, .keep_all = T)

### successful matches
sct_counts1_success <- sct_counts1 %>% 
    filter(missing == FALSE) %>% 
    select(-missing)

### no match after using primary gene name
sct_counts1_remain <-  sct_counts1 %>% 
    filter(missing == TRUE) %>% 
    select(Entry)

## prep uniprot gene names for alt. name matching
uniprot_sp_alt <- uniprot %>% 
    filter(Reviewed == "reviewed") %>% 
    select(-c(Reviewed, Entry.Name, Protein.names, Organism, Length)) %>% 
    filter(Entry %in% sct_counts1_remain$Entry) %>% 
    unite(c(Gene.Names, Gene.Names..ordered.locus.,
            Gene.Names..ORF., Gene.Names..synonym.),
          col = alt_gene_name,
          sep = "; ",
          na.rm = T) %>%
    separate_longer_delim(cols = alt_gene_name,
                          delim = "; ") %>%
    separate_longer_delim(cols = alt_gene_name,
                          delim = " ") %>% 
    # filter(alt_gene_name != Gene.Names..primary.) %>% #these names have already failed to match but will leave them so they end up in sct_counts2_remain, for convinience
    select(-Gene.Names..primary.) %>% 
    distinct(Entry, alt_gene_name) %>% 
    mutate(alt_gene_name = sub(";", "", alt_gene_name))


## second match, using only unmatched accession & alt. gene names
sct_counts2 <- as.matrix(fullref@assays[["SCT"]]@counts) %>% #  fullref@assays[["SCT"]]@counts[, 1:50]
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    `rownames<-`(NULL) %>% 
    right_join(uniprot_sp_alt, by = c("gene" = "alt_gene_name")) %>% 
    relocate(c(Entry, gene), .before = colnames(.)[1]) %>% 
    mutate(missing = if_else(is.na(rowSums(.[-c(1:2)])), TRUE, FALSE), .after = gene) %>% 
    arrange(missing) %>% 
    distinct(Entry, .keep_all = T)

### successful matches
sct_counts2_success <- sct_counts2 %>% 
    filter(missing == FALSE) %>% 
    select(-missing)

### no match after using primary & alt gene names ()
sct_counts2_remain <-  sct_counts2 %>% 
    filter(missing == TRUE) %>% 
    select(Entry, gene)

if (!file.exists("RD4_1-Proteins_w_no_RNASeq.RData")) {
    save(sct_counts2_remain, file = "RD4_1-Proteins_w_no_RNASeq.RData")
}


## combine matching results
sct_counts <- bind_rows(sct_counts1_success, sct_counts2_success, sct_counts2_remain) %>% 
    select(-gene)


    
# relate RNA-Seq counts to meta data
sct_counts <- sct_counts %>% 
    t() %>%
    as.data.frame() %>% 
    `colnames<-`(.[1,]) %>% 
    mutate(cell = rownames(.), .before = colnames(.)[1]) %>% #mult uniprot genes are correlated to single proteins. After checking RNA_seq, I will keep only the first observation in Gene col. from uniprot_clean
    slice(2:nrow(.)) %>% 
    `rownames<-`(NULL)

sct_counts_meta <- sct_counts %>% 
    left_join(meta, by = "cell") %>% 
    relocate(celltype, .before = cell) %>% 
    select(-cell) %>% 
    mutate(across(!matches("celltype"), as.numeric))

protein_distributions <- sct_counts_meta %>% 
    group_by(celltype) %>% 
    summarise(across(everything(), sum)) %>% 
    ungroup() %>% 
    filter(celltype %in% c("acinar", "alpha", "beta"))

# sct_counts_meta %>% 
#     group_by(celltype) %>%
#     summarise(across(everything(), sum)) %>%
#     ungroup() %>%
#     write_rds(file = "RD4_3-RNA-Seq_counts_all_cell_types.rds")

RNA_Seq_NA_proteins <-  protein_distributions[, colSums(is.na(protein_distributions)) != 0] %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    select(-c("V1", "V2", "V3")) %>% 
    rename(missing_proteins = rowname)

protein_distributions <-  protein_distributions[, colSums(is.na(protein_distributions)) == 0]
# NOTE: after removing missing RNA-Seq data, we are left with 1514 proteins (originally 1560 proteins)

protein_ratios <- protein_distributions %>% 
    mutate(across(where(~ !is.character(.x) && sum(.x) != 0), \(x){x / sum(x)}))


# CHECK islet markers: Pro-glucagon (P01275); Insulin (P01308); Secretogranin-2 (P13521);
prot_islet <- protein_ratios %>% 
    select(c("celltype", "P01275", "P01308", "P13521"))

# CHECK acinar markers: Lithostathine-1-alpha (P05451); Lithostathine-1-beta (P48304);
# Carboxypeptidase A1 (P15085)
prot_acinar <- protein_ratios %>% 
    select(c("celltype", "P05451", "P48304", "P15085"))


# Save
if(!file.exists("RD4_2-RNA-Seq_cell_protein_ratios.RData")){
    save(protein_ratios, file = "RD4_2-RNA-Seq_cell_protein_ratios.RData")
}




# Covariation between cell number per pixel & protein intensity
# NOTE: I do not think this will work because we do not have reliable cell count data
# ms_abs <- ms_abs[!is.na(ms_abs$pixel), ]

# cell_summ <- as.data.frame(ms_abs) %>% 
#     group_by(pixel, cell_type) %>% 
#     tally() %>% 
#     rename(cell_count = n)
#  
# cor <- as.data.frame(ms_abs) %>% 
#     select(-c(cell.roi, geometry)) %>% 
#     select(-contains("cell_count")) %>% 
#     filter(!is.na(pixel))

# cor <- full_join(cor, cell_summ, by = c("pixel", "cell_type")) %>% 
#     relocate(cell_count, .after = cell_type)









# missing proteins --------------------------------------------------------

main <- readxl::read_xlsx(path = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Data\MS_output\processed\3D_mapping_no_imputation.xlsx)",
                          sheet = "noinputation_combined_r")

subset <- readxl::read_xlsx(path = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Data\MS_output\processed\3D_mapping.xlsx)",
                          sheet = "log_inpu_r")



main <- main %>% 
    select(PROTID) %>% 
    mutate(PROTID = sub("sp\\|(.*?)\\|.*", "\\1", PROTID))

subset <- subset %>% 
    select(PROTID)

# proteins in subset w/o match in main (23) 
missing <- anti_join(subset, main)














