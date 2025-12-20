library(readxl)
library(MSnSet.utils)
library(tidyverse)
library(cluster)
library(factoextra)
library(proteinminion)

# Import data -------------------------------------------------------------


m <- readMSnSet2(file = "../../../Protein_Data/MS_output/processed/3D_mapping_no_imputation.csv",
            ecol = 3:179,
            fnames = 1)

pixel_pheno <- readxl::read_xlsx(path = r"(C:\Users\kell343\OneDrive - PNNL\Documents\11 HuBMAP\Protein_Data\MS_output\processed\3D_clustering_result.xlsx)",
                                 sheet = "pixel_phenotype_r")


# Prep data ---------------------------------------------------------------

m1 <- m

pixel_pheno <- pixel_pheno[match(sampleNames(m1), pixel_pheno$Pixel),]

pData(m1)$pixel <- pixel_pheno$Pixel
pData(m1)$phenotypeID <- pixel_pheno$Leiden
pData(m1)$phenotype <- case_when(pData(m1)$phenotypeID == 1 ~ "Endocrine",
                                 pData(m1)$phenotypeID == 2 ~ "Peri-endocrine",
                                 pData(m1)$phenotypeID == 3 ~ "Exocrine")
pData(m1)$phenotype <- as.factor(pData(m1)$phenotype)

featureNames(m1) <- sub("_HUMAN", "", featureNames(m1))

# check sample (pixel) missingness
barplot(is.na(exprs(m1))/nrow(exprs(m1)))
plotNA(m1)

# Convert to z-score
exprs(m1) <- sweep(exprs(m1),
            MARGIN = 1,
            STATS = apply(exprs(m1), 1, mean, na.rm = TRUE),
            FUN = "-")

exprs(m1) <- sweep(exprs(m1),
            MARGIN = 1,
            STATS = apply(exprs(m1), 1, sd, na.rm = TRUE),
            FUN = "/")

# missingness filter
m1 <- m1[rowMeans(!is.na(exprs(m1))) >= 0.95,]


euc.dist <- stats::dist(exprs(m1), method = "euclidean")

# nclust <- fviz_nbclust(exprs(m1),
#                        FUNcluster = cluster::pam,
#                        method = "silhouette", #2 clust
#                        # method = "gap_stat",
#                        diss = euc.dist,
#                        k.max = 10,
#                        verbose = T)


# clusgap_res <- clusGap2(exprs(m1), FUN = cluster::pam, nstart = 5, K.max = 64, B = 64)
# Note: clusGap results - monotonically increasing

# plot(clusgap_res)


# K-Medoids Clustering ------------------------------------------------------


# res_clara <- clara(x = spearman.dist, k = 100, diss = TRUE)


# PAM with Euclidian distance
if (!file.exists("RD6_1-pam_k_50.RData")) {
    set.seed(0)
    res_pam2 <- pam(x = euc.dist, k = 50, diss = TRUE)
    save(res_pam2, file = "RD6_1-pam_k_50.RData")
} else {
    load("RD6_1-pam_k_50.RData")
}


clust <-res_pam2$clustering %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(cluster = ".",
           feature = rowname)

fData(m1) <- fData(m1) %>% 
    rownames_to_column("rowname") %>% 
    left_join(clust, by = c("rowname" = "feature")) %>% 
    column_to_rownames("rowname") %>% 
    select(-PROTID)

# save cluster results
# fData(m1) %>%
#     rownames_to_column("Feature") %>%
#     arrange(cluster) %>%
#     group_by(cluster) %>%
#     mutate(Number_of_Members = n()) %>%
#     writexl::write_xlsx(path = "RD6_2-PAM_Results_k50.xlsx")

# visualize cluster size
ggplot(fData(m1), aes(x= cluster))+ 
    geom_bar()



## PAM Cluster Heatmaps ----

# Set global ComplexHeatmap variable
library(ComplexHeatmap) #must load to adjust column anno. gap
ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "cm")
ht_opt$message = FALSE

save_hmaps <- F


hmaps <- \(i) {
    m_floop <- m1
    
    m_floop <- m_floop[fData(m_floop)$cluster == i,]
    
    height <- 1.5 + 0.18*nrow(fData(m_floop))
    
    if (save_hmaps) {
        pdf(glue::glue("./Slide61_PAM_Clustering_Heatmaps/Cluster_{i}.pdf"),
            height = height,
            width = 15,
            onefile = F)
    }
    
    complex_heatmap(eset = m_floop,
                    anno_column = "phenotype",
                    anno_column_colors = list(phenotype = c("blue", "gold", "red")),
                    clustering_method = "complete",
                    cluster_columns = T,
                    cluster_rows = T,
                    show_row_names = T,
                    heatmap_title = glue::glue("Cluster: {i}"))
    
    if (save_hmaps) {
        dev.off()
    }
}


hmaps_v <- Vectorize(hmaps)

#save heatmaps
hmaps_v(1:50)




## Pixel phenotype ANOVA (exocrine, peri-endocrine, endocrine) ----


m2 <- filter_by_occurrence(m1, 1)

fData(m2)$cluster_name <- paste0("Cluster_", fData(m2)$cluster)

m_clust <- MSnbase::combineFeatures(m2, 
                                    groupBy = fData(m2)$cluster_name,
                                    method = "mean",
                                    na.rm = T) 

library(magrittr)
limma_res <- limma_gen(m_clust, model.str = "~ phenotype", coef.str = "phenotype") %>% 
    arrange(adj.P.Val) %>% 
    filter(adj.P.Val < 0.05)


fit_res <- lapply(split(exprs(m_clust), f = featureNames(m_clust)),
                  \(y){aov(y~phenotype, data = pData(m_clust))})

# tukey_res <- lapply(fit_res, \(x){TukeyHSD(x)[["phenotype"]] %>% 
#                                   rownames_to_column("contrast") %>% 
#                                   mutate()})
# did not finish here

# tukey_res[[1]]$phenotype




# Investigate Cluster 4 ---------------------------------------------------

# Get cluster 4
clus_4 <- m1[fData(m1)$cluster == 4, ]


## Pathway enrichment ----
fData(clus_4)$uniprot <- sub(".*?\\|", "", rownames(fData(clus_4)))
fData(clus_4)$uniprot <- sub("\\|.*", "", fData(clus_4)$uniprot)

fData(m1)$uniprot <- sub(".*?\\|", "", rownames(fData(m1)))
fData(m1)$uniprot <- sub("\\|.*", "", fData(m1)$uniprot)


library(proteinminion)


if(!file.exists("Ontology_Tables/Uniprot_Homo_sapiens_proteome_UP000005640_2024_01_23.fasta")){   
    download_UniProt_Table(proteomeID = "UP000005640",
                           reviewed = T,
                           export = TRUE, 
                           file="Ontology_Tables/Uniprot_Homo_sapiens_proteome_UP000005640_2024_01_23.fasta")  
    UniProtFasta_info_Human<-Uniprot_Fasta_Parse(file = "Ontology_Tables/Uniprot_Homo_sapiens_proteome_UP000005640_2024_01_23.fasta")   
    write.csv(UniProtFasta_info_Human, "Ontology_Tables/UniProtFasta_info_human_2024_01_23.csv") } 

if(file.exists("Ontology_Tables/UniProtTables_2024_01_23.rda")){   
    load(file = "Ontology_Tables/UniProtTables_2024_01_23.rda")    
}else{   generate_Ontology_Tables(proteomeID = "UP000005640",reviewed = T)  
        save(UniProtTable,UniProtTable_GO,UniProtTable_KEGG,UniProtTable_REACTOME, 
             file ="Ontology_Tables/UniProtTables_2024_01_23.rda")}



enrich_c4 <- proteinminion::Enrich_EASE(query = fData(clus_4)$uniprot, 
                           universe = fData(m1)$uniprot)


enrich_c4_filt <- proteinminion::enrichment_filter(enrich_c4, 
                                                   Test_p_type = "p",
                                                   Test_p = 0.05, 
                                                   foldchange = 1) %>% 
    arrange(Test_p)

enrich_c4_filt2 <- enrich_c4_filt %>% 
    mutate(Uniprot_Accession = Proteins_in_query) %>% 
    separate_longer_delim(Uniprot_Accession, delim = ";") %>% 
    left_join(UniProtTable, by = c("Uniprot_Accession" = "Uniprot_Accession")) %>% 
    select(-c(Organism:Reactome_List))


# writexl::write_xlsx(x = enrich_c4_filt, path = "RD6_3-PAM_cluster_4_enrichment_results.xlsx")
# writexl::write_xlsx(x = enrich_c4_filt2, path = "RD6_3-PAM_cluster_4_enrichment_results_long.xlsx")



## Dominant cell type ----

cell_counts <- readRDS("RD4_3-RNA-Seq_counts_all_cell_types.rds") %>% 
    filter(!is.na(celltype)) %>% 
    column_to_rownames(var = "celltype") %>% 
    t()


# colnames(cell_counts) <- cell_counts[1,]
cell_counts <- as_tibble(cell_counts[-1,], rownames = "uniprot") %>% 
    mutate(across(!matches("uniprot"), as.integer))



fData(clus_4) <- fData(clus_4) %>% 
    rownames_to_column() %>% 
    left_join(cell_counts, by = "uniprot") %>% 
    column_to_rownames("rowname")


fData(clus_4) %>% 
    select(-c(uniprot, GENEID, cluster)) %>% 
    colMeans(na.rm = T)



c4 <- fData(clus_4) %>% 
    select(-c(uniprot, GENEID, cluster)) %>% 
    as.matrix()

c4 <- sweep(c4,
            MARGIN = 1, 
            STATS = apply(c4, 1, \(x){ if_else(sum(x) == 0, 0.1, max(x, na.rm = TRUE)) }), 
            FUN = "/")

cluster_expression <- colMeans(c4, na.rm = T) %>% 
    as_tibble_row()

# writexl::write_xlsx(cluster_expression, path = "RD6_4-PAM_Cluster_4_avg_cell_type_expression_by_RNASeq_counts.xlsx")


