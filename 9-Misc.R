library(MSnSet.utils)
library(tidyverse)
library(correlation)
library(see)


# Task 1 ------------------------------------------------------------------

# Goal - 
# want to make correlation heatmaps for INS & GCG peptides


# peptide-level data
d <- read_tsv("../../../Protein_Data/MS_output/unprocessed_fragpipe_raw_output/combined_modified_peptide.tsv")


# only msnset data
d2 <- d %>% 
    select(!matches("Spectral Count")) %>% 
    select(!matches("MaxLFQ")) %>% 
    filter(grepl("INS|GCG", Gene)) %>% 
    column_to_rownames("Modified Sequence") %>% 
    arrange(`Peptide Sequence`)

fdata <- d2 %>% 
    select(!contains("Intensity"))

edata <- d2 %>% 
    select(contains("Intensity"))

m <- MSnSet(as.matrix(edata), fdata)

# process missingness
exprs(m)[near(exprs(m), 0)] <- NA
plotNA(m)
m <- filter_by_occurrence(m, 0.01)

# log2 and center
m <- log2_zero_center(m)



## correlation heatmaps ----

# INS
m_ins <- m[grepl("INS", fData(m)$Gene), ]
c <- correlation(as.data.frame(t(exprs(m_ins))))
s <- summary(c)
p_ins <- plot(s) +
    ggtitle("INS Modified Peptides")+
    scale_fill_gradientn(colours = c("blue", "white", "red")
                         , limits = c(-1, 1))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave("RD9_1-INS_observed_peptides_correlation_heatmap.png"
#        , plot = p_ins
#        , width = 8
#        , height = 6
#        , units = "in"
#        , dpi = 200)

# GCG
m_gcg <- m[grepl("GCG", fData(m)$Gene), ]
c <- correlation(as.data.frame(t(exprs(m_gcg))))
s <- summary(c)
p_gcg <- plot(s) +
    ggtitle("GCG Modified Peptides")+
    scale_fill_gradientn(colours = c("blue", "white", "red")
                         , limits = c(-1, 1))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)
          , axis.text = element_text(size = 12))

# ggsave("RD9_1-GCG_observed_peptides_correlation_heatmap.png"
#        , plot = p_gcg
#        , width = 15
#        , height = 15
#        , units = "in"
#        , dpi = 200)


## heatmaps ----

# INS
feat_fraction <- length(featureNames(m_ins)) * 0.6
to_keep <- apply(exprs(m_ins), 2, \(x){sum(!is.na(x)) > feat_fraction })
m_ins2 <- m_ins[, to_keep]

m_ins3 <- filter_by_occurrence(m_ins, 0.1)
pData(m_ins3)$INS_GFFYTPK <- exprs(m_ins3)[featureNames(m_ins) == "GFFYTPK", ]

png("RD9_1-INS_observed_peptides_heatmap.png"
    , height = 5
    , width = 10
    , units = "in"
    , res = 200)

image_msnset(m_ins3, sOrderBy = "INS_GFFYTPK")

dev.off()

# complex_heatmap(m_ins3)
# too much missingness.


# GCG
feat_fraction <- length(featureNames(m_gcg)) * 0.7
to_keep <- apply(exprs(m_gcg), 2, \(x){sum(!is.na(x)) > feat_fraction })
m_gcg2 <- m_gcg[, to_keep]

m_gcg3 <- filter_by_occurrence(m_gcg2, 0.7)
image_msnset(m_gcg3)

png("RD9_1-GCG_observed_peptides_heatmap.png"
    , height = 5
    , width = 10
    , units = "in"
    , res = 200)

complex_heatmap(m_gcg3
                , show_row_names = T
                , color_range = c(-2, 0, 4)
                , heatmap_title = "GCG Observed Peptides")

dev.off()










