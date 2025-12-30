#Library
myPaths <- .libPaths()   # get the paths
myPaths <- c(myPaths[2], myPaths[1])  # switch them
.libPaths(myPaths)  # reassign them
library(tidyverse)
#Used library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
install.packages("ggpubr")
#1_Installation
install.packages("remotes")
remotes::install_github("RubD/Giotto")
library(Giotto)


#2_Environment setup
installGiottoEnvironment()


#3_Downloading dataset
data_directory<-"C:/Users/kwon421/Documents/Data/R/220829_HubMAP_3D_Giotto/data"
save_directory<-"C:/Users/kwon421/Documents/Data/R/220829_HubMAP_3D_Giotto/output"


#4_Creating a Giotto object
my_instructions <- createGiottoInstructions(save_plot = TRUE,
                                            show_plot = TRUE,
                                            return_plot = FALSE,
                                            save_dir = save_directory)

my_giotto_object<-createGiottoObject(raw_exprs=paste0(data_directory,
                                                      "/expression.csv"),
                                     spatial_locs = paste0(data_directory,
                                                           "/spatial_locs.csv"),
                                     instruction=my_instructions)


spatPlot3D(gobject = my_giotto_object,
           axis_scale = "real",
           point_size = 5)

###############################################
pDataDT(my_giotto_object)
fDataDT(my_giotto_object)

my_giotto_image <- createGiottoImage(gobject = my_giotto_object, 
                                     mg_object = paste0(data_directory,"/i6_image.png"))

my_giotto_object <- addGiottoImage(gobject = my_giotto_object,images = list(my_giotto_image))


my_giotto_object <- updateGiottoImage(gobject = my_giotto_object, 
                                      image_name = "image",
                                      xmax_adj = 25,
                                      xmin_adj = 25,
                                      ymax_adj = 25,
                                      ymin_adj = 25)
spatPlot2D(gobject = my_giotto_object,
           show_image = TRUE,
           point_alpha = 0.5)


#5_Preprocessing
##5-1_Filtering
filterCombinations(gobject = my_giotto_object,
                   expression_thresholds = c(1E+4, 1E+5),
                   gene_det_in_min_cells = c(3, 5, 10),
                   min_det_genes_per_cell = c(1000, 2000, 3000))
?filterCombinations
##expression_thresholds	;all thresholds to consider a gene expressed
##gene_det_in_min_cells	;minimum number of cells that should express a gene to consider that gene further
##min_det_genes_per_cell	;minimum number of expressed genes per cell to consider that cell further





my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                                 expression_threshold = 3,
                                 gene_det_in_min_cells = 3,
                                 min_det_genes_per_cell = 100)


#5-2_Normalization
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,
                                    norm_methods = "standard",
                                    scalefactor = 6000,
                                    scale_order = "first_genes")
?normalizeGiotto()
#5-3_Statistics
my_giotto_object <- addStatistics(gobject = my_giotto_object)
# view gene and cell stats respectively
head(fDataDT(my_giotto_object))
head(pDataDT(my_giotto_object))

###number of genes
spatPlot2D(gobject = my_giotto_object,
           show_image = TRUE,
           point_alpha = 1,
           point_size = 5,
           cell_color = 'nr_genes',
           color_as_factor = F)

spatPlot3D(gobject = my_giotto_object,
           axis_scale = "cube",
           point_size = 5,
           cell_color = 'nr_genes')

?spatPlot3D

#6_Clustereing and cell-type identification #HVG

##6-1_Feature selection
my_giotto_object <- calculateHVG(gobject = my_giotto_object,
                                 expression_values = "normalized",
                                 method = "cov_groups",
                                 nr_expression_groups = 20,
                                 zscore_threshold = 1.5)

##6-2_Dimensionality reduction
my_giotto_object <- runPCA(gobject = my_giotto_object,
                           expression_values = "normalized",
                           genes_to_use = "hvg")



plotPCA(gobject = my_giotto_object,
        point_size=5,
        show_legend=TRUE,
        show_plot=TRUE,
        legend_text=10)
?runPCA()
?plotPCA()
screePlot(gobject = my_giotto_object,
          expression_values = "normalized",
          genes_to_use = "hvg",
          ncp = 20,
          ylim = c(0,12.5))

?screePlot()

my_giotto_object<-runUMAP(gobject = my_giotto_object,
                          dimensions_to_use = 1:10,
                          n_neighbors = 40,
                          n_components = 2,
                          min_dist = 0.01)


plotUMAP(my_giotto_object,
         point_size=5,
         label_size=5)

?plotUMAP()

##6-3_sNN
my_giotto_object<-createNearestNetwork(gobject = my_giotto_object,
                                       dimensions_to_use = 1:10)

##6-4_Leiden cluster
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    name = "leiden_clus")
?doLeidenCluster()
plotUMAP(gobject = my_giotto_object,
         cell_color = 'leiden_clus',
         point_size = 5)

plotPCA(gobject = my_giotto_object,
        point_size=5,
        show_legend=TRUE,
        show_plot=TRUE,
        legend_text=10,
        cell_color='leiden_clus',
        color_as_factor = TRUE)

plotPCA(gobject = my_giotto_object,
        point_size=5,
        show_legend=TRUE,
        show_plot=TRUE,
        legend_text=10,
        cell_color='P00441_SOD1',
        color_as_factor = FALSE)



write.table(my_giotto_object@dimension_reduction[["cells"]][["pca"]][["pca"]][["coordinates"]], file="PCA.txt", row.names = TRUE, sep = "\t")
view(my_giotto_object@dimension_reduction[["cells"]][["pca"]][["pca"]][["coordinates"]])

write.table( my_giotto_object@cell_metadata, file="PCA.txt", row.names = TRUE, sep = "\t")
write.table(my_giotto_object@dimension_reduction[["cells"]][["pca"]][["pca"]][["misc"]][["loadings"]], file="PCA_loading.txt", row.names = TRUE, sep = "\t")
view(my_giotto_object@dimension_reduction[["cells"]][["pca"]][["pca"]][["misc"]][["loadings"]])
view( my_giotto_object@cell_metadata)

?plotPCA()

?doLeidenCluster()
##6-5_DEGs, Heatmap
ST_scran_markers_subclusters = findMarkers_one_vs_all(gobject = my_giotto_object,
                                                      method = 'scran',
                                                      expression_values ='normalized',
                                                      cluster_column = 'leiden_clus')
?findMarkers_one_vs_all()

unique(my_giotto_object)
ST_top3genes = ST_scran_markers_subclusters[, head(.SD, 3), by = 'cluster']$genes
ST_topNgenes = ST_scran_markers_subclusters[, head(.SD, 30), by = 'cluster']$genes

plotMetaDataHeatmap(gobject = my_giotto_object,
                    selected_genes = ST_topNgenes,
                    metadata_cols = c('leiden_clus'),
                    gradient_color = c("#440154", "#21908C", "#FDE725"),
                    custom_cluster_order = c("1", "3", "4", "2"))


##############################################################
#Create metagenes from cluster modules and visualize_Leiden###
###############################################################

top40_per_leiden=cluster_genes_DT[,head(.SD, 40), by = clus]
topN_per_leiden=ST_scran_markers_subclusters[,head(.SD,40), by='cluster']

Leiden_cluster_genes=topN_per_leiden$cluster; names(Leiden_cluster_genes)= topN_per_leiden$genes

my_giotto_object= createMetagenes(my_giotto_object,
                                  gene_clusters = Leiden_cluster_genes,
                                  name = 'Leiden_cluster_metagene')


spatCellPlot(my_giotto_object,
             spat_enr_names='Leiden_cluster_metagene',
             cell_annotation_values=as.character(c(1:8)),
             point_size=5,
             cow_n_col=2)

write.table(my_giotto_object@spatial_enrichment[["Leiden_cluster_metagene"]], file="ST_scran_metagene.txt", row.names = TRUE, sep = "\t")
Leiden_meta<-read.delim("clipboard")
Leiden_meta%>%
  group_by(., Z)%>%
  ggplot(aes(x=X, y= Y, fill=Cluster4))+
  geom_tile(color="white", width=1, height=1)+
  scale_fill_viridis_c(option = "viridis")+
  theme_minimal()+
  theme(axis.title = element_blank(),axis.text = element_blank(),
        strip.text = element_blank(),axis.ticks = element_blank())+
  coord_fixed()+facet_wrap(~PCA_grad$Z, ncol = 1)



                    
?plotMetaDataHeatmap
view(ST_top3genes)

#7_Spatial structure analysis tools
##7_1_Spatial grid
my_giotto_object<-createSpatialGrid(gobject = my_giotto_object,
                                    sdimx_stepsize = 50,
                                    sdimy_stepsize = 50,
                                    sdimz_stepsize = 50,
                                    minimum_padding = 1)

spatPlot(gobject = my_giotto_object,
         cell_color = 'leiden_clus',
         point_size = 5,
         point_alpha = 0.8,
         show_grid = T,
         show_image= F,
         grid_color = 'grey',
         spatial_grid_name = 'spatial_grid')
spatPlot3D(gobject = my_giotto_object,
           cell_color = 'leiden_clus',
           point_size = 20,
           axis_scale = 'cube')

spatPlot2D(gobject = my_giotto_object,
           point_size = 5,
           point_alpha = 0.5,
           coord_fix_ratio = 1,
           cell_color = 'leiden_clus')
##7-2_Spatial network
plotStatDelaunayNetwork(gobject = my_giotto_object,
                        method = c("delaunayn_geometry"),
                        maximum_distance = 100,
                        dimensions = 3)

?plotStatDelaunayNetwork()

sessionInfo()
install.packages('geometry')
?createSpatialNetwork()
my_giotto_object <- createSpatialNetwork(gobject = my_giotto_object,
                                         minimum_k = 8,
                                         maximum_distance_delaunay = 1,
                                         delaunay_method = "delaunayn_geometry")

?spatPlot3D()
spatPlot3D(
  gobject = my_giotto_object,
  sdimx = "sdimx",
  sdimy = "sdimy",
  sdimz = "sdimz",
  spat_enr_names = NULL,
  point_size = 10,
  cell_color = "leiden_clus",
  cell_color_code = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 10,
  other_cell_alpha = 0.5,
  show_network = T,
  spatial_network_name = "Delaunay_network",
  network_color = "black",
  network_alpha = 0.5,
  show_grid = T,
  spatial_grid_name = "spatial_grid",
  grid_color = NULL,
  grid_alpha = 1,
  title = "",
  show_legend = T,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spat3D"
)

spatPlot(gobject = my_giotto_object,
         show_network = T,
         show_image =F,
         point_shape = "no_border",
         network_color = 'black',
         spatial_network_name = 'Delaunay_network',
         point_size = 3,
         cell_color = "leiden_clus",
         coord_fix_ratio = 1)
?spatPlot()




##7-3_Spatial gene expression pattern
binspect_k <- binSpect(gobject = my_giotto_object,
                       bin_method = "kmeans",
                       expression_values = "normalized",
                       spatial_network_name = "Delaunay_network")
?binSpect()
interesting_genes = ST_topNgenes


spatGenePlot(gobject = my_giotto_object,
             expression_values = "normalized",
             genes = ST_top3genes,
             cow_n_col = 2,
             point_size = 2.5)


##7-4_Spatial gene co-expression modules
ext_spatial_genes = binspect_k[1:1000]$genes
spat_cor_netw_DT = detectSpatialCorGenes(my_giotto_object,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_genes = ext_spatial_genes)
?detectSpatialCorGenes()
?clusterSpatialCorGenes()

spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT,
                                          name = 'spat_netw_clus',
                                          k=3)
                                          

heatmSpatialCorGenes(gobject = my_giotto_object,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus')



cluster_genes_DT= showSpatialCorGenes(spat_cor_netw_DT,
                                      use_clus_name = 'spat_netw_clus',
                                      show_top_genes = 1)

top40_per_module=cluster_genes_DT[,head(.SD, 40), by = clus]


#Create metagenes from cluster modules and visualize
cluster_genes=top40_per_module$clus; names(cluster_genes)= top40_per_module$gene_ID
my_giotto_object= createMetagenes(my_giotto_object,
                                  gene_clusters = cluster_genes,
                                  name = 'cluster_metagene')

?createMetagenes()
spatCellPlot(my_giotto_object,
             spat_enr_names='cluster_metagene',
             cell_annotation_values=as.character(c(1:8)),
             point_size=5,
             cow_n_col=2)
spatDimPlot3D()
write.table(cluster_genes_DT, file="cluster_1K.txt", row.names = FALSE, sep = "\t")
#8_Spatial domain detection by using a HMRF
##8_1_Implementation
hmrf_folder <- paste0(path.expand(save_directory),'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

# input is the top 40 genes per co-expression module
remotes::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='master')
HMRF_spat_genes = doHMRF(gobject = my_giotto_object,
                         expression_values = "normalized",
                         spatial_genes = names(cluster_genes),
                         spatial_network_name = "Delaunay_network",
                         zscore = "none",
                         k = 8,
                         betas = c(0,5,6),
                         output_folder = paste0(hmrf_folder, '/', 'HMRF_output2'))

## add HMRF of interest to giotto object 
my_giotto_object = addHMRF(gobject = my_giotto_object,
                           HMRFoutput = HMRF_spat_genes,
                           k = 8,
                           betas_to_add = c(0, 10, 15, 20),
                           hmrf_name = 'HMRF')

spatPlot2D(my_giotto_object,
           cell_color = 'HMRF_k8_b.0',
           show_image = TRUE,
           point_size = 4.75,
           coord_fix_ratio = 1)
spatPlot2D(my_giotto_object,
           cell_color = 'HMRF_k8_b.20',
           show_image = TRUE,
           point_size = 4.75,
           coord_fix_ratio = 1)

spatPlot3D(
  my_giotto_object,
  sdimx = "sdimx",
  sdimy = "sdimy",
  sdimz = "sdimz",
  spat_enr_names = NULL,
  point_size = 3,
  cell_color = NULL,
  cell_color_code = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 0.5,
  other_cell_alpha = 0.5,
  show_network = F,
  spatial_network_name = "Delaunay_network",
  network_color = NULL,
  network_alpha = 1,
  show_grid = F,
  spatial_grid_name = "spatial_grid",
  grid_color = NULL,
  grid_alpha = 1,
  title = "",
  show_legend = T,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spat3D"
)
#9_Spatial Proximity
##9-1_Cell proximity enrichment
cell_proximities = cellProximityEnrichment(gobject = my_giotto_object,
                                           cluster_column = 'leiden_clus',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)


cellProximityBarplot(gobject = my_giotto_object,
                     CPscore = cell_proximities,
                     min_orig_ints = 3,
                     min_sim_ints = 3)

cellProximityNetwork(gobject = my_giotto_object,
                     CPscore = cell_proximities,
                     remove_self_edges = F,
                     self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1,2),
                     edge_weight_range_enrichment = c(2,5))

##9-2_Interaction-changed genes (ICGs)
## select top 25th highest expressing genes 
gene_metadata = fDataDT(my_giotto_object) 
high_expressed_genes = gene_metadata[mean_expr_det > quantile(gene_metadata$mean_expr_det)[4]]$gene_ID 

## identify genes that are associated with proximity to other cell types 


ICGscoresHighGenes = findInteractionChangedGenes(gobject = my_giotto_object,
                                                 selected_genes = high_expressed_genes,
                                                 spatial_network_name = 'Delaunay_network',
                                                 cluster_column = 'HMRF_k8_b.20',
                                                 diff_test = 'permutation',
                                                 adjust_method = 'fdr',
                                                 nr_permutations = 2000,
                                                 do_parallel = TRUE,
                                                 cores = 4)

## visualize 
plotCellProximityGenes(my_giotto_object,cpgObject = ICGscoresHighGenes,method = 'dotplot')

## filter genes 
ICGscoresFilt = filterInteractionChangedGenes(ICGscoresHighGenes)

##===============================================================================
## visualize subset of interaction changed genes (ICGs)



ICGscoresFilt$CPGscores[type_int == 'hetero'][cell_type == '3']
ICG_genes = c('LAMC2', 'CXCL10', 'PIP', 'WIPI2', 'PI16')
ICG_genes_types = c('7', '7', '2', '2', '5')
names(ICG_genes) = ICG_genes_types



#8_Spatial proximity-associatede cell-cell interactions
cell_type_vector <- readRDS(paste0(data_directory,"/", "cell_type_vector.RDS"))
cell_type_vector