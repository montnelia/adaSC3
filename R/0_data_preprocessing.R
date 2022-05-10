require("SingleCellExperiment")
require("diffusionMap")
require("SC3")
require("mclust")
require("scater")
require("vegan") 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
biocLite("scater")
BiocManager::install(("destiny"))
biocLite("destiny")


library(destiny) # Diffusion Maps optimal sigma 
library(scater)
require(SC3)
require(SingleCellExperiment)

#############
# Read data # 
#############

biase <- readRDS("inst/biase.Rds")
goolam <- readRDS("inst/goolam.Rds")

data_list_sc3_ok <- list("biase" = logcounts(biase),
                      "goolam" = logcounts(goolam))

remove_cat_biase <- which(biase$cell_type1=="blast")
data_list_sc3_ok[[1]] <- data_list_sc3_ok[[1]][,-remove_cat_biase]

###################
### Gene filter ###
###################
gene_filter_alg_sc3 <- function(data){
  
  data_sc3_org_alg <- data
  
  for(i in 1:length(data)){
    dropouts <- rowSums((data[[i]])==0)/ncol(data[[i]])*100
    gene_filter_2 <- dropouts < 90 & dropouts > 10 
    data_sc3_org_alg[[i]] <- data[[i]][which(gene_filter_2 ==TRUE),]
  }
  data_sc3_org_alg
  #saveRDS(data_sc3_org_alg, "tmp/data_diverse_gene_filter/data_filter_sc3_R.Rds")
}
data_filter_sc3_R <- gene_filter_alg_sc3(data_list_sc3_ok)




gene_filter_paper_sc3 <- function(data){
  data_neu <- data
  remove_gene <- c()
  for(i in 1:length(data)){
    for(j in 1:dim(data[[i]])[1]){
      if(sum(data[[i]][j,] > log(2+1))/ dim(data[[i]])[2] < 0.06){
        remove_gene <- c(remove_gene, j) 
      }
      if(sum(data[[i]][j,] > log(0+1))/ dim(data[[i]])[2] >= 0.94){
        remove_gene <- c(remove_gene, j) 
      }
      
    }
    data_neu[[i]] <- data[[i]][-(remove_gene), ]
  }
  data_neu
  #saveRDS(data_neu, "tmp/data_diverse_gene_filter/data_filter_paper.Rds")
}
data_filter_paper <- gene_filter_paper_sc3(data_list_sc3_ok)










########################################
### Calculation of distance matrices ###
########################################

#################
### All genes ###
#################

for(i in 1:length(data_list_sc3_ok)){
  data_ok <- data_list_sc3_ok[[i]]
  data_t <- t(data_ok)
  exp_count_data <- exp(data_t) - 1

  D_euclidean <- dist(exp_count_data, method = "euclidean")
  D_eucl <- as.data.frame(as.matrix(D_euclidean))
  D_pearson <- 1 - cor(data_ok, use = "everything", method = "pearson")
  D_spearman <- 1 - cor(data_ok, use = "everything", method = "spearman")
  
  saveRDS(D_pearson, paste0("tmp/all_genes/pear/D_pear_", names(data_list_sc3_ok)[i], ".Rds"))
  saveRDS(D_spearman, paste0("tmp/all_genes/spear/D_spear_", names(data_list_sc3_ok)[i], ".Rds"))
  saveRDS(D_eucl, paste0("tmp/all_genes/eukl_Distanz/D_eucl_", names(data_list_sc3_ok)[i], ".Rds"))
}


##################################################
### Filter from the paper - own implementation ###
##################################################


for(i in 1:length(data_list_sc3_ok)){


  data_ok <- data_filter_paper[[i]]
  data_t <- t(data_ok)
  exp_count_data <- exp(data_t) -1
  
  D_euclidean <- dist(exp_count_data, method = "euclidean")
  D_eucl <- as.data.frame(as.matrix(D_euclidean))
  D_pearson <- 1 - cor(data_ok, use = "everything", method = "pearson")
  D_spearman <- 1 - cor(data_ok, use = "everything", method = "spearman")
  
  saveRDS(D_pearson, paste0("tmp/filter_paper/pear/D_pear_", names(data_filter_paper)[i], ".Rds"))
  saveRDS(D_spearman, paste0("tmp/filter_paper/spear/D_spear_", names(data_filter_paper)[i], ".Rds"))
  saveRDS(D_eucl, paste0("tmp/filter_paper/eukl_Distanz/D_eucl_", names(data_filter_paper)[i], ".Rds"))
}


#####################################
### Filter of origina SC3 Package ###
#####################################
data_filter_sc3_R[[1]]
for(i in 1:length(data_list_sc3_ok)){
  

  data_ok <- data_filter_sc3_R[[i]]
  data_t <- t(data_ok)
  
  exp_count_data <- exp(data_t)-1
  
  D_euclidean <- dist(exp_count_data, method = "euclidean")
  D_eucl <- as.data.frame(as.matrix(D_euclidean))
  D_pearson <- 1 - cor(data_ok, use = "everything", method = "pearson")
  D_spearman <- 1 - cor(data_ok, use = "everything", method = "spearman")
  
  saveRDS(D_pearson, paste0("tmp/filter_sc3_R/pear/D_pear_", names(data_list_sc3_ok)[i], ".Rds"))
  saveRDS(D_spearman, paste0("tmp/filter_sc3_R/spear/D_spear_", names(data_list_sc3_ok)[i], ".Rds"))
  saveRDS(D_eucl, paste0("tmp/filter_sc3_R/eukl_Distanz/D_eucl_", names(data_list_sc3_ok)[i], ".Rds"))
}


SC3_function <- function(data, k){
  data_prep <- SingleCellExperiment(assay=list(counts=exp(data)-1, logcounts=data)) #needs matrix
  rowData(data_prep)$feature_symbol <- rownames(data_prep)
  sc3_clust <- sc3(data_prep, ks = k, gene_filter=TRUE, biology = FALSE, rand_seed = 3)
  return(sc3_clust)
  
}

SC3_org <- function(data, k){
  SC3_function(as.matrix(data), k)
}

remove_cat_biase <- which(biase$cell_type1=="blast")
data_list_sc3_ok[[1]] <- data_list_sc3_ok[[1]][,-remove_cat_biase]
data_sc3_1 <- SC3_org(data_list_sc3_ok[[1]], 3)
truth_biase<- readRDS(paste0("tmp/truth/truth_biase.Rds"))  
adjustedRandIndex(data_sc3_1$sc3_3_clusters,as.factor(truth_biase)) # 0.95


