if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("destiny")
BiocManager::install("SC3")
BiocManager::install("SingleCellExperiment")
require("destiny")
require(SC3)
require(SingleCellExperiment)
require(diffusionMap)
require(pcaMethods)
require(lle)
require(vegan) #IM
require(RANN)
require(mclust)


#############
# Read data # 
#############


getwd()
data_base <- list.files("inst")

biase <- readRDS(paste0("inst/",data_base[1]))
goolam <- readRDS(paste0("inst/",data_base[2]))


data_list_sc3_ok <- list("biase" = logcounts(biase),
                         "goolam" = logcounts(goolam))                        

remove_cat_biase <- which(biase$cell_type1=="blast")
data_list_sc3_ok[[1]] <- data_list_sc3_ok[[1]][,-remove_cat_biase]

data_list_sc3_ok <- readRDS("C:/Users/ri89por/Desktop/ECDA/R Revision/tmp/data_list_sc3_ok.RDS")
#####################################################################
### Transformation PCA, Laplacian, Diffusion Maps, Isomaps, MDS #####
#####################################################################

getwd()






names(data_list_sc3_ok)
# h = 1: biase: K = 1
# h = 5: goolam, K= 30
# h = 4: fan: k= 40
# h= 2: darmanis: k= 40
# h=3: deng: k= 40
# h =11: yan k= 30


folder_str = "tmp/filter_sc3_R/"
h = 1
name =""
kmin_dataset = 1

dataset_eukl <- list.files(paste0(folder_str, "eukl_Distanz"))
dataset_pear <- list.files(paste0(folder_str, "pear"))
dataset_spear <- list.files(paste0(folder_str, "spear"))




D_eucl <- readRDS(paste0(folder_str, "eukl_Distanz_sc3/", dataset_eukl[h]))
D_spear <- readRDS(paste0(folder_str,"spear/", dataset_spear[h]))
D_pear <- readRDS(paste0(folder_str,"pear/", dataset_pear[h]))

truth <- as.factor(readRDS(paste0("tmp/truth/",list.files("tmp/truth/")[h])))



k_true  <- length(levels(truth))
k_true

output <- data.frame(matrix(NA, nrow= 10, ncol = 17))



colnames(output) <- c( "P", "L", "D", "P + L", "P + D", "L + D", 
                       "I", "P + I", "L+ I",
                       "LLE", "P + LLE", "L + LLE",
                       "MDS", "P+ MDS", "L + MDS", "n_sc", "boot_sample")









for(l in 1:10){
  
  set.seed(l)
  
  D_eucl <- readRDS(paste0(folder_str, "eukl_Distanz_sc3/", dataset_eukl[h]))
  D_spear <- readRDS(paste0(folder_str,"spear/", dataset_spear[h]))
  D_pear <- readRDS(paste0(folder_str,"pear/", dataset_pear[h]))
  
  
  truth <- as.factor(readRDS(paste0("tmp/truth/",list.files("tmp/truth/")[h])))
  
  sc <- 1:dim(D_eucl)[2]
  
  boot_sample <- sc
  boot_sample
  #}
  n_sc_boot <- length(sc) # number of single-cells 
  
  
  
  truth <- truth[boot_sample]
  
  
  
  D_eucl <- D_eucl[boot_sample, boot_sample]
  D_spear <- D_spear[boot_sample, boot_sample]
  D_pear <- D_pear[boot_sample, boot_sample]
  
  #PCA_eucl <- prcomp(D_eucl[boot_sample,boot_sample], scale = TRUE, center = TRUE)
  PCA_eucl <- prcomp(D_eucl, scale = TRUE, center = TRUE)
  #PCA_eucl <- readRDS(paste0("tmp/filter_sc3_R/eukl_Distanz_sc3/P_eucl_", names(data_list_sc3_ok)[i], ".Rds"))
  #PCA_eucl <- PCA_eucl@metadata$sc3$transformations$euclidean_pca
  
  PCA_spear <- prcomp(D_spear, scale = TRUE, center = TRUE)
  PCA_pear <- prcomp(D_pear, scale = TRUE, center = TRUE)
  
  #Lapl_eucl <- readRDS(paste0("tmp/filter_sc3_R/eukl_Distanz_sc3/L_eucl_", names(data_list_sc3_ok)[i], ".Rds"))
  #Lapl_eucl <- Lapl_eucl@metadata$sc3$transformations$euclidean_laplacian
  Lapl_eucl <- norm_laplacian(as.matrix(D_eucl)) # no scaling no centering, bzw in funktion bereits 
  Lapl_spear <- norm_laplacian(as.matrix(D_spear)) 
  Lapl_pear <- norm_laplacian(as.matrix(D_pear)) 
  
  
  n_DC <- floor(0.04 * dim(D_spear)[2]):ceiling(0.07* dim(D_spear)[2]) 
  
  #set.seed(3)
  
  
  
  
  n_DC  <- sort(n_DC)
  
  
  Diffus_eucl <- as.data.frame(as.matrix(D_eucl))# no scaling no centering, bzw in funktion bereits 
  
  opt_sigma <- destiny::find_sigmas(Diffus_eucl, verbose =FALSE)
  Diffmap_eucl <- destiny::DiffusionMap(Diffus_eucl, n_eigs = n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma)
  opt_sigma <- destiny::find_sigmas(D_spear, verbose =FALSE)
  Diffmap_spear = destiny::DiffusionMap(D_spear, n_eigs= n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma) 
  opt_sigma <- destiny::find_sigmas(D_pear, verbose =FALSE)
  Diffmap_pear = destiny::DiffusionMap(D_pear, n_eigs= n_DC[length(n_DC)], sigma = opt_sigma@optimal_sigma) 
  

  
  knn_test <- calc_k(D_eucl,  n_DC[length(n_DC)], kmin =kmin_dataset, kmax =dim(D_eucl)[1], plotres=FALSE, parallel = TRUE, cpus = 6)
  k_neigh <- knn_test[which(knn_test$rho==min(knn_test$rho)),"k"]
  iso_eucl <- isomap(D_eucl, k=k_neigh, ndim = n_DC[length(n_DC)])
  lle_eucl <- lle(D_eucl, m= n_DC[length(n_DC)], k=k_neigh, reg=2, ss=FALSE, id=TRUE, v=0.9)
  mds_eucl <- cmdscale(D_eucl, eig=TRUE, k=n_DC[length(n_DC)]) # k is the number of dim
  
  
  
  knn_test <- calc_k(D_spear,  n_DC[length(n_DC)], kmin =kmin_dataset, kmax =dim(D_eucl)[1], plotres=FALSE, parallel = TRUE, cpus = 6)
  k_neigh <- knn_test[which(knn_test$rho==min(knn_test$rho)),"k"]
  iso_spear <- isomap(D_spear, k=k_neigh, ndim = n_DC[length(n_DC)])
  lle_spear <- lle(D_spear, m= n_DC[length(n_DC)], k=k_neigh, reg=2, ss=FALSE, id=TRUE, v=0.9)
  mds_spear <- cmdscale(D_spear, eig=TRUE, k=n_DC[length(n_DC)]) # k is the number of dim
  
  knn_test <- calc_k(D_pear,  n_DC[length(n_DC)], kmin =kmin_dataset, kmax =dim(D_eucl)[1], plotres=FALSE, parallel = TRUE, cpus = 6)
  k_neigh <- knn_test[which(knn_test$rho==min(knn_test$rho)),"k"]
  iso_pear <- isomap(D_pear, k=k_neigh, ndim = n_DC[length(n_DC)])
  lle_pear <- lle(D_pear, m= n_DC[length(n_DC)], k=k_neigh, reg=2, ss=FALSE, id=TRUE, v=0.9)
  mds_pear <- cmdscale(D_pear, eig=TRUE, k=n_DC[length(n_DC)]) # k is the number of dim
  
  
  
  
  
  res_PCA_eucl<- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_PCA_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_PCA_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  
  res_Lapl_eucl <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_Lapl_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_Lapl_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  
  res_diff_eucl<- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_diff_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  res_diff_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot))
  
  res_iso_eucl <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_iso_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_iso_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  
  
  res_lle_eucl <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_lle_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_lle_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  
  
  res_mds_eucl <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_mds_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  res_mds_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = n_sc_boot)) # NAs entstehen wenn res 36x80
  
  tmp <- c()
  
  
  for(j in n_DC[1]:n_DC[length(n_DC)]){ #pc
    
    tmp <- kmeans(PCA_eucl$rotation[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(PCA_pear$rotation[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(PCA_spear$rotation[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_PCA_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    
    tmp <- kmeans(eigen(Lapl_eucl)$vectors[, order(eigen(Lapl_eucl)$values)][,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(eigen(Lapl_pear)$vectors[, order(eigen(Lapl_pear)$values)][,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(eigen(Lapl_spear)$vectors[, order(eigen(Lapl_spear)$values)][,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_Lapl_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    
    
    tmp <- kmeans(Diffmap_eucl@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_diff_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(Diffmap_pear@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_diff_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(Diffmap_spear@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_diff_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    
    
    # Individual k-means clustering based on Isomaps 
    tmp <- kmeans(iso_eucl$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_iso_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(iso_pear$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_iso_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(iso_spear$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_iso_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    # Individual k-means clustering based on LLE 
    tmp <- kmeans(lle_eucl$Y[,(n_DC[length(n_DC)]-j):n_DC[length(n_DC)]], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_lle_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(lle_pear$Y[,(n_DC[length(n_DC)]-j):n_DC[length(n_DC)]], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_lle_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(lle_spear$Y[,(n_DC[length(n_DC)]-j):n_DC[length(n_DC)]], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_lle_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    # Individual k-means clustering based on MDS
    tmp <- kmeans(mds_eucl$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_mds_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(mds_pear$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_mds_pear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
    tmp <- kmeans(mds_spear$points[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
    res_mds_spear[n_DC[length(n_DC)]-j+1, ] <- tmp$cluster
    
  }
  
  ###########################################################################################
  #### Create combinations that are targeted to be determined with regard to performance ####
  ###########################################################################################
  
  #### MDS -> ari_mds #### 
  neu_diff <- rbind(res_mds_eucl, res_mds_pear, res_mds_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_mds <- adjustedRandIndex(truth, res_2)
  
  #### P + MDS -> ari_pmds #### 
  neu_diff <- rbind(res_mds_eucl, res_mds_pear, res_mds_spear, res_PCA_eucl, res_PCA_pear, res_PCA_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_pmds <- adjustedRandIndex(truth, res_2)
  
  #### L + MDS -> ari_lmds #### 
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_mds_eucl, res_mds_pear, res_mds_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_lmds <- adjustedRandIndex(truth, res_2)
  
  #### LLE -> ari_lle #### 
  neu_diff <- rbind(res_lle_eucl, res_lle_pear, res_lle_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_lle <- adjustedRandIndex(truth, res_2)
  
  #### P + LLE -> ari_plle #### 
  neu_diff <- rbind(res_PCA_eucl, res_PCA_pear, res_PCA_spear,
                    res_lle_eucl, res_lle_pear, res_lle_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_plle <- adjustedRandIndex(truth, res_2)
  
  #### L + LLE -> ari_llle #### 
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_lle_eucl, res_lle_pear, res_lle_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_llle <- adjustedRandIndex(truth, res_2)
  
  ### P + I -> ari_piso### 
  
  
  neu_diff <- rbind(res_iso_eucl, res_iso_pear, res_iso_spear,
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_piso <- adjustedRandIndex(truth, res_2)
  
  #### I -> ari_iso #### 
  
  neu_diff <- rbind(res_iso_eucl, res_iso_pear, res_iso_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_iso <- adjustedRandIndex(truth, res_2)
  
  ### P + I -> ari_piso### 
  
  
  neu_diff <- rbind(res_iso_eucl, res_iso_pear, res_iso_spear,
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_piso <- adjustedRandIndex(truth, res_2)
  
  
  ### L + I -> ari_liso### 
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_iso_eucl, res_iso_pear, res_iso_spear)
  
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_liso <- adjustedRandIndex(truth, res_2)
  
  
  
  ### L + P  + D -> ari_dipl### 
  
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_diff_eucl, res_diff_pear, res_diff_spear, 
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_dipl <- adjustedRandIndex(truth, res_2)
  
  
  
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_sc3 <- adjustedRandIndex(truth, res_2)
  
  neu_diff <- rbind(res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_p <- adjustedRandIndex(truth, res_2)
  
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear)
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_l <- adjustedRandIndex(truth, res_2)
  
  neu_diff <- rbind(res_diff_eucl, res_diff_pear, res_diff_spear)
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_d <- adjustedRandIndex(truth, res_2)
  
  
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_pl <- adjustedRandIndex(truth, res_2)
  
  neu_diff <- rbind(res_diff_eucl, res_diff_pear, res_diff_spear, 
                    res_PCA_eucl,res_PCA_pear,res_PCA_spear)
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_pd <- adjustedRandIndex(truth, res_2)
  
  neu_diff <- rbind(res_Lapl_eucl, res_Lapl_pear, res_Lapl_spear,
                    res_diff_eucl, res_diff_pear, res_diff_spear)
  
  
  sc3_diff <- rbind(neu_diff)
  cl_kmeans_dmap_pca_lapl <- na.omit(sc3_diff)
  dim(cl_kmeans_dmap_pca_lapl)
  
  ### Creation of consensus matrix 
  consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
  for(i in 1:dim(D_eucl)[2]){
    for(j in 1:dim(D_eucl)[2]){
      consens_mat[i,j]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
      consens_mat[j,i]<- sum(cl_kmeans_dmap_pca_lapl[,i]== cl_kmeans_dmap_pca_lapl[,j])
    }
  }
  
  ### Complete Linkage for final clustering decision
  ### based on the scores of previous clustering decisions
  hc <- hclust(dist(consens_mat))
  res_2 <- cutree(hc, k=k_true)
  ari_ld <- adjustedRandIndex(truth, res_2)
  ari_ges <- c(ari_p, ari_l, ari_d, ari_pl, ari_pd, ari_ld,
               ari_iso, ari_piso, ari_liso,
               ari_lle, ari_plle, ari_llle,
               ari_mds, ari_pmds, ari_lmds, n_sc_boot, l)
  
  
  
  output[l,] <- ari_ges
  print(l) 

  output
  getwd()
}
output

