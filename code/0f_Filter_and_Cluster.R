#############################################################
###Program name: 0e_Filter_and_Cluster.R
###Created by:  Originally Sarah Ryan. Edited by
###             Nichole E Carlson, and William Lippitt
###Purpose: This code trims the radiomic
###         panel using a decorrelation filter and clusters the radiomic 
###         measures using robust and sparse k-means.
###         
###         The output of this code are cluster labels and a dataset called 
###         dat_clean_add_clusters.rds in the CIDA drive. It also outputs a 
###         short file of the names of the 5 most important variables 
###         according to the RSKC weights, saved as important_vars_radiomics_paper.rds
###################################################################

############################################# Libraries
library(tidyverse) #Data Manipulation
library(cluster) 
library(sparcl) #Sparse K-means
library(factoextra)
library(RSKC) #Robust and Sparse K-Means
library(ggfortify)
library(mclust)

############################################# Paths
##Location of your data
path_data <- "/DataLibrary/"
path_rad <- paste0(path_data,"radiomics_data/")
##Location of code files and output files
path_repo <- "/radiomics/"

############################################# Helper functions
path_functions <- paste0(path_repo,"code/functions/")
source(paste0(path_functions,'f_featsel.R')) # decorrelation filter itself
source(paste0(path_functions,'greedy_then_random.R')) # greedy rank function for decorrelation filter
source(paste0(path_functions,'f_clean.R')) # helper function for f_featsel
source(paste0(path_functions,'RSKC.std.permute.R')) # Calculate standardized gap stat, return max for given k, wbounds, permutation
rm(path_functions)

############################################# Data
dat <- read_csv(paste0(path_rad,"dat_clean.csv")) %>% # Primary dataset
  dplyr::select(-BDDLCO) %>%
  rename(BDDLCO = BDDLCOcorAltitude)
rad <- dat %>%  # Subset to radiomic features and compute correlation
  select(c(mask2.first_order.Mean:mask1.glcm.ep_b16_d1_mean.Uniformity))
rad_cor <- cor(rad)

############################################# Filter
set.seed(110421)
## f_featsel: returns indices of features selected using tolerance,
##            specified rank function, and a priori features to keep
## f_featsel_firstorder: applies f_featsel to firstorder features only.
##            then applies f_featsel to all features, adding selected 
##            first order features to be kept a priori. first order 
##            features are identified from fullcor column names by
##            string detection for "first_order"
## greedy_then_random: ranks candidate features according as the number 
##            of candidate features it represents. returns index of the 
##            feature with maximal rank (ties broken at random) as next 
##            feature to be kept by f_featsel
filtered_columns <- f_featsel_firstorder(fullcor=rad_cor, tol=0.9, 
                                         rank_f=greedy_then_random)  
filtered_rad <- rad[,filtered_columns]

############################################# Cluster
### Define parameters
nperm = 100 # number of reference permutated data sets for computing gap stat
nstart = 100 # initializations of clustering algorithm
K = 2:8 # number of groups considered
wbounds_candidates <- seq(1.5,floor(2*sqrt(dim(filtered_rad)[2]))/2,.5)
alpha = 0.1 # 10% trimmed mean

### Cluster for each combo of K, wbounds
mywbounds <- sapply(1:length(K),function(ki){ # takes about 10 min per ki
  print(K[ki])
  km.perm <- RSKC.std.permute(x = as.matrix(filtered_rad), K=K[ki], 
                              wbounds = wbounds_candidates, 
                              alpha = alpha, nperms=nperm, silent = TRUE) 
  saveRDS(km.perm, paste0(paste0(path_repo,'results/RSKC_perm_k',K[ki],'.RDS')))
  return(km.perm$bestw)
})
saveRDS(mywbounds, paste0(path_repo,'results/mywbounds_RSKCstd_filtered.RDS'))  

### Choose optimal combo of K, wbounds
gapStat <- lapply(K, function(ki){
  km.perm <- readRDS(paste0(path_repo,'results/RSKC_perm_k',ki,'.RDS'))
  km.perm$gaps
}) %>% do.call(rbind.data.frame, .)
colnames(gapStat) <- paste0('wbound_', wbounds_candidates)
rownames(gapStat) <- paste0('K_', K)
gapStat <- as.matrix(gapStat)
max.ind <- which(gapStat == max(gapStat), arr.ind = TRUE)
wbound_optimal <- colnames(gapStat)[max.ind[2]] %>%
  str_replace("wbound_","") %>%
  as.numeric()
k_optimal <- rownames(gapStat)[max.ind[1]] %>%
  str_replace("K_","") %>%
  as.numeric()
saveRDS(list(k = k_optimal,wbound = wbound_optimal),
        paste0(path_repo,'results/RSKC_optimalKandWbound.RDS'))

### Cluster optimally
set.seed(918273645)
optimal_cluster <- RSKC(filtered_rad, ncl = k_optimal, L1 = wbound_optimal, 
                  alpha = alpha, nstart = nstart)
rad$rskc_cluster <- optimal_cluster$labels
dat$rskc_cluster <- optimal_cluster$labels

### Retrieve features with highest weights for further analysis 
important_vars <- tail(optimal_cluster$weights[order(optimal_cluster$weights)],n=5)
names_important_vars <- names(important_vars)
# first first order measures is ___th most important
optimal_cluster$weights %>%
  sort(decreasing=T) %>%
  names() %>%
  str_which("first") %>%
  head(n=1)

############################################# Save
write_rds(names_important_vars,
          paste0(path_rad,"important_vars_radiomics_paper.rds"))

write_rds(dat,paste0(path_rad,"dat_clean_add_clusters.rds"))
write_csv(dat,paste0(path_rad,"dat_clean_add_clusters.csv"))

