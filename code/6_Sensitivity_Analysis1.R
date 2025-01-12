#############################################################
###Program name: 6_Sensitivity_Analysis1_WilliamsCopy.R
###Created by: William Lippitt, Editing from Zachary Buchalski's work
###Purpose: A sensitivity analysis for the primary analysis. We assess 
###         sensitivity of cluster labels and model significance to 
###         randomness in the decorrelation filter, in the cluster fitting
###         initializations, and in the permutations used in estimating the
###         L1 parameter and number of clusters. 
###Requirements for running: This code file should be locally run from the 
###         Code folder of the radiomics repo. Before running, dat_clean.csv
###         and all sensitivity analysis RDS objects should be dropped into
###         the repo DataProcessed folder. All of these files can be found
###         in the p-drive project folder, in particular with all sensitivity
###         analysis RDS objects in their own folder for ease of locating.
###
###         Produces:
###         figures/Figure_SA1_ARI.png
###         results/sensitivity1_pvaluesummary.RDS
###################################################################

################################################# Libraries needed
library("parallel") #Parallelization
library("doParallel") #Parallelization
library("foreach") #Parallelization
library(tidyverse) #Data Manipulation
library(cluster) #Gap stat calculation/l1 optimization
library(RSKC) #Cluster label creation
library(mclust) #Adjusted Rand Indices

####################################### paths
##Location of your data
path_pdrive <- "/DataLibrary/radiomics_data/"
##Location of your code and output
path_repo <- "radiomics/"
path_results <- paste0(path_repo,"results/")
path_figures <- paste0(path_repo,"figures/")
path_functions <- paste0(path_repo,"code/Functions/")

################################################# Helper functions
source(paste0(path_functions,'f_featsel.R')) #correlation filter itself
source(paste0(path_functions,'greedy_then_random.R')) # rank function
source(paste0(path_functions,'f_clean.R')) # helper for featsel
source(paste0(path_functions,'RSKC.std.permute.R')) # standardized gap statistic

################################################# Data
dat <- read_csv(paste0(path_pdrive,"dat_clean.csv"))
rad <- dat %>% select(c(mask2.first_order.Mean:mask1.glcm.ep_b16_d1_mean.Uniformity))
rad_cor <- cor(rad)

################################################# Wrapper for parallelizing
simWrapper <- function(n_sim,f_sim,TF_parallel=F,
                       n_cores=2,list_export=NULL,
                       list_package=c(),
                       f_clusterCall=function(){}){
  v_seed <- sample(2^31-1,n_sim)
  if(TF_parallel){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    try({
      clusterExport(cl,list_export)
      clusterCall(cl,function() sapply(list_package,library,character.only=T))
      clusterCall(cl,f_clusterCall)
      output_full <- foreach(i =1:n_sim,.combine = rbind) %dopar% {
        set.seed(v_seed[i])
        return(f_sim(i))
      }
    })
    stopCluster(cl)
  }else{
    output_full <- foreach(i =1:n_sim,.combine = rbind) %do% {
      set.seed(v_seed[i])
      return(f_sim(i))
    }
  }
  return(output_full)
}

################################################# Function to apply decorrelation filter and fit clusters
#### We will apply this function 512 times to produce 512 analyses for comparison
f_crudeCV <- function(K=2:8,rad,rad_cor,tol=0.9,alpha=0.1,
                      nperms=10, nstart=10, rank_f=greedy_then_random){
  mywbounds <- numeric(length(K))
  bestgapstat <- numeric(length(K))
  # decorrelate
  filtered_columns <- f_featsel_firstorder(fullcor=rad_cor,tol=tol, rank_f=rank_f)
  filtered_rad <- rad[,filtered_columns]
  # 10 equally spaced values 1.5 to sqrt(num_features)
  wbounds_candidates <- seq(1.5,floor(2*sqrt(length(filtered_columns)))/2,0.5)
  # estimate number of clusters by computing gap stat for each candiate Kval
  for(ki in 1:length(K)){
    km.perm <- RSKC.std.permute(x = as.matrix(filtered_rad), K=K[ki],
                                wbounds = wbounds_candidates,
                                alpha = alpha, nperms=nperms, silent = TRUE)
    mywbounds[ki] <- km.perm$bestw
    bestgapstat[ki] <- max(unlist(km.perm))
  } 
  # extract optimal K and L1 bound
  indie <- which.max(bestgapstat)
  nclust <- K[indie]
  l1_parm <- mywbounds[indie]
  # use optimal parameters to fit clusters
  cluster_ki <- RSKC(filtered_rad,ncl = nclust,L1=l1_parm,alpha=alpha,nstart=nstart)
  return(cluster_ki$labels)
}


################################################# Settings
n_sim <- 32 # This is the number of analyses run per batch. We run 32 batches
#             for a total of 512 analyses.
n_batch <- 16
f_sim <- function(i) f_crudeCV(rad=rad,rad_cor=rad_cor,rank_f=greedy_then_random)
TF_parallel <- T
n_cores=8
list_export=list("f_featsel","f_featsel_firstorder","greedy_then_random","list_package","f_crudeCV",
                 "rad_cor","rad","f_clean","RSKC.std.permute")
list_package=c("RSKC","stringr")

################################################# The 512 analyses!
#Set seeds for each loop
set.seed(122221)
seeds <- sample(2^30-1,n_batch)
#### Takes a while to run. Saved and read back in.
for(i in 1:n_batch){
  current_save_path <- # batch save location
    paste0(paste0(path_pdrive,"SensitivityAnalyses/sensitivity1results_",i,".RDS"))
  if(file.exists(current_save_path)){next} # don't run batches that have already run
  set.seed(seeds[i])
  simulated_clusters <- # run batch
    simWrapper(n_sim=n_sim,
               f_sim=f_sim,
               TF_parallel=TF_parallel,
               n_cores=n_cores,
               list_export=list_export,
               list_package=list_package)
  saveRDS(simulated_clusters,current_save_path)
  rm(simulated_clusters)
  gc()
}
#### read in
simulated_clusters <-
  do.call(
    rbind,
    lapply(
      1:n_batch,
      function(i){
        current_save_path <-
          paste0(path_pdrive,"SensitivityAnalyses/sensitivity1results_",i,".RDS")
        return(readRDS(current_save_path))
      }
    )
  )

################################################# Compute pairwise ARI
#### between analyses to assess sensitivity of cluster labels
aris_par <-
  expand.grid(i=1:512,j=1:512) %>%
  as_tibble() %>%
  filter(i<j) %>%
  mutate(
    ARI =
      apply(
        cbind(i,j),
        1,
        function(ij){
          return(
            adjustedRandIndex(
              simulated_clusters[ij[1],],simulated_clusters[ij[2],]
            )
          )
        }
      )
  )
saveRDS(aris_par,paste0(path_pdrive,"SensitivityAnalyses/sensitivity1aris.RDS"))
aris_par <- readRDS(paste0(path_pdrive,"SensitivityAnalyses/sensitivity1aris.RDS"))
png(paste0(path_figures,"Figure_SA1_ARI.png"))
hist(aris_par$ARI,breaks=30,
     xlab="Adj. Rand Index", xlim=0:1, ylim=c(0,10000),
     main="Pairwise Adj. Rand Indices from Decorrelation \n Filter and Estimating # of clusters")
dev.off()

################################################# Significance in predicting FVC
#### Compute p-values and ARI with and without scad to assess sensitivity of cluster
#### significance and assess overlap of identified clusters with Scadding
base_mod <- lm(BDFVC ~ AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat)
base_scadding <- lm(BDFVC ~ factor(SCADDING)+AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat)
sensitivity1_summary <-
  sapply(
    1:dim(simulated_clusters)[1],
    function(i){
      dat$cluster <- simulated_clusters[i,]
      full_mod <- lm(BDFVC ~ factor(SCADDING)+ factor(cluster) + AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat)
      cluster_mod <-lm(BDFVC ~ factor(cluster) + AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat)
      pscad <- (anova(full_mod, base_scadding)$`Pr(>F)`)[2]
      pbase <- (anova(base_mod,cluster_mod)$`Pr(>F)`)[2]
      ariscad <- adjustedRandIndex(dat$SCADDING,dat$cluster)
      return(
        c(
          pscad = pscad,
          pbase = pbase,
          ariscad = ariscad
        )
      )
    }
  ) %>%
  t() %>%
  as_tibble()

###fit clusters and Scadding do not agree at all. Good!
hist(sensitivity1_summary$ariscad)

###max p values are tiny. significance insensitive.
pmax_base <- max(sensitivity1_summary$pbase)
pmax_scad <- max(sensitivity1_summary$pscad)

###proportion below 0.05
prop_005_base <- sum(sensitivity1_summary$pbase < 0.05)/nrow(sensitivity1_summary)
prop_005_scad <- sum(sensitivity1_summary$pscad < 0.05)/nrow(sensitivity1_summary)

###proportion below 0.01
prop_001_base <- sum(sensitivity1_summary$pbase < 0.01)/nrow(sensitivity1_summary)
prop_001_scad <- sum(sensitivity1_summary$pscad < 0.01)/nrow(sensitivity1_summary)

sensitivity1_pvals <-
  c(
    maxariscad = max(sensitivity1_summary$ariscad),
    pmax_base = pmax_base,
    pmax_scad = pmax_scad,
    prop_001_base = prop_001_base,
    prop_001_scad = prop_001_scad
  )

saveRDS(sensitivity1_pvals,paste0(path_pdrive,"SensitivityAnalyses/sensitivity1_pvaluesummary.RDS"))


