#############################################################
###Program name: 6_Sensitivity_Analysis2_WilliamsCopy.R
###Created by: William Lippitt
###Purpose: A sensitivity analysis for the primary analysis. We assess 
###         sensitivity of cluster labels and model significance to 
###         randomness in the decorrelation filter, in the cluster fitting
###         initializations, and in the observations (bootstrapping). Distinct
###         from sensitivity analysis 1, we fix the number of cluster and the 
###         L1 parameter.
###Requirements for running: This code file should be locally run from the 
###         Code folder of the radiomics repo. Before running, dat_clean.csv
###         and all sensitivity analysis RDS objects should be dropped into
###         the repo DataProcessed folder. All of these files can be found
###         in the p-drive project folder, in particular with all sensitivity
###         analysis RDS objects in their own folder for ease of locating.
###
###         Produces:
###         figures/Figure_SA2_ARI.png
###         results/sensitivity2_pvalueNnbootsummary.RDS
###################################################################

################################################# Libraries needed
library(tidyverse) #Data Manipulation
library(cluster) #Gap stat calculation/l1 optimization
library(RSKC) #Cluster label creation
library(mclust) #Adjusted Rand Indices

####################################### paths
##Location of your data
path_pdrive <- "/DataLibrary/radiomics_data/"
##Location of your code and output
path_repo <- "/radiomics/"
path_results <- paste0(path_repo,"results/")
path_figures <- paste0(path_repo,"figures/")
path_functions <- paste0(path_repo,"code/Functions/")

################################################# Helper functions

# source('functions/my_fviz_gap_stat.R') 
# source('functions/clusGapKMSC.R')
# source('functions/myRSKC.R')
# source('functions/myGap.R')
source(paste0(path_functions,'f_featsel.R')) #correlation filter itself
source(paste0(path_functions,'greedy_then_random.R')) #rank function
source(paste0(path_functions,'f_clean.R')) #option in featsel
source(paste0(path_functions,'RSKC.std.permute.R')) #standardized gap statistic

################################################# Read data in
dat <- read_csv(paste0(path_pdrive,"dat_clean.csv"))
rad <- dat %>% select(c(mask2.first_order.Mean:mask1.glcm.ep_b16_d1_mean.Uniformity))
rad_cor <- cor(rad)
list_optimal <- readRDS(paste0(path_repo,'results/RSKC_optimalKandWbound.RDS'))

################################################# Settings
alpha <- 0.1
nstart <- 10
nsim <- 500   # we'll be sampling 500 bootstrap datasets

################################################# Bootstrapping!
n=dim(rad)[1]
rad$obs_labels<-1:n          # number the observations
vec_of_names <- 
  c(                         # to compare bootstraps, we need to save both
    paste0("index",1:n),     # the observation indices sampled and
    paste0("label",1:n)      # the labels fit to the observations
  )
set.seed(1)
boot_label <-
  sapply(
    1:nsim,
    function(i){
      boot_samps <- rad[sample(n,n,replace=T),]      # bootstrap
      boot_cor <- cor(boot_samps[,-567])
      filtered_boot_cols <- f_featsel_firstorder(fullcor=boot_cor,tol=0.9, rank_f=greedy_then_random)    # decorrelate
      filtered_boots <- boot_samps[,filtered_boot_cols]
      boot_cluster_i <- RSKC(filtered_boots,ncl = list_optimal$k,L1=list_optimal$wbound,alpha=alpha,nstart=nstart)     # cluster
      return_vec <- c(boot_samps$obs_labels,boot_cluster_i$labels)    # extract obs indices and cluster labels
      names(return_vec) <- vec_of_names
      return(return_vec)
    }
  ) %>%
  t() %>%
  as_tibble()
saveRDS(boot_label,paste0(path_pdrive,"SensitivityAnalyses/sensitivity2results.RDS"))
boot_label <- readRDS(paste0(path_pdrive,"SensitivityAnalyses/sensitivity2results.RDS"))

################################################# Compute pairwise ARI
#### between analyses to assess sensitivity of cluster labels. Because
#### we're bootstrapping, an overlapping subsample between analyses
#### must be computed before computing ARI.
aris <-
  expand.grid(i = 1:nsim, j = 1:nsim) %>%
  filter(i<j) %>%
  apply(
    1,
    function(ij){
      i <- ij[1]
      j <- ij[2]
      boot_i <-      # pull bootstrap, remove repetition
        matrix(as.numeric(boot_label[i,]),ncol=2) %>%
        unique()
      boot_j <-
        matrix(as.numeric(boot_label[j,]),ncol=2) %>%
        unique()
      boot_i <- boot_i[order(boot_i[,1]),]   # ordered bootstrap for ease of comparison
      boot_j <- boot_j[order(boot_j[,1]),]
      intersect_values <- intersect(boot_i[,1],boot_j[,1])   # overlapping subsample
      ari <-       # ari on overlapping subsample
        adjustedRandIndex(
          boot_i[which(boot_i[,1] %in% intersect_values),2],
          boot_j[which(boot_j[,1] %in% intersect_values),2]
        )
      return(c(ij,ari=ari))
    }
  ) %>%
  t() %>%
  as_tibble()
saveRDS(aris,paste0(path_pdrive,"SensitivityAnalyses/sensitivity2aris.RDS"))
aris <- readRDS(paste0(path_pdrive,"SensitivityAnalyses/sensitivity2aris.RDS"))
png(paste0(path_figures,"Figure_SA2_ARI.png"))
hist(
  aris$ari,breaks=30, xlim=0:1, ylim=c(0,10000),
  main="Pairwise Adj. Rand Indices from Decorrelation \n Filter and Bootstrapping Error",
  xlab="Adj. Rand Index")
dev.off()

################################################# Significance in predicting FVC
#### Compute p-values and ARI with and without scad to assess sensitivity of cluster
#### significance and assess overlap of identified clusters with Scadding
sensitivity2_summary <-
  sapply(
    1:nsim,
    function(i){
      boot_i <- 
        matrix(as.numeric(boot_label[i,]),ncol=2)
      dat_i <- dat[boot_i[,1],]
      dat_i$cluster <- boot_i[,2]
      base_mod <- lm(BDFVC ~ AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat_i)
      base_scadding <- lm(BDFVC ~ factor(SCADDING)+AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat_i)
      full_mod <- lm(BDFVC ~ factor(SCADDING)+ factor(cluster) + AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat_i)
      cluster_mod <-lm(BDFVC ~ factor(cluster) + AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT, data = dat_i)
      pscad <- (anova(full_mod, base_scadding)$`Pr(>F)`)[2]
      pbase <- (anova(base_mod,cluster_mod)$`Pr(>F)`)[2]
      ariscad <- adjustedRandIndex(dat_i$SCADDING,dat_i$cluster)
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


### Unique sample sizes and overlap in bootstraps
n_unique_boot <-
  apply(
    boot_label[,1:321],
    1,
    function(v)length(unique(v))
  )
summary(n_unique_boot)
n_unique_overlap_boot <-
  expand.grid(i = 1:nsim, j = 1:nsim) %>%
  filter(i<j) %>%
  apply(
    1,
    function(ij){
      boot_i_obs_unique <- unique(unlist(boot_label[ij[1],1:321]))
      boot_j_obs_unique <- unique(unlist(boot_label[ij[2],1:321]))
      n_unique_overlap_ij <- length(intersect(boot_i_obs_unique,boot_j_obs_unique))
      return(n_unique_overlap_ij)
    }
  )
saveRDS(n_unique_overlap_boot,paste0(path_pdrive,"SensitivityAnalyses/sensitivity2uniqueoverlap.RDS"))
n_unique_overlap_boot <- readRDS(paste0(path_pdrive,"SensitivityAnalyses/sensitivity2uniqueoverlap.RDS"))
summary(n_unique_overlap_boot)

###Fit clusters and Scadding have slightly more agreement than in previous
#  sensitivity analysis, but is still terrible (which is good).
hist(sensitivity2_summary$ariscad)

###We can see that additional overlap between Scadding and the fit clusters
#  in these p-values. It still looks pretty good, though!
pmax_base <- max(sensitivity2_summary$pbase)
pmax_scad <- max(sensitivity2_summary$pscad)
hist(-log(sensitivity2_summary$pscad,10),xlim=c(0,30))
hist(-log(sensitivity2_summary$pbase,10),xlim=c(0,30))

###proportion below 0.05
prop_005_base <- sum(sensitivity2_summary$pbase < 0.05)/nrow(sensitivity2_summary)
prop_005_scad <- sum(sensitivity2_summary$pscad < 0.05)/nrow(sensitivity2_summary)

###proportion below 0.01
prop_001_base <- sum(sensitivity2_summary$pbase < 0.01)/nrow(sensitivity2_summary)
prop_001_scad <- sum(sensitivity2_summary$pscad < 0.01)/nrow(sensitivity2_summary)


sensitivity2_pvals_N_nboot <-
  list(
    n_unique_boot = summary(n_unique_boot),
    n_unique_overlap_boot = summary(n_unique_overlap_boot),
    maxariscad = max(sensitivity2_summary$ariscad),
    pmax_base = pmax_base,
    pmax_scad = pmax_scad,
    prop_001_base = prop_001_base,
    prop_001_scad = prop_001_scad
  )

saveRDS(sensitivity2_pvals_N_nboot,paste0(path_pdrive,"SensitivityAnalyses/sensitivity2_pvalueNnbootsummary.RDS"))


