#############################################################
###Program name: 3_Figure_1.R
###Created by: Zachary Buchalski modified from code by Sarah Ryan
###Purpose: This code contains the heat map in the radiomics report, saved as 
###fig1_heatmap.png. It requires only the dat_clean.csv file 
###and does not require any code files be run before it.
###################################################################

############################################# Libraries
library(tidyverse) #Data manipulation
library(gplots) #Heatmap function with legend
library(viridis) #Colorblind friendly color palette

############################################# Paths
##Location of code and output
path_repo <- "/radiomics/"
##Locataion of your data
path_pdrive <- "/DataLibrary/"

############################################# Data
dat <- read_csv(paste0(path_pdrive,"radiomics_data/dat_clean.csv"))
rad <- dat %>% select(c(mask2.first_order.Mean:mask1.glcm.ep_b16_d1_mean.Uniformity)) #Only radiomics panel variables
###Figure 1
cormat <- rad %>% as.matrix() %>% cor()

htmp <- heatmap.2(cormat, distfun = function(c) as.dist(1-c), symm = TRUE,scale="none",
        dendrogram='both',labRow = FALSE,labCol = FALSE,trace="none",col=viridis)
png(filename=paste0(path_repo,"figures/fig1_heatmap_v2.png"),width=800,height=800)
eval(htmp$call)
dev.off()
pdf(file=paste0(path_repo,"figures/fig1_heatmap_v2.pdf"),width=10,height=10)
eval(htmp$call)
dev.off()
