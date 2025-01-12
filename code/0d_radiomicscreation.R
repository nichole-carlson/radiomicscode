##Program name: 0d_radiomicscreation.R
##Purpose: To read in the CT imaging data and create radiomics profile
##

#################################
# Libraries
#################################
library(tidyverse)
library(magrittr)
library(RIA)
library(oro.nifti)
#################################
# Paths
#################################

##Paths of data
path_pdrive <- "/DataLibrary/"
raddir<- paste0(pagth_pdrive,"radiomics_data/")
ctdir <- paste0(path_pdrive,"GRADS_CT/")

##This is imaging data that has been resampled to the desired resolution (1x1x1 in our case)
path_img2 <- paste0(path_pdrive,"GRADS_CT/DataProcessed/2_resampled_centered/")

##This is the corresponding lung mask that you will have created in your image processing pipeline
path_img3 <- paste0(path_pdrive,"GRADS_CT/DataProcessed/3_masked_centered/")

##This is a curated list of who has data available to be processed. GRADS provides an image quality file.
##We also aligned the IDS with the available clinical data.
df <- read_csv(paste0(ctdir,'Subjects_to_Include.csv'))

##Create a data structure for making full file of radiomics data
allradiomicpaths<-matrix(nrow=2*nrow(df),ncol=3)
colnames(allradiomicpaths)<-c("gradsid","lung","savepath")

##Loop through for all patients right lung.  Could probably make this an lapply but didn't for now.
for (i in 1:nrows(df)) {
    img.file<-paste0(path_img2,str_to_lower(df$gradsid[i]),".nii.gz")
    mask.file<-paste0(path_img3,str_to_lower(df$gradsid[i]),".nii.gz")
    save.file<-paste0(path_localdata,"radiomics_right_",df$gradsid[i])
    
    ##Use later for data processing to make a clean file
    allradiomicpaths[i,1]<-df$gradsid[i]
    allradiomicpaths[i,2]<-"right"
    allradiomicpaths[i,3]<-save.file
    
    ##The mask file has 1 for right and 2 for left.  Will need to run this 2 times for each lung.
    image.dat<-load_nifti(img.file,mask_filename = mask.file, keep_mask_values = 1)
    rm(img.file)
    rm(mask.file)
    
    rad.sum<-radiomics_all(image.dat,bins_in = c(16),fo_discretized = FALSE, equal_prob = TRUE, geometry_discretized = FALSE)
    
    ##Note with a version update in RIA the variable names have a different ordering in their components.
    ##But, the interpretation compared to the paper and our dataset is easy to map.
    save_RIA(rad.sum,save_to = save.file,stats = c("stat_fo", "stat_glcm_mean"))
}

##Loop through for all patients left lung.  Could probably make this an lapply but didn't for now.
for (i in 1:nrows(df)) {
  img.file<-paste0(path_img2,str_to_lower(df$gradsid[i]),".nii.gz")
  mask.file<-paste0(path_img3,str_to_lower(df$gradsid[i]),".nii.gz")
  save.file<-paste0(path_localdata,"radiomics_left_",df$gradsid[i])
  
  ##Use later for data processing to make a clean file
  allradiomicspaths[(i+nrows(df)),1]<-df$gradsid[i]
  allradiomocspaths[(i+nrows(df)),2]<-"left"
  allradiomicspaths[(i+nrows(df)),3]<-save.file
  
  ##The mask file has 1 for right and 2 for left.  Will need to run this 2 times for each lung.
  image.dat<-load_nifti(img.file,mask_filename = mask.file, keep_mask_values = 2)
  rm(img.file)
  rm(mask.file)
  
  rad.sum<-radiomics_all(image.dat,bins_in = c(16),fo_discretized = FALSE, equal_prob = TRUE, geometry_discretized = FALSE)
  
  save_RIA(rad.sum,save_to = save.file,stats = c("stat_fo", "stat_glcm_mean"))
}

########################################### Read and aggregate results
df_results <- as.data.frame(allradiomicpaths) %>%
  mutate(savepath = paste0(savepath,"RIA_stat.csv")) %>%
   filter(file.exists(savepath)) %>% # this should only be uncommented during debugging
  mutate(analyses = lapply(savepath,read_csv)) %>%
  rowwise 

analyses <- df_results$analyses%>%
  do.call(bind_rows,.)

df_results <- cbind(df_results,analyses) %>%
              dplyr::select(gradsid,lung,savepath,Mean__orig:Skew__ep_b16_d1_mean)

# Write results to file
write.csv(df_results, paste0(raddir,"radiomics.csv"))
