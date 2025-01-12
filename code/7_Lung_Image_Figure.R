#############################################################
###Program name: 7_Lung_Image_Figure.R
###Created by: 
###Purpose: This code produces the images of lungs of individuals with the
### minimum, maximum, and median values of our "important" radiomics variables
### as far as informing our RSKC clusters. It is saved as "Lung_Comparison_Plot.png"
### This file requires the data files from running code files 0 and 1.
###################################################################

############################################# Pick what images to use
# true for old images in SlicerMaskPipeline/2_CT_Masked_R/
# false for new images in 2_resampled_centered/
# we probably want true
OldPipelineTF <- T

############################################# Libraries
# library(oro.nifti) #Get slice images
library(tidyverse)
library(ANTsR)
# library(cowplot) #Create panel of plots
# library(imager)
# library(png)

############################################# Paths
##Location of your code and output
path_repo <- "/radiomics/"
path_figures <- paste0(path_repo,"figures/")
##Location of your data
path_pdrive <- "/DataLibrary/"
path_rad <- paste0(path_pdrive,"radiomics_data/")
if(OldPipelineTF){
  path_img <- paste0(path_pdrive,"GRADS_CT/DataProcessed/SlicerMaskPipeline/2_CT_Masked_R/")
}else{
  path_img <- paste0(path_pdrive,"GRADS_CT/DataProcessed/")
}


############################################# Helper functions
source(paste0(path_repo,"code/7_Lung_Image_Figure_helperfunctions.R"))

############################################# Data
dat <- readRDS(paste0(path_rad,"dat_clean_FVCsorted_clusters_VAS.rds"))
importantVars <- readRDS(paste0(path_rad,"important_vars_radiomics_paper.rds"))
if(!OldPipelineTF){img_filenames <- list.files(paste0(path_img,"2_resampled_centered/"))}

############################################# Select images to visualize
#### automatically select GRADSIDS. Can be manually modified in next bit
pick_ids <-
  expand.grid(
    Order = c("min","median","max"),
    Variable = importantVars[c(2,3,5)], # choose what vars you want to plot
    slice = NA # slice selection defaults to median 
  ) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(GRADSID = f_pickid(Order,Variable)) %>%
  mutate(Variable2 = str_split(Variable,fixed("."))[[1]][
    length(str_split(Variable,fixed("."))[[1]])
  ]) %>%
  ungroup()

#### if you want to handpick a different person or slice, do so now!
# pick_ids$gradsids[] <- 
# pick_ids$slice[] <- 
# pick_ids

#### Locate images
if(OldPipelineTF){
  pick_ids <- pick_ids %>%
    mutate(pathIMG = paste0(path_img,"Sarc_",GRADSID,"_right.RDS"),
           pathMASK = NA)
}else{
  pick_ids <- pick_ids %>%
    rowwise() %>%
    mutate(
      pathIMG = paste0(path_img,
                       "2_resampled_centered/",
                       img_filenames[str_detect(toupper(img_filenames),GRADSID)])
    ) %>%
    mutate(pathMASK = str_replace(pathIMG,"2_resampled","3_masked")) %>%
    ungroup()
}

############################################# Fetch and process images
# read in images
# crop images down
# select slice
# minimally expand images back up to have the same dimension
# convert images to long format
# add pick_ids info back in
dat_img <- f_doAllTheThings(pick_ids,OldPipelineTF=OldPipelineTF)


############################################# Plot Kurtosis Histogram
g_kurtosis <- dat_img %>%
  filter(!is.na(HU)) %>%
  mutate(Order = str_to_title(Order)) %>%
  filter(Variable2 == "Kurtosis") %>% ggplot(aes(x=HU, color = Order)) + 
  geom_density(size=2) +
  scale_color_viridis_d() + 
  theme_bw() + 
  theme(legend.position="bottom")
  
height = 4.8
width = 8

pdf(paste0(path_figures,"Kurtosis_Plot.pdf"),
    height=height,width=width)
g_kurtosis
dev.off()

png(paste0(path_figures,"Kurtosis_Plot.png"),
    height=height*75,width=width*75)
g_kurtosis
dev.off()

############################################# Plot

#### Warning! aspect ratio is manually controlled during 
#### save as automatic specification in ggplot conflicts 
#### with free scales, which reduces white space. If you 
#### choose to change the number of variables for which
#### lungs are display, please check the dimension of
#### the saved image to ensure a good aspect ratio
g_exampleimgs <- dat_img %>%
  filter(!is.na(HU)) %>%
  mutate(Order = str_to_title(Order)) %>%
  ggplot(aes(x=x,y=y,fill=HU)) +
  geom_tile() +
  scale_fill_distiller(type = "seq",
                        direction = -1,
                        palette = "Greys",
                       na.value = "white") +
  facet_grid(Variable2 ~ Order,space="free",scales="free") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        plot.background=element_blank(),
        legend.position="none",
        panel.spacing = unit(0, "lines"))

height = 6
width = 4.8

pdf(paste0(path_figures,"Lung_Comparison_Plot.pdf"),
    height=height,width=width)
g_exampleimgs
dev.off()

png(paste0(path_figures,"Lung_Comparison_Plot.png"),
    height=height*75,width=width*75)
g_exampleimgs
dev.off()




