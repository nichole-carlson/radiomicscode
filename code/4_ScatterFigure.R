#############################################################
###Program name: 4_ScatterFigure.R
###Created by: William Lippitt
###Purpose: This code produces a two-panel plot with scatter plots
###         of the two most important radiomics features colored by
###         either overall Oberstin score or by Scadding stage.
###################################################################

############################################# Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)

############################################# Paths
##Location of your code and output
path_figures <- "/radiomics/figures/"
##Location of your data
path_pdrive <- "/DataLibrary/radiomics_data/"

############################################# Data
dat_raw <- readRDS(paste0(path_pdrive,"dat_clean_FVCsorted_clusters_VAS.rds"))
important_vars <- readRDS(paste0(path_pdrive,"important_vars_radiomics_paper.rds"))
dat <- dat_raw %>%
  dplyr::select(Oberstein_Overall,SCADDING,SCMODEL,FVCordered_rskc,all_of(important_vars[1:5]))
names(dat) <- names(dat) %>% 
  str_replace("mask1.glcm.ep_b16_d1_mean.","GLCM_") %>% 
  str_replace("mask2.glcm.ep_b16_d1_mean.","GLCM_") %>%
  str_replace("mask2.first_order.","")
names(dat)

g_ober <- dat %>%
  filter(!is.na(Oberstein_Overall)) %>%
  ggplot(aes(x=GLCM_Min,y=Kurtosis,col=Oberstein_Overall)) +
  geom_point() +
  # scale_color_viridis_b(breaks = 1:5*2) +
  # scale_color_viridis_b(breaks = c(4,8)) +
  scale_color_viridis_b(name = "Oberstein Overall") +
  theme_bw() + 
  theme(legend.position="bottom")

g_scad <- dat %>%
  filter(!is.na(SCADDING)) %>%
  # mutate(SCADDING = factor(SCADDING)) %>%
  ggplot(aes(x=GLCM_Min,y=Kurtosis,col=SCADDING)) +
  geom_point() +
  # scale_color_viridis_d(option = "H",begin = 0.1,end = 0.9) + 
  scale_color_viridis_d(name = "Scadding") + 
  theme_bw() + 
  theme(legend.position="bottom")

g_combo <- ggarrange(g_ober,g_scad,ncol=2,align="h")
g_combo

height = 5
width = 8

pdf(paste0(path_figures,"RadiomicScatterPlot_v2.pdf"),
    height=height,width=width)
g_combo
dev.off()

png(paste0(path_figures,"RadiomicScatterPlot_v2.png"),
    height=height*75,width=width*75)
g_combo
dev.off()


##Other plots
dat %>%
  mutate(FVCordered_rskc = factor(FVCordered_rskc)) %>%
  ggplot(aes(x=GLCM_Sum_energy,y=GLCM_Min,col=FVCordered_rskc)) +
  geom_point() +
  scale_color_viridis_d() + 
  theme(legend.position="bottom")
dat %>%
  mutate(FVCordered_rskc = factor(FVCordered_rskc)) %>%
  ggplot(aes(x=GLCM_Sum_energy,y=Kurtosis,col=FVCordered_rskc)) +
  geom_point() +
  scale_color_viridis_d() + 
  theme(legend.position="bottom")
dat %>%
  mutate(FVCordered_rskc = factor(FVCordered_rskc)) %>%
  ggplot(aes(x=GLCM_Min,y=Kurtosis,col=FVCordered_rskc)) +
  geom_point() +
  scale_color_viridis_d() + 
  theme(legend.position="bottom")

dat %>%
  ggplot(aes(x=GLCM_Sum_energy,y=GLCM_Min,col=SCMODEL)) +
  geom_point() +
  scale_color_viridis_d(option = "H") + 
  theme(legend.position="bottom")
