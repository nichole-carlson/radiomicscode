#############################################################
###Program name: 0c_BuildComprehensive GRADSData.R
###Created by: William Lippitt
###Purpose: This file is used to produce a single GRADS data set which should
###         contain all available demographic, GRADS VAS, and PFT data as well 
###         as Scadding and scanner model for the full cohort (N=368). Ideally, 
###         the file produced should contain all relevant information for all 
###         GRADS (non-longitudinal) projects excluding imaging data and VAS 
###         data from other readers.
###Produces: dat_GRADS368.RDS, dat_GRADS368.csv
#################################


library(tidyverse)

################ Paths
pdrive <- "/DataLibrary/"
data_processed <- paste0(pdrive,"GRADS_CTandbaseline_data/DataProcessed/")

loadpath_pft <- paste0(data_processed,"dat_pft.csv")
loadpath_pro <- paste0(data_processed,"dat_pro.csv")
loadpath_lab <- paste0(data_processed,"dat_lab.csv")
loadpath_ctGRADS <- paste0(data_processed,"dat_ct_GRADS.csv")
loadpath_baseline <- paste0(data_processed,"dat_baseline.csv")
savepath_GRADS368_RDS <- paste0(data_processed,"dat_GRADS368.RDS")
savepath_GRADS368_csv <- paste0(data_processed,"dat_GRADS368.csv")


dat_baseline<- read.csv(loadpath_baseline)
dat_ctGRADS<- read.csv(loadpath_ctGRADS) %>% rename(gradsid = GRADSID)
dat_pft <- read.csv(loadpath_pft)
dat_pro <- read.csv(loadpath_pro)
dat_lab <- read.csv(loadpath_lab)

dat <- # order of loading is important. data should have 368 unique GRADSIDs
  ##Incorrect line below
  dat_baseline %>% left_join(.,dat_ctGRADS,by='gradsid') %>%
  left_join(.,dat_pft,by='gradsid') %>%
  left_join(.,dat_pro,by='gradsid') %>%
  left_join(.,dat_lab,by='gradsid') %>%
  as_tibble()

########### Site Elevation Data for Correcting DLCO
site_elevation_vec <- 
  c(1223,144,5279,427,39,597,NA,59,2389,20)# site elevations in feet above sea level
site_Pb_vec <- # bariometric pressure in kPa
  site_elevation_vec %>%
  sapply(function(a)return(101.325*(1-a*0.3048*2.25577*10^(-5))^5.25588))
names(site_elevation_vec) <- 
  c("UPitt","JH","NJH","UCSF","UPenn","Vand","NA","Yale","UA","MUSC")# site names

dat2 <-
  dat %>%
  mutate(RACEnew =
           case_when(
             ethnic=="Yes" ~ "Hispanic",
             race=="White" ~ "White",
             race=="Black or African American" ~ "Black",
             race %in% c("American Indian or Alaska Native","Asian",
                         "Multi-racial or No Primary Race","Unknown") ~ "Other"
           )) %>%
  mutate(RACEnew = apply(cbind(RACEnew,ethnic),1,function(v){if(is.na(v[2])){NA}else{v[1]}})) %>%
  mutate(
    RACEnew = factor(RACEnew, levels=c("White","Black","Hispanic","Other")),
    gender = factor(gender,levels = c("Female","Male")),
    Site = names(site_elevation_vec)[site_id],
    Elevation = site_elevation_vec[site_id],
    Pb = site_Pb_vec[site_id]
  ) %>%
  mutate(
    BDDLCOcorAltitude = bddlco*(0.505+0.00488*Pb) # BDDLCO corrected for altitude
  ) %>%
  dplyr::select(-bddla) %>% # BDDLA has very poor quality and should not be used
  as_tibble
  

saveRDS(dat2,savepath_GRADS368_RDS)
write.csv(dat2,savepath_GRADS368_csv,row.names = F)



