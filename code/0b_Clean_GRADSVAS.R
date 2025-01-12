#############################################################
###Program name: 0b_Clean_GRADSVAS.R
###Created by: William Lippitt
###Purpose: This file is used to process VAS reads of the GRADS
###         cohort by the GRADS radiologist. Data quality is not assessed.
###Produces: dat_ct_GRADS.csv
#################################


library(tidyverse)
library(readxl)
library(stringr)

################ Paths
##Location of your data
pdrive <- "/DataLibrary/"
data_raw <- paste0(pdrive,"GRADS_CTandbaseline_data/DataRaw/")
data_processed <- paste0(pdrive,"GRADS_CTandbaseline_data/DataProcessed/")

loadpath_vas <- paste0(data_raw,"GRADSCTReads_091317.csv")
savepath_vas <- paste0(data_processed,"dat_ct_GRADS.csv")

################ Remove repeat and make Oberstein Overall
vas <- 
  read.csv(loadpath_vas) %>% 
  unique() %>%
  as_tibble() %>%
  mutate(
    Oberstein_Overall = apply(cbind(Oberstein_BVB,Oberstein_PC,Oberstein_ND,
                                    Oberstein_LS,Oberstein_PLT,Oberstein_LN),
                              1,sum)
  ) %>%
  dplyr::select(-Visit,-DOB)

write.csv(vas,savepath_vas,row.names=FALSE)
