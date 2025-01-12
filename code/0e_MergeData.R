#########################
##Program Name: 0d_MergeData.R
##
##Purpose: To integrate the radiomics imaging data together with the GRADS clinical data.
##        1) You will need to obtain the HRCT data from GRADS to make the radiomics data. 
##           File 0_CreateRadiomics.R reads in masked Left and Right CT data.
##        2) GRADS clinical data need to be obtained directly from GRADS. We used data from XX 
##           to obtain the demographic data, PFT, VAS, and CT parameters data sources.
##        3) When interleaving these data N=321 will be the final sample size for formal analysis.
##########################

#################################
# Libraries
#################################
library(tidyverse)
library(magrittr)
library(ez.combat)

#################################
# Paths
#################################

##Location of your data
path_pdrive <- "/DataLibrary/"
datadir <- paste0(path_pdrive,"GRADS_CTandBaseline_data/DataProcessed/")
raddir <- paste0(path_pdrive,"radiomics_data/")
ctdir <- paste0(path_pdrive,"GRADS_CT/")

#################################
# Clean Data
#################################

###Read in the clinical data
gradsdat<-readRDS(paste0(datadir,"dat_GRADS368.RDS"))


# Read in subjects that we want to include from the clinical dataset
# Select out the individuals with viable CT data. This file is available with the paper code.
# There is code that goes with this.  0_Find_IDs_To_Include.R in GRADS CT pipeline code.
df <- read_csv(paste0(ctdir, 'Subjects_to_Include.csv'))
ids <- df$gradsid %>% toupper()
gradsdat %<>% filter(GRADSID %in% ids) # 321 obs


# Also read in radiomic data, note: GLCM means are all the same, so have to remove
dc <- read_csv(paste0(raddir,'dat_raw_radiomics.csv'))%>%
  select(id, mask2.first_order.Mean:mask2.glcm.ep_b16_d1_mean.Entropy, 
         mask2.glcm.ep_b16_d1_mean.Median:mask1.glcm.ep_b16_d1_mean.Entropy, 
         mask1.glcm.ep_b16_d1_mean.Median:mask1.glcm.ep_b16_d1_mean.Uniformity) %>%
  rename(GRADSID = id)

# Combine the datasets!
gradsdat <- inner_join(gradsdat, dc, by = 'GRADSID')
rm(dc)

# GRADS368 scanner data have been updated using header data, and so are no longer
# read in from elsewhere

###################################
# Harmonize the data 
###################################
###Needed to fix race variable: data_clean_add_clusters has white and hispanic together
gradsdat <-gradsdat %>% 
  mutate(raceethnic = case_when(RACE=="White" & (ETHNIC == "No" | is.na(ETHNIC)) ~ 1,
                                ETHNIC == "Yes" ~ 2,
                                ETHNIC == "No" & RACE=="Black or African American" ~ 3,
                                (ETHNIC == "No" & RACE %in% c("American Indian or Alaska Native","Asian", "Multi-racial or No Primary Race", "Unknown")) | is.na(RACE) ~ 4))

##make race variable white vs. other
gradsdat <- gradsdat %>% mutate(racewhite = case_when(raceethnic == 1 ~ 1,
                                                 raceethnic != 1 ~ 0))

##Remove missing people if missing variables needing to batch on.  This is really one individual wiht no gender information
gradsNA<- gradsdat%>%filter(!is.na(GENDER))

# Select only the radiomic data and scale
raddf <- gradsNA %>% dplyr::select(mask2.first_order.Mean:mask1.glcm.ep_b16_d1_mean.Uniformity) %>% scale()
dim(raddf) 

#raddf<-cbind(gradsNA$GRADSID,raddf)

mod <- as.data.frame(model.matrix(~ -1 + gradsNA$AGE_YR + gradsNA$PEHGHT + gradsNA$BMI + factor(gradsNA$GENDER) + factor(gradsNA$raceethnic) + factor(gradsNA$PHENGRP)))
mod <- mod %>%
          dplyr::select(-`factor(gradsNA$GENDER)Male`)
names(mod) <- paste0("ez",1:ncol(mod))
fmla <- paste0("~ ",paste(names(mod),collapse = "+"))
df_toharm <- cbind(raddf,mod,SCMODEL = gradsNA$SCMODEL)
df_harm <- 
  ez.combat(df = df_toharm, adjust.var = colnames(raddf),
            batch.var = "SCMODEL", output="overwrite", 
            verbose = TRUE, model = fmla)$df %>%
  as_tibble()

df_harm <- df_harm %>%
  dplyr::select(all_of(colnames(raddf))) %>%
  scale
check_beforeharm <- colnames(raddf) %>% # association with SCMODEL and BMI
  sapply(function(nm){
    dummy <- gradsNA
    dummy$rad <- unlist(dummy[,nm])
    c(BMI = summary(lm(rad ~ BMI,data=dummy))$r.squared,
      SCMODEL = summary(lm(rad ~ factor(SCMODEL),data=dummy))$r.squared)
  })
dat <- 
  bind_cols(
    dplyr::select(gradsNA,-all_of(colnames(raddf))),
    df_harm
  ) %>%
  as_tibble
check_afterharm <- colnames(raddf) %>%  # association with SCMODEL and BMI
  sapply(function(nm){
    dummy <- dat
    dummy$rad <- unlist(dummy[,nm])
    c(BMI = summary(lm(rad ~ BMI,data=dummy))$r.squared,
      SCMODEL = summary(lm(rad ~ factor(SCMODEL),data=dummy))$r.squared)
  })
# Before and after r squared values
cbind(
  apply(check_beforeharm,1,summary),
  apply(check_afterharm,1,summary)
) %>%
  round(5) 
# looks good! BMI about the same, SCMODEL much lower (nonzero bc emp bayes)
# 
# check weird duplicate values issue from earlier:
colnames(raddf) %>%
  sapply(function(nm){
    dat[,nm] %>% unlist %>% table %>% table
  }) %>%
  t %>% unname %>% c %>% unique # if just nrow(dat), we're good! all values unique


###################################
# Save the data 
###################################

write_csv(dat, paste0(raddir,"dat_clean.csv"))


