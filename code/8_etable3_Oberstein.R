#############################################################
###Program name: 8_etable3_Oberstein.R
###Created by: William Lippitt
###Purpose: Make a table summarizing Oberstein scores stratified
###         by radiomics grouping.
###################################################################
####################################### libraries
library(tidyverse) #Data manipulation
library(arsenal) # for tableby
library(kableExtra)

####################################### paths
##Location of your data
path_pdrive <- "/radiomics_data/"
##Location of code and output
path_results <- "/radiomics/results/"

####################################### read in data
df_rad <- readRDS(paste0(path_pdrive,"dat_clean_FVCsorted_clusters_VAS.rds")) %>%
  dplyr::select("FVCordered_rskc","Oberstein_BVB","Oberstein_PC",
                "Oberstein_ND","Oberstein_LS","Oberstein_PLT","Oberstein_LN") %>%
  mutate_all(factor) %>%
  filter(!is.na(Oberstein_BVB))

####################################### Make table
set.seed(1)
t <- tableby(FVCordered_rskc ~ .,
        data = df_rad,
        control=tableby.control(digits=2,
                                digits.p=2,
                                cat.simplify=T,
                                cat.test="fe",
                                simulate.p.value = T,
                                B=20000)) %>%
  summary(title = "Distribution of each Oberstein component measure by radiomics group") %>%
  kable() %>%
  kable_styling(bootstrap_options = "condensed", font_size = 11) %>%
  add_header_above(c(" "=1,"Radiomics Group"=4," "=2))
saveRDS(t,paste0(path_results,"etable_Oberstein_v2.RDS"))

