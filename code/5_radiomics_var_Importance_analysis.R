#############################################################
###Program name:5_radiomics_var_importance_analysis
###Created by: Zachary Buchalski modified from code by Sarah Ryan
###Purpose: This code creates tables for the analysis using the variables
###         with the most weight in our RSKC clustering. It requires running code files
###         0 and 1 beforehand to get the correct data file and the names of our 
###         important variables. This creates two tables, Table_5 and Table_5_alt
###         I'm pretty sure Table_5_alt is deprecated but leaving it until I'm sure.
###################################################################

############################################# Libraries
library(tidyverse) #Data manipulation
library(magrittr) #Pipe

############################################# Paths
##Location of code and output
path_repo <- "/radiomics/"
##Location of your data
path_pdrive <- "/DataLibrary/"

############################################# Data
dat_read <- readRDS(paste0(path_pdrive,"radiomics_data/dat_clean_FVCsorted_clusters_VAS.rds"))
names_important_vars <- readRDS(paste0(path_pdrive,"radiomics_data/important_vars_radiomics_paper.rds"))
dat <- dat_read %>%
  mutate(agec = AGE_YR - mean(AGE_YR),
         BMIc = BMI - mean(BMI),
         heightc = PEHGHT - mean(PEHGHT))

############################################# Radiomics Models
# Are our radiomic measures associated with lung function?
# Put all radiomics in the same model
# Adjusted rsquare -- using all five radiomic features
adj_formula <- paste(c("AGE_YR","GENDER","factor(raceethnic)","BMI","PEHGHT"),collapse=" + ")
rad_formula <- 
  paste0("outcome ~ ", paste(names_important_vars,collapse = " + "), " + ",
         adj_formula) %>%
  as.formula
rad_formula_c <- 
  paste0("outcome ~ ", paste(names_important_vars,collapse = " + "), " + ",
         adj_formula) %>%
  as.formula
base_formula <- paste0("outcome ~ ", adj_formula) %>%
  as.formula
base_formula_c <- paste0("outcome ~ ", adj_formula) %>%
  as.formula
vars_outcome <- c('BDFVC', 'BDFEV1', 'BDFEVFVC', 'BDDLCO',
                  'fas.overall', 'gerdq.overall', 'cfq.overall', 'sobq.overall', 'promis.overall', 
                  'sf12.overall.pcs', 'sf12.overall.mcs','logCD4','logcrp')

mymodel <- function(outcome){
  dat_temp <- dat
  dat_temp$outcome <- dat_temp[,outcome] %>% unlist()
  mod1 <- lm(rad_formula, data = dat_temp)
  mod1_c <- lm(rad_formula_c, data = dat_temp)
  red <- lm(base_formula, data = dat_temp)
  red_c <- lm(base_formula_c, data = dat_temp)
  sum.mod1 <- summary(mod1)
  betaci <- paste0(sprintf("%.2f", sum.mod1$coefficients[2:6,1]), 
                   " (", sprintf("%.2f", sum.mod1$coefficients[2:6,2]), ") ", 
                   signif(sum.mod1$coefficients[2:6,4],2))
  overall.pvalue <- round(anova(red,mod1)$`Pr(>F)`[2],3)
  if(overall.pvalue == "0"){overall.pvalue = "<0.0001"}
  r2 <- sum.mod1$r.squared
  r2base <- summary(red)$r.squared
  result <- c(betaci, overall.pvalue, 
              round(r2,3), round(summary(mod1_c)$adj.r.squared,3),
              round(r2base,3), round(summary(red_c)$adj.r.squared,3))
  return(result)
}
myresults <- sapply(vars_outcome,mymodel) %>% t %>%
  as.data.frame
colnames(myresults) <- c(names_important_vars, "Pvalue", "Rsq", "AdjRsq", "Rsq - base", "AdjRsq - base")
rownames(myresults) <- c('FVC','FEV1','FEV1/FVC','DLCO','FAS','GERDQ','CFQ',
                         'SOBQ','Promis','SF12 - Physical','SF12 - Mental',
                         'CD4 (log)','CRP (log)') 
saveRDS(myresults,paste0(path_repo,"results/Table_5_v2.RDS"))

############################################# Cluster Models
cl_formula <- paste0("outcome ~  factor(FVCordered_rskc) + factor(SCADDING) + ",
                     adj_formula) %>%
  as.formula()

mymodel2 <- function(outcome){
  dat_temp <- dat
  dat_temp$outcome <- dat_temp[,outcome] %>% unlist()
  mod2 <- lm(cl_formula, data = dat_temp)
  pv <- round(car::Anova(mod2)$`Pr(>F)`[1:2],4)
  if(pv[1] == "0"){pv[1] = "<0.001"}
  if(pv[2] == "0"){pv[2] = "<0.001"}
  r2 <- summary(mod2)$adj.r.squared
  result <- c(pv, round(r2,3))
  return(result)
}
myresults2 <- sapply(vars_outcome,mymodel2) %>% t %>%
  as.data.frame
rownames(myresults2) <- c('FVC','FEV1','FEV1/FVC','DLCO','FAS','GERDQ','CFQ','SOBQ','Promis','SF12 - Physical','SF12 - Mental','CD4 (log)','CRP (log)')
colnames(myresults2) <- c('RSKC', 'Scad','AdjRsq_CLSCAD')
myresults2 <- cbind(myresults2,ADJRsq_rad = myresults[rownames(myresults2),"AdjRsq"],
                    ADJRsq_base = myresults[rownames(myresults2),"AdjRsq - base"])
saveRDS(myresults2,file=paste0(path_repo,"results/Table_5_alt_v2.RDS"))

