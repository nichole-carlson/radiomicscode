#############################################################
###Program name: 1_Tables_1_and_2.R
###Created by: Originally Sarah Ryan. Edited by Zac Buchalski, 
###             Nichole E Carlson, and William Lippitt
###Purpose: This code makes the table 1-3 in the radiomics paper. It 
###         requires running code file 0f_filter_and_Cluster.R to create
###         cluster labels, which were appended to the original dataset.
###         
###         Output: Three RDS files are created, two which correspond to alternate
###         versions of Table 1 that Nichole can choose between depending on how
###         she wants to handle race and one corresponding to table 2, a table
###         with information on our outcomes stratified by the clustering group.
###################################################################

############################################# Libraries
library(tableone) #Used to create tables
library(xtable) #Used to export some tables
library(tidyverse) #Data manipulation
library(arsenal) # for t3
library(kableExtra) # for t3
library(gtsummary)

############################################# Paths
#Location of code and output
path_repo <- "/radiomics/"
#Location of your data
path_pdrive <- "/DataLibrary/"

############################################# Helper functions
suppressWarnings(source(paste0(path_repo,"code/f_variable_names.R"))) # for t3. in code file of repo

############################################# Data Processing

####Read in cluster/radiomic data, pre-processing
dat <- read_rds(paste0(path_pdrive,"radiomics_data/dat_clean_add_clusters.rds")) %>%
  as_tibble()

###  Note Gender Female =0
###       Race 1=white, 2=black, 3 other
###       Ethnic, 0 = non-hispanic, 1 = hispanic
###       RaceEthnic, 1 = non-hispanic white, 2 = hispanic, 3 = black, 4 = other (Other = NA, Asian, American Indian,Alaska Native, or not
##identifying a single primary race, NA)

k <- length(table(dat$rskc_cluster))-1
FVC.ordering <- dat %>%
  dplyr::select(BDFVC,rskc_cluster) %>%
  group_by(rskc_cluster) %>% summarize(meanFVC = mean(BDFVC,na.rm=T)) %>%
  arrange(meanFVC) %>% mutate(new_cluster = k:0) %>%
  arrange(rskc_cluster) %>%
  dplyr::select(new_cluster) %>% unlist %>% unname
dat$FVCordered_rskc <- FVC.ordering[dat$rskc_cluster]
dat %>%
  dplyr::select(BDFVC,FVCordered_rskc) %>%
  group_by(FVCordered_rskc) %>%
  summarize(meanFVC = mean(BDFVC,na.rm=T))


##Create the obstructive lung PFT that Lisa wanted
dat$BDFEVFVC7 <- dat$BDFEVFVC <= 0.7

###Read in visual scoring data. Note this dataset is from some cleaning done in the ctvscxr repo in the 0_datClean.R code file###
vas.dat <- read.csv(paste0(path_pdrive,"GRADS_VAS/dat_CTvCXR_descriptive.csv"))
vas.dat <- vas.dat %>%  dplyr::select(
    "GRADSID","Med_Lymphadenopathy","Hilar_Lymphadenopathy","Micronodule","Micronodule_dist",
    "Micronodule_conglomerate","BVB_Distortion","Traction_Bronchiectasis",
    "Ground_Glass","Honeycombing","Reticular_Abnormality","Consolidation",
    "Mosaic_Attenuation","Interlobular_Septal_Thick",
    "Distribution_Cranial_Caudal","Distribution_Axial","Oberstein_BVB",
    "Oberstein_PC","Oberstein_ND","Oberstein_LS","Oberstein_PLT","Oberstein_LN",
    "Oberstein_Overall"
  )

###Join the VAS with the radiomics data. N=320 if missing SCADDING stage patients in the dataset
dat <- dat %>%
  dplyr::select(GRADSID,all_of(setdiff(names(.),names(vas.dat)))) %>%
  left_join(vas.dat,by=c("GRADSID"))

###Create some y/n variables
dat <- dat %>% mutate(
        med_lymph_yn =
          case_when(
            Med_Lymphadenopathy==0 ~ "No",
            Med_Lymphadenopathy>0 ~ "Yes"
          ),
        hilar_lymph_yn =
          case_when(
            Hilar_Lymphadenopathy==0 ~ "No",
            Hilar_Lymphadenopathy>0 ~ "Yes"
          ),
        septal_thick_yn =
          case_when(
            Interlobular_Septal_Thick==0 ~ "No",
            Interlobular_Septal_Thick>0 ~ "Yes"
          ),
        cranial_caudal_yn =
          case_when(
            Distribution_Cranial_Caudal==0 ~ "No",
            Distribution_Cranial_Caudal>0 ~ "Yes"
          ),
        axial_yn =
          case_when(
            Distribution_Axial==0 ~ "No",
            Distribution_Axial>0 ~ "Yes"
          )
      ) %>%
  mutate(
    med_lymph_yn = factor(med_lymph_yn, levels=c("No","Yes")),
    hilar_lymph_yn = factor(hilar_lymph_yn, levels=c("No","Yes")),
    septal_thick_yn = factor(septal_thick_yn, levels=c("No", "Yes")),
    cranial_caudal_yn = factor(cranial_caudal_yn, levels=c("No", "Yes")),
    axial_yn = factor(axial_yn, levels=c("No", "Yes"))
  )

write_rds(dat,file=paste0(path_pdrive,"radiomics_data/dat_clean_FVCsorted_clusters_VAS.rds"))

############################################# Tables 1
#### New Table 1 (with missing)
set.seed(1)
old_rownames <- c('Sample Size', 'Age (yr) ','Male, N(%)','White, N(%)','Race/Ethnicity, N(%)','Black','Hispanic','Other','Non-Hispanic White',
                  'Height (in)','BMI (kg/m^2)','Scadding','0','1','2','3','4','N/A',
                  'FVC (L)','Scanner Model','Discovery','iCT128','LightSpeed',
                  'Optima','Sensation','SOM Def','SOM Def AS','SOM Def Flash')
wt1 <- dat %>%
  dplyr::select(FVCordered_rskc,AGE_YR,GENDER,racewhite,raceethnic,PEHGHT,BMI, SCADDING,BDFVC, 
                PHENGRP, YrsSinceDX, SCMODEL) %>%
  tbl_summary(by = FVCordered_rskc,
              label = list(
                AGE_YR ~ "Age (yr)",
                GENDER ~ "Sex",
                racewhite ~ "White",
                raceethnic ~ "Race/Ethnicity",
                PEHGHT ~ "Height (in)",
                BMI ~ "BMI (kg/m^2)",
                SCADDING ~ "Scadding",
                BDFVC ~ "FVC (L)",
                PHENGRP ~ "Enrollment Phenotype",
                YrsSinceDX ~ "Years Since DX",
                SCMODEL ~ "Scanner Model"
              ),
              missing = "no",
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 2
  ) %>%
  add_n(statistic = "{N_miss}",col_label = "**Nmiss**") %>%
  add_overall(last=T) %>%
  add_p(test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "fisher.test"),
        test.args = all_categorical() ~ list(simulate.p.value = T, B=20000))
saveRDS(wt1,paste0(path_repo,"figures/Table1_withMissing.RDS"))

############################################# Tables 2
##Create table of outcome data overall and by sorted cluster labels
SE <- function(x) sd(x,na.rm=T)/sqrt(sum(!is.na(x)))
set.seed(1)
wt2 <- dat %>%
  dplyr::select(FVCordered_rskc,BDFEV1,BDFVC,BDFEVFVC,BDFEVFVC7,BDDLCO,
                fas.overall,gerdq.overall,cfq.overall,sobq.overall,
                promis.overall,sf12.overall.pcs,sf12.overall.mcs) %>%
  tbl_summary(by = FVCordered_rskc,
              label = list(
                BDFEV1 ~ "FEV1 (L)",
                BDFVC ~ "FVC (L)",
                BDFEVFVC ~ "FEV1/FVC (%)",
                BDFEVFVC7 ~ "Obstructive Type^a N (%)",
                BDDLCO ~ "DLCO (mL/min/mmHg)",
                fas.overall ~ "FAS",
                gerdq.overall ~ "GERDQ",
                cfq.overall ~ "CFQ",
                sobq.overall ~ "SOBQ",
                promis.overall ~ "PROMIS Fatigue",
                sf12.overall.pcs ~ "SF12-Physical",
                sf12.overall.mcs ~ "SF12-Mental"
              ),
              missing = "no",
              statistic = list(
                all_continuous() ~ "{mean} ({SE})", # SE put here
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 2
  ) %>%
  add_n(statistic = "{N_miss}",col_label = "**Nmiss**") %>%
  add_overall(last=T) %>%
  add_p(test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "fisher.test"),
        test.args = all_categorical() ~ list(simulate.p.value = T, B=20000))
saveRDS(wt2,paste0(path_repo,"figures/Table2_revision_withMissing_v2.RDS"))

############################################# Tables 3
##Create table of outcome data overall and by sorted cluster labels

dichotomizer <- function(v){
  case_when(
    v>0 ~ "Yes",
    v=="0" ~ "No"
  ) %>%
    return()
}


vars_t3 <- c("Med_Lymphadenopathy","Hilar_Lymphadenopathy","Micronodule",
             "Micronodule_dist","BVB_Distortion",
             "Traction_Bronchiectasis","Ground_Glass","Honeycombing",
             "Reticular_Abnormality","Mosaic_Attenuation",
             "Interlobular_Septal_Thick","Distribution_Cranial_Caudal",
             "Distribution_Axial","Oberstein_BVB","Oberstein_PC","Oberstein_ND",
             "Oberstein_LS","Oberstein_PLT","Oberstein_LN","Oberstein_Overall")
vars_t3_red <- c("Med_Lymphadenopathy","Hilar_Lymphadenopathy","Micronodule",
             "Micronodule_dist","BVB_Distortion",
             "Traction_Bronchiectasis","Ground_Glass","Honeycombing",
             "Reticular_Abnormality","Mosaic_Attenuation",
             "Interlobular_Septal_Thick","Distribution_Cranial_Caudal",
             "Distribution_Axial","Oberstein_Overall")
dat_t3 <-
  dat %>%
  dplyr::select(all_of(vars_t3),FVCordered_rskc) %>%
  mutate(
    Hilar_Lymphadenopathy = dichotomizer(Hilar_Lymphadenopathy),
    Med_Lymphadenopathy = dichotomizer(Med_Lymphadenopathy),
    Interlobular_Septal_Thick = dichotomizer(Interlobular_Septal_Thick),
    Micronodule = dichotomizer(Micronodule),
    BVB_Distortion = dichotomizer(BVB_Distortion),
    Traction_Bronchiectasis = dichotomizer(Traction_Bronchiectasis),
    Ground_Glass = dichotomizer(Ground_Glass),
    Honeycombing = dichotomizer(Honeycombing),
    Reticular_Abnormality = dichotomizer(Reticular_Abnormality),
    Mosaic_Attenuation = dichotomizer(Mosaic_Attenuation),
    Distribution_Axial =
      f_variable_names(paste0("Distribution_Axial",Distribution_Axial),
                       from_var="namesinmodels",to_var="levelstoprint"),
    Distribution_Cranial_Caudal =
      f_variable_names(paste0("Distribution_Cranial_Caudal",Distribution_Cranial_Caudal),
                       from_var="namesinmodels",to_var="levelstoprint"),
    Micronodule_dist =
      f_variable_names(paste0("Micronodule_dist",Micronodule_dist),
                       from_var="namesinmodels",to_var="levelstoprint")
  ) %>%
  mutate(
    Distribution_Cranial_Caudal = factor(Distribution_Cranial_Caudal,
                                         levels=c("None","Upper","Lower")),
    Distribution_Axial = factor(Distribution_Axial,
                                levels=c("None","Central","Peripheral"))
  ) %>%
  mutate(
    across(all_of(setdiff(vars_t3,"Oberstein_Overall")),factor)
  ) %>% filter(complete.cases(.))

dat_t3_red <-
  dat_t3 %>%
  dplyr::select(all_of(vars_t3_red),FVCordered_rskc)

set.seed(1)
vars_t3_red2 <- f_variable_names(vars_t3_red,"namesindata","namestoprint") %>%
  as.list()
names(vars_t3_red2) <- vars_t3_red
newt3 <- dat_t3_red %>%
  tbl_summary(by = FVCordered_rskc,
              label = vars_t3_red2,
              missing = "no",
              statistic = list(
                all_continuous() ~ "{mean} ({SE})",
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 2
  ) %>%
  # add_n(statistic = "{N_miss}",col_label = "**Nmiss**") %>%
  add_overall(last=T) %>%
  add_p(test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "fisher.test"),
        test.args = all_categorical() ~ list(simulate.p.value = T, B=20000)) %>%
  as_kable() %>%
  pack_rows("Lymphadenopathy",1,2) %>%
  pack_rows("Micronodules",3,9) %>%
  pack_rows("Airway and Vascular Distortion",10,11) %>%
  pack_rows("Parenchymal Opacity and Distortion",12,16) %>%
  pack_rows("Distribution and Pattern of Distortion",17,24) %>%
  kable_styling(bootstrap_options = "condensed", font_size = 11) %>%
  add_header_above(c(" "=1,"Radiomics Group"=k+1," "=2)) %>%
  add_indent(c(5:9,18:20,22:24), level_of_indent = 2, all_cols = FALSE) %>%
  column_spec(1,bold=T)

set.seed(1)
t3 <- tableby(FVCordered_rskc ~ ., 
              data=dat_t3,
              control=tableby.control(digits=2,
                                      digits.p=2,
                                      numeric.stats = c("Nmiss","meansd"),
                                      cat.simplify=T,
                                      cat.test="fe",
                                      simulate.p.value = T))
set.seed(1)
t3_red <- tableby(FVCordered_rskc ~ ., 
              data=dat_t3_red,
              control=tableby.control(digits=2,
                                      digits.p=2,
                                      numeric.stats = c("Nmiss","meansd"),
                                      cat.simplify=T,
                                      cat.test="fe",
                                      simulate.p.value = T))
t3_labs <- f_variable_names(vars_t3,from_var="namesindata",to_var="namestoprint")
t3_labs_red <- f_variable_names(vars_t3_red,from_var="namesindata",to_var="namestoprint")
names(t3_labs) <- vars_t3
names(t3_labs_red) <- vars_t3_red
arsenal::labels(t3) <- t3_labs
arsenal::labels(t3_red) <- t3_labs_red
s3 <-
  summary(
    t3,
    title = "Table 6: Distribution of each CT VAS measure by Radiomics Group (N;%)."
  )
s3_red <-
  summary(
    t3_red,
    title = "Table 6: Distribution of each CT VAS measure by Radiomics Group (N;%)."
  )
s3 <-
  s3 %>%
  kable(
    escape=FALSE
    # caption = "Table 3: Distribution of each CT VAS measure by Radiomics Group (N;%)."
  ) %>%
  pack_rows("Lymphadenopathy",1,2) %>%
  pack_rows("Micronodules",3,9) %>%
  pack_rows("Airway and Vascular Distortion",10,11) %>%
  pack_rows("Parenchymal Opacity and Distortion",12,16) %>%
  pack_rows("Distribution and Pattern of Distortion",17,24) %>%
  pack_rows("Oberstein Components",25,52) %>%
  kable_styling(bootstrap_options = "condensed", font_size = 11) %>%
  add_header_above(c(" "=1,"Radiomics Group"=k+1," "=2))

s3_red <-
  s3_red %>%
  kable(
    escape=FALSE
    # caption = "Table 3: Distribution of each CT VAS measure by Radiomics Group (N;%)."
  ) %>%
  pack_rows("Lymphadenopathy",1,2) %>%
  pack_rows("Micronodules",3,9) %>%
  pack_rows("Airway and Vascular Distortion",10,11) %>%
  pack_rows("Parenchymal Opacity and Distortion",12,16) %>%
  pack_rows("Distribution and Pattern of Distortion",17,24) %>%
  kable_styling(bootstrap_options = "condensed", font_size = 11) %>%
  add_header_above(c(" "=1,"Radiomics Group"=k+1," "=2))

saveRDS(s3,paste0(path_repo,"figures/Table3_VAS_withOberComponents_v2.RDS"))
saveRDS(s3_red,paste0(path_repo,"figures/Table3_VAS_v2.RDS"))
saveRDS(newt3,paste0(path_repo,"figures/Table3_VAS_OverOberSE_v2.RDS"))
