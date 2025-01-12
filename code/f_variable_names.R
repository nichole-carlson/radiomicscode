library(stringr)
library(tidyverse)

coding_error <- "ERRORINCODING"


namesinmodels <-
  c("GENDERFemale","GENDERMale","RACEWhite","RACEBlack","RACEHispanic",
    "RACEOther","PEHGHT","BMI","AGE_YR","Oberstein_Overall",
    "BDFVC","POSTFVC","BDFEV1","POSTFEV1","BDFEVFVC","POSTFEVFVC","BDDLCO",
    "Base",# not actually a variable. just for base models
    "SCADDING0","SCADDINGI","SCADDINGII","SCADDINGIII","SCADDINGIV",
    "Med_Lymphadenopathy1","Hilar_Lymphadenopathy1","Micronodule1",
    "Micronodule_dist1","Micronodule_dist2","Micronodule_dist3",
    "Micronodule_dist4","Micronodule_conglomerate1","BVB_Distortion1",
    "Traction_Bronchiectasis1","Ground_Glass1","Honeycombing1",
    "Reticular_Abnormality1","Consolidation1","Mosaic_Attenuation1",
    "Interlobular_Septal_Thick1","Cystic_Changes1",
    "Distribution_Cranial_Caudal1","Distribution_Cranial_Caudal2",
    "Distribution_Axial1","Distribution_Axial2","Oberstein_BVB1",
    "Oberstein_BVB2","Oberstein_BVB3","Oberstein_PC1","Oberstein_PC2",
    "Oberstein_PC3","Oberstein_ND1","Oberstein_ND2","Oberstein_ND3",
    "Oberstein_LS1","Oberstein_LS2","Oberstein_PLT1","Oberstein_LN1",
    "Oberstein_LN2","Oberstein_LN3","Emphysema_Severity1","Pulmonary_Art1",
    "ThoracicSurgery1") 
namesinmodels <-
  namesinmodels %>%
  sapply(function(v){
    if(str_sub(v,-4,-1)=="FEV1"){return(v)}
    if(!is.na(as.numeric(str_sub(v,-1)))){
      return(paste0(str_sub(v,1,-2),0))
    }else{return(v)}
  }) %>%
  c(namesinmodels) %>%
  unique
levelsinmodels <- 
  c("Female","Male","White","Black","Hispanic","Other",
    rep(NA,12),c("0","I","II","III","IV"),
    str_sub(namesinmodels[24:87],-1))
# cbind(namesinmodels,levelsinmodels)


namesindata <- 
  apply(cbind(namesinmodels,levelsinmodels),1,
        function(v){
          if(is.na(v[2])){return(v[1])}
          return(str_sub(v[1],1,-1-nchar(v[2])))
        }) %>%
  unname
# cbind(namesinmodels,namesindata)


namestoprint <- 
  namesindata %>%
  sapply(function(v){
    v2 <- str_replace_all(v,"_"," ")
    if(v %in% c(
      "Micronodule","Oberstein_Overall","BMI","BVB_Distortion",
      "Traction_Bronchiectasis","Ground_Glass","Honeycombing",
      "Reticular_Abnormality","Consolidation","Mosaic_Attenuation",
      "Cystic_Changes","Distribution_Cranial_Caudal","Distribution_Axial",
      "Oberstein_BVB","Oberstein_PC","Oberstein_ND","Oberstein_LS",
      "Oberstein_PLT","Oberstein_LN","Base"
    ))return(v2)
    if(v %in% c(
      "GENDER","RACE","SCADDING"
    ))return(str_to_title(v2))
    if(v=="PEHGHT"){return("Height in inches")}
    if(v=="AGE_YR"){return("Age in years")}
    if(v=="Med_Lymphadenopathy"){return("Mediastinal LA")}
    if(v=="Hilar_Lymphadenopathy"){return("Hilar LA")}
    if(v=="Micronodule_dist"){return("Micronodule Distribution")}
    if(v=="Micronodule_conglomerate"){return("Conglomerate Mass")}
    if(v=="Interlobular_Septal_Thick"){return("Interlobular Septal Thickening")}
    if(v=="Emphysema_Severity"){return("Emphysema")}
    if(v=="Pulmonary_Art"){return("Pulmonary Artery Enlargement")}
    if(v=="BDFVC"){return("Pre-BD FVC")}
    if(v=="POSTFVC"){return("Post-BD FVC")}
    if(v=="BDFEV1"){return("Pre-BD FEV1")}
    if(v=="POSTFEV1"){return("Post-BD FEV1")}
    if(v=="BDFEVFVC"){return("Pre-BD FEV1/FVC")}
    if(v=="POSTFEVFVC"){return("Post-BD FEV1/FVC")}
    if(v=="BDDLCO"){return("Post-BD DLCO")}
    if(v=="ThoracicSurgery"){return("Prior Thoracic Surgery")}
    return(coding_error)
  })
# sum(namestoprint==coding_error)
# cbind(namesindata,namestoprint)


levelstoprint <- 
  cbind(namesindata,levelsinmodels) %>%
  apply(1,function(v){
    if(is.na(v[2])){return(NA)}
    if(v[1] %in% 
       c("GENDER","RACE","Oberstein_BVB","Oberstein_PC","Oberstein_ND",
         "Oberstein_LS","Oberstein_PLT","Oberstein_LN","SCADDING"
       )){return(v[2])}
    if(v[1] %in%
       c(
         "Med_Lymphadenopathy","Hilar_Lymphadenopathy","Micronodule",
         "Micronodule_conglomerate",
         "BVB_Distortion","Traction_Bronchiectasis","Ground_Glass",
         "Honeycombing","Reticular_Abnormality","Consolidation",
         "Mosaic_Attenuation","Interlobular_Septal_Thick","Cystic_Changes",
         "Emphysema_Severity","Pulmonary_Art","ThoracicSurgery"
       )){
      if(v[2]=="1"){return("Yes")}else{return("No")}}
    if(v[1]=="Micronodule_dist"){
      return(c("No","Perilymphatic","Peribronchovascular",
               "Both","Random")[as.numeric(v[2])+1])}
    if(v[1]=="Distribution_Axial"){
      return(c("None","Central","Peripheral")[as.numeric(v[2])+1])}
    if(v[1]=="Distribution_Cranial_Caudal"){
      return(c("None","Upper","Lower")[as.numeric(v[2])+1])}
    return(coding_error)
  }) %>%
  unname
# sum(levelstoprint==coding_error,na.rm=T)
# cbind(namesindata,levelsinmodels,levelstoprint)

levelsreference <-
  cbind(levelsinmodels,levelstoprint,namesindata) %>%
  apply(1,function(v){
    if(is.na(v[1])){return(NA)}
    if(v[1] %in% c("Female","Male")){return("Female")}
    if(v[1] %in% c("White","Black","Hispanic","Other")){return("White")}
    if(v[2] %in% c("Yes","No")){return("No")}
    if(v[3] %in% c("Distribution_Cranial_Caudal",
                   "Distribution_Axial",
                   "Micronodule_dist")){return("No")}
    if(v[1] %in% as.character(0:10)){return("0")}
    if(v[3]=="SCADDING"){return("0")}
    return(coding_error)
  })
# sum(levelsreference==coding_error,na.rm=T)
# cbind(namesindata,levelstoprint,levelsreference)

groups <-
  namesindata %>%
  sapply(function(v){
    if(v=="Base"){return("Base")}
    if(v %in% c("GENDER","RACE","PEHGHT","BMI","AGE_YR")){return("Base Factors")}
    if(v %in% c("BDFVC","POSTFVC","BDFEV1","POSTFEV1","BDFEVFVC","POSTFEVFVC",
                "BDDLCO")){return("PFT")}
    if(v=="SCADDING"){return("Scadding")}
    if(v %in% c("Micronodule","Micronodule_dist","Micronodule_conglomerate")){
      return("Nodules")}
    if(v %in% c("BVB_Distortion","Traction_Bronchiectasis")){
      return("Airway and Vascular Distortion")}
    if(v %in% c("Med_Lymphadenopathy","Hilar_Lymphadenopathy")){
      return("Lymphadenopathy")}
    if(v %in% c("Ground_Glass","Honeycombing","Reticular_Abnormality",
                "Consolidation","Mosaic_Attenuation",
                "Interlobular_Septal_Thick","Cystic_Changes")){
      return("Parenchymal Opacity and Distortion")}
    if(v %in% c("Distribution_Cranial_Caudal","Distribution_Axial")){
      return("Distribution and Pattern of Distortion")}
    if(v %in% c("Oberstein_BVB","Oberstein_PC","Oberstein_ND","Oberstein_LS",
                "Oberstein_PLT","Oberstein_LN","Oberstein_Overall")){
      return("Oberstein Components")}
    if(v %in% c("Emphysema_Severity","Pulmonary_Art","ThoracicSurgery")){
      return("Other Findings")}
    return(coding_error)
  })
# sum(groups==coding_error)
# cbind(namesindata,groups) %>% unique
  

var_coefficients <-
  cbind(namesindata,namesinmodels,namestoprint,
        levelsinmodels,levelstoprint,levelsreference,
        groups)
rownames(var_coefficients) <- NULL
rm(namesindata,namesinmodels,namestoprint,
   levelsinmodels,levelstoprint,levelsreference,
   groups)
gplevels <- c(
  "Base","Base Factors","PFT","Lymphadenopathy","Nodules",
  "Airway and Vascular Distortion",
  "Parenchymal Opacity and Distortion",
  "Distribution and Pattern of Distortion",
  "Oberstein Components","Other Findings","Scadding"
)


f_variable_names <- function(vec,from_var,to_var,groupsfactorTF=T){
  v1 <- 
    var_coefficients[,c(to_var,from_var)] %>%
    unique
  v <- v1[,to_var]
  names(v) <- v1[,from_var]
  r <- unname(unlist(v[vec]))
  if(to_var=="groups" & groupsfactorTF){
    r <- factor(r,levels=gplevels)
  }
  return(r)
}



rm(coding_error)


