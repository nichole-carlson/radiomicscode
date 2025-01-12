#############################################################
###Program name: 0a_Clean_GRADSclinical.R
###Created by: William Lippitt
###Purpose: This file is used to produce: patient and visit data to contain
###         only variables with data; pft data from enrollment visit; baseline
###         demographics as well as scadding and scanner. Note that data quality
###         is not considered here, though perhaps it should be.
###Produces: dat_visit.csv, dat_patient.csv, dat_pft.csv, dat_baseline.csv
#################################

library(tidyverse)

################ Paths
##location of all the data
pdrive <- "/DataLibrary/"
data_raw <- paste0(pdrive,"GRADS_CTandbaseline_data/DataRaw/")
data_processed <- paste0(pdrive,"GRADS_CTandbaseline_data/DataProcessed/")

loadpath_visit <- paste0(data_raw,"VisitData.csv")
loadpath_patient <- paste0(data_raw,"PatientData.csv")
savepath_visit <- paste0(data_processed,"dat_visit.csv")
savepath_patient <- paste0(data_processed,"dat_patient.csv")
savepath_pft <- paste0(data_processed,"dat_pft.csv")
savepath_lab <- paste0(data_processed,"dat_lab.csv")
savepath_pro <- paste0(data_processed,"dat_pro.csv")
savepath_baseline <- paste0(data_processed,"dat_baseline.csv")

################ Load Data
dat_visit <- read.csv(loadpath_visit)
dat_patient <- read.csv(loadpath_patient)

################ Functions
dichotomizer <- function(x,cutoff_low,cutoff_hi=Inf){
  if(is.na(x)){return(NA)}
  else{return(1*((x>=cutoff_low) & (x<=cutoff_hi)))}}
f_empty <- function(x){is.na(x) | x %in% c(""," ","N/A")} # TF is NA?
f_na <- function(x){ # make missing data NA
  if(f_empty(x)){return(NA)}
  else if(is.character(x)){return(x)}
  else if(as.numeric(x)==-2 | as.numeric(x)==-1 | as.numeric(x)==-3){return(NA)}
  else{return(as.numeric(x))}}
f_na_vec <- function(v){unname(sapply(unname(c(as.matrix(v))),f_na))} # vectorize f_na for df
form_questions <- function(dat,form_name){ # extract variable names for a given form
  colnames(dat)[which(form_name==sub("\\..*","",colnames(dat)))]}
f_visit <- function(dat){ # crosstabs for having data on first and 6m follow up visits
  mydat <- dat %>% as.data.frame() %>% select("gradsid","visit")
  dummy <- spread(mydat,visit,visit)[,c("gradsid","Enrollment / Clinic Visit","FU6M")]
  mytab <- table(dummy[,2],dummy[,3],useNA="ifany")
  return(mytab)}
f_temp <- function(temp,type){ # put all temps in fahrenheit
  if(type=="Celsius"){return(temp*1.8+32)}
  else{return(temp)}}

################ Nonempty Visit and Patient
vars_notempty <- 
  colnames(dat_visit)[which(nrow(dat_visit)!=
                              apply(dat_visit,2,function(x) sum(f_empty(x))))]
dat_visit <- dat_visit[,vars_notempty]
write.csv(dat_visit,savepath_visit,row.names=F)
rm(vars_notempty)
f_visit(dat_visit) # how many have first and/or second visit data

vars_notempty <- 
  colnames(dat_patient)[which(nrow(dat_patient)!=
                                apply(dat_patient,2,function(x) sum(f_empty(x))))]
dat_patient <- dat_patient[,vars_notempty]
write.csv(dat_patient,savepath_patient,row.names=F)
rm(vars_notempty)

################ PFT
dat_pft <- 
  dat_visit[,c("gradsid","visit",form_questions(dat_visit,"pft"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
f_visit(dat_pft) # how many have pft from first and/or second visit data
dat_pft <- dat_pft %>% 
  filter(visit=="Enrollment / Clinic Visit") %>%
  select("gradsid","pft.bdfvc","pft.postfvc","pft.fvcpred",
         "pft.bdfev1","pft.postfev1","pft.fev1pred","pft.bdsvc","pft.bdtlc",
         "pft.bdrv","pft.bdfrcp","pft.bdtgv","pft.bderv","pft.hemog","pft.bddlco",
         "pft.preddlco","pft.bddla","pft.bddva","pft.predva","pft.bddlv") %>%
  mutate(across(all_of(colnames(.)),f_na_vec))
colnames(dat_pft) <- colnames(dat_pft) %>%
  str_replace("pft\\.","") 
  
dat_pft <- dat_pft %>%
  mutate(
    bdfevfvc = bdfev1/bdfvc, 
    postfevfvc = postfev1/postfvc,
    fvcpred_01 = fvcpred %>% sapply(dichotomizer,cutoff_low = 80),
    fev1pred_01 = fev1pred %>% sapply(dichotomizer,cutoff_low = 80),
    preddlco_01 = preddlco %>% sapply(dichotomizer,cutoff_low = 80),
    bdfevfvc_01 = bdfevfvc %>% sapply(dichotomizer,cutoff_low = .7,cutoff_hi=.8),
    postfevfvc_01 = postfevfvc %>% sapply(dichotomizer,cutoff_low = .7,cutoff_hi=.8)
  )
write.csv(dat_pft,savepath_pft,row.names=F)

################ CD4&CRP lab variables####
dat_lab <- 
  dat_visit[,c("gradsid","visit",form_questions(dat_visit,"smplc"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
f_visit(dat_lab) # how many have pft from first and/or second visit data
dat_lab <- dat_lab %>% 
  filter(visit=="Enrollment / Clinic Visit") %>%
  select("gradsid","smplc.cd4","smplc.crp") %>%
  mutate(across(all_of(colnames(.)),f_na_vec))
colnames(dat_lab) <- colnames(dat_lab) %>%
  str_replace("smplc\\.","")

dat_lab <- dat_lab%>%
  mutate(logCD4 = log(cd4),
         logcrp = log(crp))
write.csv(dat_lab,savepath_lab,row.names=F)

################ PRO
dat_fas <- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"fas"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
dat_cfq <- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"cfq"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
dat_sobq<- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"sobq"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
dat_promis <- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"pf"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
dat_gerd <- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"gerdq"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
dat_sf <- dat_visit[,c("gradsid","visit",form_questions(dat_visit,"sf12"))] %>%
  as_tibble() %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )

dat_pro <- dat_fas %>% full_join(.,dat_cfq,by = c("gradsid","visit")) %>%
          full_join(.,dat_sobq, by = c("gradsid","visit")) %>%
          full_join(.,dat_gerd, by = c("gradsid","visit")) %>%
          full_join(.,dat_promis, by = c("gradsid","visit")) %>%
          full_join(.,dat_sf, by = c("gradsid","visit"))

  
###Variable names for converting scales properly including reverse coding
vars_fas <- c("fas.botherd",  "fas.tirequik", "fas.donothg",
                "fas.exhaustd", "fas.nostart",
                "fas.nothink",  "fas.nodesire", "fas.feeltird")
vars_fas_flip <- c("fas.enrgetic", "fas.concntrt")
vars_gerdq <- c("gerdq.burning",  "gerdq.regurg",  "gerdq.heartbrn", "gerdq.tums")
vars_gerdq_flip <- c("gerdq.stompain", "gerdq.nausea7")
vars_cfq <- dat_pro %>% select(cfq.reread:cfq.dksay) %>% colnames()
vars_sobq <- dat_pro %>% select(sobq.atrest:sobq.fearshbt) %>% colnames()
vars_promis <- dat_pro %>% select(c(pf.fatigue, pf.avgtire:pf.washdout, pf.limitact)) %>% colnames()
vars_promis2 <- dat_pro %>% select(c(pf.pushself, pf.difdecsn)) %>% colnames()

f_visit(dat_pro) # how many have pft from first and/or second visit data
dat_pro <- dat_pro %>% 
  filter(visit=="Enrollment / Clinic Visit")

# Fatigue assessment scale scoring
# Put values on Q1-10, flipping direction of Q4 and Q10.
dat_pro2 <- dat_pro %>% 
  mutate_at(vars_fas, function(x){
    case_when(
      x == 'Never'  ~ 1,
      x == 'Sometimes'  ~ 2,
      x == 'Regularly'  ~ 3,
      x == 'Often' ~ 4,
      x == 'Always' ~ 5
    )
  }) %>%
  mutate_at(vars_fas_flip, function(x){
    case_when(
      x == 'Never'  ~ 5,
      x == 'Sometimes'  ~ 4,
      x == 'Regularly'  ~ 3,
      x == 'Often' ~ 2,
      x == 'Always' ~ 1
    )
  }) %>%
  # GERD Questionaire
  # Assign values to questions, flipping direction of "gerdq.stompain", "gerdq.nausea7"
  mutate_at(vars_gerdq, function(x){
    case_when(
      x == '0 Days'  ~ 0,
      x == '1 Day'  ~ 1,
      x == '2 - 3 Days'  ~ 2,
      x == '4 - 7 Days' ~ 3
    )
  }) %>%
  mutate_at(vars_gerdq_flip, function(x){
    case_when(
      x == '0 Days'  ~ 3,
      x == '1 Day'  ~ 2,
      x == '2 - 3 Days'  ~ 1,
      x == '4 - 7 Days' ~ 0
    )
  }) %>% 
  # Cognitive Failure Questionairre
  # Assign values to questions
  mutate_at(vars_cfq, function(x){
    case_when(
      x == 'Never'  ~ 0,
      x == 'Very rarely'  ~ 1,
      x == 'Occasionally'  ~ 2,
      x == 'Quite often' ~ 3,
      x == 'Very often' ~ 4
    )
  }) %>% 
  # Shortness of breath Questionairre
  # Assign values to questions
  mutate_at(vars_sobq, function(x){
    case_when(
      x == '0 - None at all'  ~ 0,
      x == '1'  ~ 1,
      x == '2'  ~ 2,
      x == '3' ~ 3,
      x == '4 - Severe' ~ 4,
      x == '5 - Maximal' ~ 5
    )
  }) %>% 
  # Promis Questionairre
  # Assign values to questions
  mutate_at(vars_promis, function(x){
    case_when(
      x == 'Not at all'  ~ 1,
      x == 'A little bit'  ~ 2,
      x == 'Somewhat' ~ 3,
      x == 'Quite a bit' ~ 4,
      x == 'Very much' ~ 5
    )
  }) %>% 
  mutate_at(vars_promis2, function(x){
    case_when(
      x == 'Never'  ~ 1,
      x == 'Rarely'  ~ 2,
      x == 'Sometimes' ~ 3,
      x == 'Often' ~ 4,
      x == 'Always' ~ 5
    )
  }) %>%
  # SF-12 Questionairre 
  mutate_at("sf12.ratehlth", function(x){
    case_when(
      x == 'Excellent'  ~ 1,
      x == 'Very good'  ~ 2,
      x == 'Good' ~ 3,
      x == 'Fair' ~ 4,
      x == 'Poor' ~ 5
    )
  }) %>%
  mutate_at(c("sf12.modactiv", "sf12.sevstair"), function(x){
    case_when(
      x == 'Yes, Limited a Lot'  ~ 1,
      x == 'Yes, Limited a Little'  ~ 2,
      x == 'No, Not Limited at All' ~ 3
    )
  }) %>%
  mutate_at(c("sf12.lesslike", "sf12.kindwork", "sf12.accless", "sf12.lesscare"), function(x){
    case_when(
      x == 'None of the time'  ~ 5,
      x == 'A little of the time'  ~ 4,
      x == 'Some of the time' ~ 3,
      x == 'Most of the time' ~ 2,
      x == 'All of the time' ~ 1
    )
  }) %>%
  mutate_at("sf12.work4wks", function(x){
    case_when(
      x == 'Not at all'  ~ 1,
      x == 'A little bit'  ~ 2,
      x == 'Moderately' ~ 3,
      x == 'Quite a bit' ~ 4,
      x == 'Extremely' ~ 5
    )
  }) %>%
  # Note soc4wks coding does not follow the GRADS SF-12 form, also, I don't think they combined what they think they did. Can we check our values with theirs? 
  mutate_at(c("sf12.feltcalm", "sf12.lotenrg", "sf12.blue", "sf12.soc4wks"), function(x){
    case_when(
      x == 'None of the Time'  ~ 6,
      x == 'A Little of the Time'  ~ 5,
      x == 'Some of the Time' ~ 4,
      x == 'A Good Bit of the Time' ~ 3,
      x == 'Most of the Time' ~ 2,
      x == 'All of the Time' ~ 1
    )
  }) 

# Fix questions 9-12 to collapse "A Good Bit of the Time" (3) with "Most of the time"
dat_pro2$sf12.feltcalm <- ifelse(dat_pro2$sf12.feltcalm>2, dat_pro2$sf12.feltcalm - 1, dat_pro2$sf12.feltcalm)
dat_pro2$sf12.lotenrg <- ifelse(dat_pro2$sf12.lotenrg>2, dat_pro2$sf12.lotenrg - 1, dat_pro2$sf12.lotenrg)
dat_pro2$sf12.blue <- ifelse(dat_pro2$sf12.blue>2, dat_pro2$sf12.blue - 1, dat_pro2$sf12.blue)
dat_pro2$sf12.soc4wks <- ifelse(dat_pro2$sf12.soc4wks>2, dat_pro2$sf12.soc4wks - 1, dat_pro2$sf12.soc4wks)


# Flip questions 1, 8, 9, 10
table(dat_pro2$sf12.ratehlth)
dat_pro2$sf12.ratehlth <- 6 - dat_pro2$sf12.ratehlth
dat_pro2$sf12.work4wks <- 6 - dat_pro2$sf12.work4wks
dat_pro2$sf12.feltcalm <- 6 - dat_pro2$sf12.feltcalm
dat_pro2$sf12.lotenrg <- 6 - dat_pro2$sf12.lotenrg


# Compute raw scale scores
dat_pro2$pf <- dat_pro2$sf12.modactiv + dat_pro2$sf12.sevstair
dat_pro2$rp <- dat_pro2$sf12.lesslike + dat_pro2$sf12.kindwork
dat_pro2$bp <- dat_pro2$sf12.work4wks
dat_pro2$gh <- dat_pro2$sf12.ratehlth
dat_pro2$vt <- dat_pro2$sf12.lotenrg
dat_pro2$sf <- dat_pro2$sf12.soc4wks
dat_pro2$re <- dat_pro2$sf12.accless + dat_pro2$sf12.lesscare
dat_pro2$mh <- dat_pro2$sf12.feltcalm + dat_pro2$sf12.blue

# transform raw scale scores into percentages;
dat_pro2$pf=100*((dat_pro2$pf-2)/4)
dat_pro2$rp=100*((dat_pro2$rp-2)/8)
dat_pro2$bp=100*((dat_pro2$bp-1)/4)
dat_pro2$gh=100*((dat_pro2$gh-1)/4)
dat_pro2$vt=100*((dat_pro2$vt-1)/4)
dat_pro2$sf=100*((dat_pro2$sf-1)/4)
dat_pro2$re=100*((dat_pro2$re-2)/8)
dat_pro2$mh=100*((dat_pro2$mh-2)/8)

# standardize transformed scale scores using 1998 general US population means and SDs
dat_pro2$pf=(dat_pro2$pf-81.18122)/29.10558
dat_pro2$rp=(dat_pro2$rp-80.52856)/27.13526
dat_pro2$bp=(dat_pro2$bp-81.74015)/24.53019
dat_pro2$gh=(dat_pro2$gh-72.19795)/23.19401
dat_pro2$vt=(dat_pro2$vt-55.59090)/24.84380
dat_pro2$sf=(dat_pro2$sf-83.73973)/24.75775
dat_pro2$re=(dat_pro2$re-86.41051)/22.35543
dat_pro2$mh=(dat_pro2$mh-70.18217)/20.50597

# weight standardized scale scores by their 1990 physical and mental factor score coefficients
dat_pro2$pcspf=dat_pro2$pf*0.42402 
dat_pro2$mcspf=dat_pro2$pf*-.22999

dat_pro2$pcsrp=dat_pro2$rp*0.35119 
dat_pro2$mcsrp=dat_pro2$rp*-.12329

dat_pro2$pcsbp=dat_pro2$bp*0.31754 
dat_pro2$mcsbp=dat_pro2$bp*-.09731

dat_pro2$pcsgh=dat_pro2$gh*0.24954 
dat_pro2$mcsgh=dat_pro2$gh*-.01571

dat_pro2$pcsvt=dat_pro2$vt*0.02877 
dat_pro2$mcsvt=dat_pro2$vt*0.23534

dat_pro2$pcssf=dat_pro2$sf*-.00753 
dat_pro2$mcssf=dat_pro2$sf*0.26876

dat_pro2$pcsre=dat_pro2$re*-.19206 
dat_pro2$mcsre=dat_pro2$re*0.43407

dat_pro2$pcsmh=dat_pro2$mh*-.22069 
dat_pro2$mcsmh=dat_pro2$mh*0.48581

# compute aggregate scale scores
dat_pro2$pcs = dat_pro2$pcspf + dat_pro2$pcsrp + dat_pro2$pcsbp + dat_pro2$pcsgh + dat_pro2$pcsvt + dat_pro2$pcssf + dat_pro2$pcsre + dat_pro2$pcsmh
dat_pro2$mcs = dat_pro2$mcspf + dat_pro2$mcsrp + dat_pro2$mcsbp + dat_pro2$mcsgh + dat_pro2$mcsvt + dat_pro2$mcssf + dat_pro2$mcsre + dat_pro2$mcsmh

# standardize aggregate scale scores to t scores
dat_pro2$sf12.overall.pcs=50+(dat_pro2$pcs*10)
dat_pro2$sf12.overall.mcs=50+(dat_pro2$mcs*10)


# Sum values over all questions
dat_pro2 <- dat_pro2 %>% mutate(fas.overall = select(., all_of(c(vars_fas, vars_fas_flip))) %>% apply(1, sum, na.rm=FALSE),
                 gerdq.overall = select(., all_of(c(vars_gerdq, vars_gerdq_flip))) %>% apply(1, sum, na.rm=FALSE),
                 cfq.overall = select(., all_of(vars_cfq)) %>% apply(1, sum, na.rm=FALSE),
                 sobq.overall = select(., all_of(vars_sobq)) %>% apply(1, sum, na.rm=FALSE),
                 promis.overall = select(., all_of(c(vars_promis, vars_promis2))) %>% apply(1, sum, na.rm=FALSE))
dat_pro2 <- dat_pro2 %>% select(gradsid, fas.overall, gerdq.overall, cfq.overall, sobq.overall, promis.overall, sf12.overall.pcs, sf12.overall.mcs)



write.csv(dat_pro2,savepath_pro,row.names=F)

################ baseline
dat_baseline <- 
  dat_visit %>% 
  select(
    "gradsid","visit","demog.gender","demog.race","demog.mixed","demog.ethnic",
    "demog.school","demog.income","demog.formdate","demog.dob",
    "spe.temperat","spe.oxiperc","spe.pulse","spe.resprt","spe.pehght","spe.pewght",
    "spe.pesys","spe.pedia","spe.petemp"
  ) %>%
  filter(
    ncol(.)-2 != apply(.,1,function(x) sum(is.na(f_na_vec(x[-c(1:2)]))))
  )
f_visit(dat_baseline)
dat_baseline <-
  dat_baseline %>%
  filter(visit=="Enrollment / Clinic Visit") %>%
  dplyr::select(-visit) %>%
  mutate(
    across(all_of(c("demog.formdate","demog.dob")),as.Date),
    age_yr = difftime(demog.formdate,demog.dob,units="days"),
    bmi = (spe.pewght*0.453592)/(spe.pehght*0.0254)^2,
    temperat = mapply(f_temp,temp=spe.temperat,type=spe.petemp)
  ) %>%
  mutate(age_yr = as.numeric(age_yr)/365)

dat_baseline <- merge(
  dat_patient %>% select("gradsid","smdsc.scadding","smdsc.phengrp",
                         "smdsc.site_id","ctrad.scmodel"),
  dat_baseline %>% select("gradsid","demog.gender","demog.race",
                          "demog.mixed","demog.ethnic","demog.school",
                          "demog.income","spe.oxiperc",
                          "spe.pulse","spe.resprt","spe.pehght","spe.pewght",
                          "spe.pesys","spe.pedia","age_yr","bmi","temperat"),
  by = "gradsid",
  all = TRUE) %>%
  as_tibble() %>%
  mutate(across(all_of(colnames(.)),f_na_vec))
colnames(dat_baseline) <- sub(".*\\.","",colnames(dat_baseline)) 

write.csv(dat_baseline,savepath_baseline,row.names=F)
