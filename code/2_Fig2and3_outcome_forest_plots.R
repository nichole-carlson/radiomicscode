#############################################################
###Program name: 2_Fig2and3_outcome_forest_plots.R
###Created by:  Originally Zac Buchalski. Edited by Nichole E Carlson, 
###             and William Lippitt
###Purpose: This creates the forest plots comparing our outcomes stratified
###         by the cluster labels created in our radiomic grouping. It requires
###         the dataset created after running code files 0 and 1. The outputs
###         are two files called PFT_forestplot.png and PRO_forestplot.png.
###################################################################

############################################# Libraries
library(tidyverse) #Data Manipulation
library(car) #Type III Anova
library(ggtext) #Annotations
library(stringr) #Data Manipulation
library(purrr)

############################################# Paths
##Location of code files and output
path_repo <- "/radiomics/"
##Location of your data
path_pdrive <- "/DataLibrary/"

############################################# Data Processing
dat<-read_rds(paste0(path_pdrive,"radiomics_data/dat_clean_FVCsorted_clusters_VAS.rds"))

############################################# Fit and Summarize Models
vars_outcome <- c('BDFVC', 'BDFEV1', 'BDFEVFVC', 'BDDLCO',
                  'fas.overall', 'gerdq.overall', 'cfq.overall', 'sobq.overall', 
                  'promis.overall', 'sf12.overall.pcs', 'sf12.overall.mcs')

### Model fitting and summarizing functions
f_basemodel <- function(outcomename,adj_scadTF,adj_OberTF,dataf=dat){
  fmla <- paste0(outcomename,"~ ",
                 ifelse(adj_scadTF,"factor(SCADDING) + ",""),
                 ifelse(adj_OberTF,"Oberstein_Overall + ",""),
                 "AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT")
  return(lm(as.formula(fmla),data=dataf))
}
f_model <- function(outcomename,cluster = "FVCordered_rskc",
                    adj_scadTF,adj_OberTF,dataf=dat){
  fmla <- paste0(outcomename,"~ factor(",cluster,")",
                 ifelse(adj_scadTF," + factor(SCADDING)",""),
                 ifelse(adj_OberTF," + Oberstein_Overall",""),
                 " + AGE_YR + GENDER + factor(raceethnic) + BMI + PEHGHT")
  return(lm(as.formula(fmla),data=dataf))
}

f_manymodelsummary <- function(outcome, cluster = "FVCordered_rskc"){
  dat_temp <- dat
  dat_temp$outcome <- dat_temp[,outcome] %>% unlist()
  dat_temp$cluster <- dat_temp[,cluster] %>% unlist()
  list_mods <- list(
    f_model("outcome","cluster",F,F,dataf=dat_temp),
    f_model("outcome","cluster",T,F,dataf=dat_temp),
    f_model("outcome","cluster",F,T,dataf=dat_temp),
    f_model("outcome","cluster",T,T,dataf=dat_temp),
    f_basemodel("outcome",F,F,dataf=dat_temp),
    f_basemodel("outcome",T,F,dataf=dat_temp),
    f_basemodel("outcome",F,T,dataf=dat_temp),
    f_basemodel("outcome",T,T,dataf=dat_temp)
  )
  coeffs <- lapply(list_mods[1:4],function(mod){summary(mod)$coefficients[2:4,1]}) %>%
    do.call(c,.) %>% 
    c(rep(NA,4)) # NA for base models
  minmax_conf <- lapply(list_mods[1:4],function(mod){confint(mod)[2:4,1:2]}) %>%
    do.call(rbind,.) %>% 
    rbind(matrix(NA,4,2)) # NA for base models
  r_sqs <- sapply(list_mods[1:4],function(mod)summary(mod)$r.squared) %>%
    rep(each=3) %>% 
    c(sapply(list_mods[5:8],function(mod)summary(mod)$r.squared))
  p_vals <- lapply(list_mods[1:4],function(mod){summary(mod)$coefficients[2:4,4]}) %>%
    do.call(c,.) %>% 
    c(rep(NA,4)) # NA for base models
  overall_p_vals <- sapply(list_mods[1:4],function(mod)Anova(mod,type=3)$`Pr(>F)`[2]) %>%
    rep(each=3) %>% 
    c(sapply(list_mods[5:8],function(mod)Anova(mod,type=3)$`Pr(>F)`[2]))
  
 data.frame(outcome = outcome,
             comp = c(rep(c("2v1","3v1","4v1"),4),rep(NA,4)),
             adj_for_scad = c(rep("Cluster",3),rep("SCAD+Clust",3),rep("Ober+Clust",3),rep("SCAD+Ober+Clust",3),
                              c("Base","SCAD","Ober","SCAD+Ober")),
             coeff = coeffs, min_conf = minmax_conf[,1], max_conf = minmax_conf[,2],
             p_val = p_vals, r_sq = r_sqs, overall_p_val=overall_p_vals,
             model_type = c(rep("Cluster",12),rep("Base",4))) %>%
    return()
}

### Build summary df
plot_names <- c('FVC','FEV1','FEV1/FVC','DLCO','FAS','GERDQ','CFQ','SOBQ',
                'Promis','SF12 - Physical','SF12 - Mental')
names(plot_names) <- vars_outcome
plot_names

model.dat.all <- vars_outcome %>%
  lapply(f_manymodelsummary) %>%
  bind_rows(.id="column_label") %>%
  as_tibble() %>%
  mutate(
    sig = ifelse(p_val<0.05,TRUE,FALSE),
    lvl = as.numeric(str_sub(comp,1,1)),
    r_sq = round(r_sq,3),
    overall_p_val = ifelse(overall_p_val < 0.001,"p < 0.001",
                           paste0("p = ",as.character(round(overall_p_val,3)))),
    outcome_name = factor(plot_names[outcome],levels=plot_names),
    adj_for_scad = factor(adj_for_scad,
                          levels=c("Cluster","SCAD+Clust","Ober+Clust","SCAD+Ober+Clust",
                                   "Base","SCAD","Ober","SCAD+Ober"))
  )
saveRDS(model.dat.all,paste0(path_repo,"results/ClusterModelDat_v2.RDS"))
model.dat <- model.dat.all %>%
  filter(model_type=="Cluster") %>%
  dplyr::select(-model_type)
############################################# Forest Plots
facet_label <-
  expand.grid(adj_for_scad=c("Cluster","SCAD+Clust","Ober+Clust","SCAD+Ober+Clust"),
              outcome_name=unname(plot_names[1:4])) %>%
  as_tibble() %>%
  mutate_all(as.character) %>%
  mutate(labs=LETTERS[1:nrow(.)],
         x=apply(cbind(outcome_name,adj_for_scad),1,
                 function(v){
                   filter(model.dat,outcome_name==v[1],
                          adj_for_scad==v[2])$max_conf %>%
                     max
                 }),
         outcome_name = factor(outcome_name,levels=plot_names),
         adj_for_scad = factor(adj_for_scad,
                               levels=c("Cluster","SCAD+Clust","Ober+Clust","SCAD+Ober+Clust")))

facet_label_pro <-
  expand.grid(adj_for_scad=c("Cluster","SCAD+Clust","Ober+Clust","SCAD+Ober+Clust"),
              outcome_name=unname(plot_names[5:11])) %>%
  as_tibble() %>%
  mutate_all(as.character) %>%
  mutate(labs=LETTERS[1:nrow(.)],
         x=apply(cbind(outcome_name,adj_for_scad),1,
                 function(v){
                   filter(model.dat,outcome_name==v[1],
                          adj_for_scad==v[2])$max_conf %>%
                     max
                 }),
         outcome_name = factor(outcome_name,levels=plot_names),
         adj_for_scad = factor(adj_for_scad,
                               levels=c("Cluster","SCAD+Clust","Ober+Clust","SCAD+Ober+Clust")))

###Creating the PFT forest plot
g_forest <-
  model.dat %>%
  filter(outcome_name %in% c("FVC","FEV1","FEV1/FVC","DLCO")) %>%
  ggplot() +
  # geom_errorbar() +
  geom_pointrange(aes(x=lvl, y=coeff, ymin=min_conf, ymax=max_conf,color=sig),size=0.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Clustering Group (Compared to Group 1)") + ylab("Regression coefficient (95% CI)") +
  theme_bw() +  # use a white background
  facet_grid(adj_for_scad ~ outcome_name,scales="free_x") +
  # facet_wrap(vars(outcome_name),ncol=1,scales="free") +
  scale_x_continuous(limits = c(1,4.1))+
  scale_colour_manual(name="Significance",
                      values=c("#5ec962","#3b528b")) + 
  theme(legend.position = "bottom") +
  geom_text(data=facet_label, aes(label=labs,y=x), x = 4) +
  geom_richtext(color="black",aes(x=1.2,y=-Inf,hjust=0,
                                  label=paste0(overall_p_val,", R<sup>2</sup> = ",r_sq)))
png(filename=paste0(path_repo,"figures/PFT_SCADnOBER_forestplot_v2.png"),
    width=600,height=500)  
g_forest
dev.off()
pdf(file=paste0(path_repo,"figures/PFT_SCADnOBER_forestplot_v2.pdf"),
    width=8,height=6.7)  
g_forest
dev.off()

g_forest_pro <-
  model.dat %>%
  filter(outcome_name %in% c('FAS','GERDQ','CFQ','SOBQ','Promis','SF12 - Physical','SF12 - Mental')) %>%
  ggplot() +
  # geom_errorbar() +
  geom_pointrange(aes(x=lvl, y=coeff, ymin=min_conf, ymax=max_conf,color=sig),size=0.2) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Clustering Group (Compared to Group 0)") + ylab("Regression coefficient (95% CI)") +
  theme_bw() +  # use a white background
  facet_grid(adj_for_scad ~ outcome_name,scales="free_x") +
  # facet_wrap(vars(outcome_name),ncol=1,scales="free") +
  scale_x_continuous(limits = c(1,4.1))+
  scale_colour_manual(name="Significance",
                      values=c("#5ec962","#3b528b")) + 
  theme(legend.position = "bottom",strip.text = element_text(size = 9)) +
 # geom_text(data=facet_label_pro, aes(label=labs,y=x), x = 4) +
  geom_richtext(color="black",aes(x=1.2,y=-Inf,hjust=0,
                                  label=paste0(overall_p_val,", R<sup>2</sup> = ",r_sq)),size=2.5)
png(filename=paste0(path_repo,"figures/PRO_SCADnOBER_forestplot_v2.png"),
    width=600,height=500)  
g_forest_pro
dev.off()
pdf(file=paste0(path_repo,"figures/PRO_SCADnOBER_forestplot_v2.pdf"),
    width=8,height=6.7)  
g_forest_pro
dev.off()

############################################# Look at Model Data
f_viewmanymodels <- function(outcomename){
  mod1 <- f_model(outcomename,adj_scadTF=T,adj_OberTF=T)
  mod2 <- f_model(outcomename,adj_scadTF=F,adj_OberTF=T)
  mod3 <- f_model(outcomename,adj_scadTF=T,adj_OberTF=F)
  print("*****************************")
  print("Adjust for both:")
  print("*****************************")
  print(summary(mod1))
  print("*****************************")
  print("Type 3 Anova for both:")
  print("*****************************")
  print(Anova(mod1,type = 3))
  print("*****************************")
  print("Adjust for Ober not Scadding:")
  print("*****************************")
  print(summary(mod2))
  print("*****************************")
  print("Type 3 Anova for Ober not Scadding:")
  print("*****************************")
  print(Anova(mod2,type = 3))
  print("*****************************")
  print("Adjust for Scadding not Ober:")
  print("*****************************")
  print(summary(mod3))
  print("*****************************")
  print("Type 3 Anova for Scadding not Ober:")
  print("*****************************")
  print(Anova(mod3,type = 3))
  return(NULL)
}

###Just want variable names here for coding
vars_outcome

f_viewmanymodels("BDFVC")
f_viewmanymodels("BDFEV1")
f_viewmanymodels("BDDLCO")
f_viewmanymodels("BDFEVFVC")
f_viewmanymodels("promis.overall")

dat1 <- dat
dat1$FVCordered_rskc <- factor(dat2$FVCordered_rskc,levels = c(1,0,2,3))
dat2 <- dat
dat2$FVCordered_rskc <- factor(dat2$FVCordered_rskc,levels = c(2,0,1,3))
dat3 <- dat
dat3$FVCordered_rskc <- factor(dat3$FVCordered_rskc,levels = 3:0)

mod.ober <- f_model("Oberstein_Overall",adj_scadTF=F,adj_OberTF=F)
summary(mod.ober)
mod.ober1 <- f_model("Oberstein_Overall",adj_scadTF=F,adj_OberTF=F,dataf=dat1)
summary(mod.ober1)
mod.ober2 <- f_model("Oberstein_Overall",adj_scadTF=F,adj_OberTF=F,dataf=dat2)
summary(mod.ober2)
mod.ober3 <- f_model("Oberstein_Overall",adj_scadTF=F,adj_OberTF=F,dataf=dat3)
summary(mod.ober3)

rsq_dat <- model.dat.all %>%
  dplyr::select(outcome,adj_for_scad,r_sq) %>%
  unique()

ober.mods <- list(ref0 = mod.ober,ref1 = mod.ober1,ref2 = mod.ober2,ref3 = mod.ober3)
saveRDS(ober.mods,paste0(path_repo,"results/PredictingObersteinModels_v2.RDS"))
