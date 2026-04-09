# Evolving patterns of prevalence and management of hypertension phenotypes
# in Mexico: A two-decade analysis of nationally representative surveys

## Data Analysis: Carlos Alberto Fermin-Martinez
##                Omar Yaxmehen Bello-Chavolla
##                Alejandra Nuñez-Luna
##                Enrique C. Guerra

## Latest version of Analysis: March 2026

## For any question regarding analysis contact:
## Omar Yaxmehen Bello-Chavolla at oybello@inger.gob.mx

####----------------------- Preparation and data management ------------#### ----####
#### Package load---- ####
### Packages ###
pacman::p_load(
  readr, haven, tidyverse, data.table, gtsummary, flextable, logbin,
  ggh4x, lmtest, merDeriv, officer, survey, ggpubr, lme4, lmerTest,
  jtools, performance, cardx, sandwich, glmtoolbox, naniar, gtools,
  modelsummary)

### Working directories ###
# Potential base directories, each one for a different user
BAS_DIR <- c("~/Mi unidad (obello@facmed.unam.mx)/",
             "~/nantonio@facmed.unam.mx - Google Drive/Mi unidad/",
             "~/nantonio@facmed.unam.mx - Google Drive/My Drive/")

# Input directory where ENSANUT data sets are stored
INN_DIR <- "Datasets/ENSANUT"

# Output directory where figures and tables will be saved
OUT_DIR <- "Mexico City Cohort Study/Proyectos/2_Proyectos ENSANUT"

# Function to evaluate multiple wd's, setting the 1st existing one
set_existing_wd <- function(wds) {
  for (wd in wds) {
    if (file.exists(wd) && file.info(wd)$isdir) {
      setwd(wd); cat("Working directory set to:", wd, "\n")
      return(invisible())}}
  stop("None of the provided directories exist or are accessible.\n")}

# Set input directory
set_existing_wd(paste0(BAS_DIR, INN_DIR))
WD_INN <- getwd()

# Set output directory
set_existing_wd(paste0(BAS_DIR, OUT_DIR))
WD_OUT <- getwd()



#### ENSA    2000---- ####
setwd(WD_INN)
ens00 <- readRDS("ensanut2000_fin.rds")

#Recode age to avoid missingness
ens00$Age <- ens00$edada
#Recode smoking to avoid missingness
ens00$Smoking <- with(ens00, case_when(a2101==0~0, a2103==2~1, a2103==1~2, a2101==2~1))
#Clean blood pressure to avoid non-plausible or extreme values
ens00$SBP[with(ens00,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.00
ens00$DBP[with(ens00,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.00
ens00$SBP[with(ens00,(SBP>300) & !is.na(SBP))] <- NA
ens00$DBP[with(ens00,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens00 <- ens00 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens00_fin <- ens00
ens00_fin <- ens00_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens00_fin <- ens00_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens00_fin$categories <- (ens00_fin$SBP_cat.E*10+ens00_fin$DBP_cat.E) #Categories of measured blood pressure
ens00_fin$categories[ens00_fin$HX_HBP==1|ens00_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens00_fin$categories <- factor(ens00_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens00_fin$NEW_HTN <- (ens00_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens00_fin$PRE <- (ens00_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens00_fin$IDH <- (ens00_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens00_fin$ISH <- (ens00_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens00_fin$SDH <- (ens00_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens00_fin$DX <- (ens00_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens00_fin$HTN_fin <- (ens00_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens00_fin$HTN_goal<- ifelse((ens00_fin$DBP<80 | (ens00_fin$Age<65 & ens00_fin$SBP<130) |
                              (ens00_fin$Age>=65 & ens00_fin$SBP<140))==T,1,0)
ens00_fin$HTN_not_goal<- ifelse(ens00_fin$HTN_goal==0,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens00_fin <- ens00_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens00_fin <- ens00_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens00_fin$categories.A <- (ens00_fin$SBP_cat.A*10+ens00_fin$DBP_cat.A) #Categories of measured blood pressure
ens00_fin$categories.A[ens00_fin$HX_HBP==1|ens00_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens00_fin$categories.A <- ordered(ens00_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
 "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens00_fin$NEW_HTN.A <- (ens00_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens00_fin$PRE.A <- (ens00_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens00_fin$IDH.A <- (ens00_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens00_fin$ISH.A <- (ens00_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens00_fin$SDH.A <- (ens00_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens00_fin$DX.A <- (ens00_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens00_fin$HTN_fin.A <- (ens00_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens00_fin$HX_treated <- as.numeric((ens00_fin$DX==1&ens00_fin$TX_HBP==1))
ens00_fin$HX_untreated <- as.numeric((ens00_fin$DX==1&ens00_fin$TX_HBP==0))
ens00_fin$HTN_goal_untx<- ifelse(ens00_fin$HTN_goal==1 & ens00_fin$HX_untreated==1,1,0)

#Recode dataset
ens00_fin <- ens00_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),  
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))
ens00_fin$HX_CVD <- NA

#### ENSANUT 2006---- ####
setwd(WD_INN)
ens06 <- readRDS("ensanut2006_fin.rds")

#Recode smoking to avoid missingness
ens06$Smoking <- with(ens06, case_when(a1301==3~0, a1302==2~1, a1302==1~2, a1301==2~1))
#Clean blood pressure to avoid non-plausible or extrene values
ens06$SBP[with(ens06,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.06
ens06$DBP[with(ens06,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.06
ens06$SBP[with(ens06,(SBP>300) & !is.na(SBP))] <- NA
ens06$DBP[with(ens06,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens06 <- ens06 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens06_fin <- ens06
ens06_fin <- ens06_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens06_fin <- ens06_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens06_fin$categories <- (ens06_fin$SBP_cat.E*10+ens06_fin$DBP_cat.E) #Categories of measured blood pressure
ens06_fin$categories[ens06_fin$HX_HBP==1|ens06_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens06_fin$categories <- factor(ens06_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens06_fin$NEW_HTN <- (ens06_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens06_fin$PRE <- (ens06_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens06_fin$IDH <- (ens06_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens06_fin$ISH <- (ens06_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens06_fin$SDH <- (ens06_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens06_fin$DX <- (ens06_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens06_fin$HTN_fin <- (ens06_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens06_fin$HTN_goal<- ifelse((ens06_fin$DBP<80 | (ens06_fin$Age<65 & ens06_fin$SBP<130) |
                               (ens06_fin$Age>=65 & ens06_fin$SBP<140))==T,1,0)
ens06_fin$HTN_not_goal<- ifelse(ens06_fin$HTN_goal==0,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens06_fin <- ens06_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens06_fin <- ens06_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens06_fin$categories.A <- (ens06_fin$SBP_cat.A*10+ens06_fin$DBP_cat.A) #Categories of measured blood pressure
ens06_fin$categories.A[ens06_fin$HX_HBP==1|ens06_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens06_fin$categories.A <- ordered(ens06_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens06_fin$NEW_HTN.A <- (ens06_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens06_fin$PRE.A <- (ens06_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens06_fin$IDH.A <- (ens06_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens06_fin$ISH.A <- (ens06_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens06_fin$SDH.A <- (ens06_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens06_fin$DX.A <- (ens06_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens06_fin$HTN_fin.A <- (ens06_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens06_fin$HX_treated <- as.numeric((ens06_fin$DX==1&ens06_fin$TX_HBP==1))
ens06_fin$HX_untreated <- as.numeric((ens06_fin$DX==1&ens06_fin$TX_HBP==0))
ens06_fin$HTN_goal_untx<- ifelse(ens06_fin$HTN_goal==1 & ens06_fin$HX_untreated==1,1,0)

#Recode dataset
ens06_fin <- ens06_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))

#### ENSANUT 2012---- ####
setwd(WD_INN)
ens12 <- readRDS("ensanut2012_fin.rds")

#Clean blood pressure to avoid non-plausible or extrene values
ens12$SBP[with(ens12,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.12
ens12$DBP[with(ens12,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.12
ens12$SBP[with(ens12,(SBP>300) & !is.na(SBP))] <- NA
ens12$DBP[with(ens12,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens12 <- ens12 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens12_fin <- ens12
ens12_fin <- ens12_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens12_fin <- ens12_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens12_fin$categories <- (ens12_fin$SBP_cat.E*10+ens12_fin$DBP_cat.E) #Categories of measured blood pressure
ens12_fin$categories[ens12_fin$HX_HBP==1|ens12_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens12_fin$categories <- factor(ens12_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens12_fin$NEW_HTN <- (ens12_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens12_fin$PRE <- (ens12_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens12_fin$IDH <- (ens12_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens12_fin$ISH <- (ens12_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens12_fin$SDH <- (ens12_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens12_fin$DX <- (ens12_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens12_fin$HTN_fin <- (ens12_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens12_fin$HTN_goal<- ifelse((ens12_fin$DBP<80 | (ens12_fin$Age<65 & ens12_fin$SBP<130) |
                               (ens12_fin$Age>=65 & ens12_fin$SBP<140))==T,1,0)
ens12_fin$HTN_not_goal<- ifelse(ens12_fin$HTN_goal==0,1,0)


#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens12_fin <- ens12_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens12_fin <- ens12_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens12_fin$categories.A <- (ens12_fin$SBP_cat.A*10+ens12_fin$DBP_cat.A) #Categories of measured blood pressure
ens12_fin$categories.A[ens12_fin$HX_HBP==1|ens12_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens12_fin$categories.A <- ordered(ens12_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens12_fin$NEW_HTN.A <- (ens12_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens12_fin$PRE.A <- (ens12_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens12_fin$IDH.A <- (ens12_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens12_fin$ISH.A <- (ens12_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens12_fin$SDH.A <- (ens12_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens12_fin$DX.A <- (ens12_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens12_fin$HTN_fin.A <- (ens12_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens12_fin$HX_treated <- as.numeric((ens12_fin$DX==1&ens12_fin$TX_HBP==1))
ens12_fin$HX_untreated <- as.numeric((ens12_fin$DX==1&ens12_fin$TX_HBP==0))
ens12_fin$HTN_goal_untx<- ifelse(ens12_fin$HTN_goal==1 & ens12_fin$HX_untreated==1,1,0)


#Recode dataset
ens12_fin <- ens12_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))


#### ENSANUT 2016---- ####
setwd(WD_INN)
ens16 <- readRDS("ensanut2016_fin.rds")

#Recode age and sex to avoid missingness
ens16$Age <- ens16$edad.x; ens16 <- ens16 %>%
  filter(!is.na(Age), !is.na(Sex))
#Clean blood pressure to avoid non-plausible or extrene values
ens16$SBP[with(ens16,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.16
ens16$DBP[with(ens16,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.16
ens16$SBP[with(ens16,(SBP>300) & !is.na(SBP))] <- NA
ens16$DBP[with(ens16,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens16 <- ens16 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens16_fin <- ens16
ens16_fin <- ens16_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens16_fin <- ens16_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens16_fin$categories <- (ens16_fin$SBP_cat.E*10+ens16_fin$DBP_cat.E) #Categories of measured blood pressure
ens16_fin$categories[ens16_fin$HX_HBP==1|ens16_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens16_fin$categories <- factor(ens16_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens16_fin$NEW_HTN <- (ens16_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens16_fin$PRE <- (ens16_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens16_fin$IDH <- (ens16_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens16_fin$ISH <- (ens16_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens16_fin$SDH <- (ens16_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens16_fin$DX <- (ens16_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens16_fin$HTN_fin <- (ens16_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens16_fin$HTN_goal<- ifelse((ens16_fin$DBP<80 | (ens16_fin$Age<65 & ens16_fin$SBP<130) |
                               (ens16_fin$Age>=65 & ens16_fin$SBP<140))==T,1,0)
ens16_fin$HTN_not_goal<- ifelse(ens16_fin$HTN_goal==0,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens16_fin <- ens16_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens16_fin <- ens16_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens16_fin$categories.A <- (ens16_fin$SBP_cat.A*10+ens16_fin$DBP_cat.A) #Categories of measured blood pressure
ens16_fin$categories.A[ens16_fin$HX_HBP==1|ens16_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens16_fin$categories.A <- ordered(ens16_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens16_fin$NEW_HTN.A <- (ens16_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens16_fin$PRE.A <- (ens16_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens16_fin$IDH.A <- (ens16_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens16_fin$ISH.A <- (ens16_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens16_fin$SDH.A <- (ens16_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens16_fin$DX.A <- (ens16_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens16_fin$HTN_fin.A <- (ens16_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens16_fin$HX_treated <- as.numeric((ens16_fin$DX==1&ens16_fin$TX_HBP==1))
ens16_fin$HX_untreated <- as.numeric((ens16_fin$DX==1&ens16_fin$TX_HBP==0))
ens16_fin$HTN_goal_untx<- ifelse(ens16_fin$HTN_goal==1 & ens16_fin$HX_untreated==1,1,0)

#Recode dataset
ens16_fin <- ens16_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))


#### ENSANUT 2018---- ####
setwd(WD_INN)
ens18 <- readRDS("ensanut2018_fin.rds")

#Clean blood pressure to avoid non-plausible or extrene values
ens18$SBP[with(ens18,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.18
ens18$DBP[with(ens18,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.18
ens18$SBP[with(ens18,(SBP>300) & !is.na(SBP))] <- NA
ens18$DBP[with(ens18,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens18 <- ens18 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens18_fin <- ens18
ens18_fin <- ens18_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens18_fin <- ens18_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens18_fin$categories <- (ens18_fin$SBP_cat.E*10+ens18_fin$DBP_cat.E) #Categories of measured blood pressure
ens18_fin$categories[ens18_fin$HX_HBP==1|ens18_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens18_fin$categories <- factor(ens18_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens18_fin$NEW_HTN <- (ens18_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens18_fin$PRE <- (ens18_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens18_fin$IDH <- (ens18_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens18_fin$ISH <- (ens18_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens18_fin$SDH <- (ens18_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens18_fin$DX <- (ens18_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens18_fin$HTN_fin <- (ens18_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens18_fin$HTN_goal<- ifelse((ens18_fin$DBP<80 | (ens18_fin$Age<65 & ens18_fin$SBP<130) |
                               (ens18_fin$Age>=65 & ens18_fin$SBP<140))==T,1,0)
ens18_fin$HTN_not_goal<- ifelse(ens18_fin$HTN_goal==0,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens18_fin <- ens18_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens18_fin <- ens18_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens18_fin$categories.A <- (ens18_fin$SBP_cat.A*10+ens18_fin$DBP_cat.A) #Categories of measured blood pressure
ens18_fin$categories.A[ens18_fin$HX_HBP==1|ens18_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens18_fin$categories.A <- ordered(ens18_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens18_fin$NEW_HTN.A <- (ens18_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens18_fin$PRE.A <- (ens18_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens18_fin$IDH.A <- (ens18_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens18_fin$ISH.A <- (ens18_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens18_fin$SDH.A <- (ens18_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens18_fin$DX.A <- (ens18_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens18_fin$HTN_fin.A <- (ens18_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens18_fin$HX_treated <- as.numeric((ens18_fin$DX==1&ens18_fin$TX_HBP==1))
ens18_fin$HX_untreated <- as.numeric((ens18_fin$DX==1&ens18_fin$TX_HBP==0))
ens18_fin$HTN_goal_untx<- ifelse(ens18_fin$HTN_goal==1 & ens18_fin$HX_untreated==1,1,0)

#Recode dataset
ens18_fin <- ens18_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))


#### ENSANUT 2020---- ####
setwd(WD_INN)
ens20 <- readRDS("ensanut2020_fin.rds")

#Clean blood pressure to avoid non-plausible or extreme values
ens20$SBP[with(ens20,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.20
ens20$DBP[with(ens20,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.20
ens20$SBP[with(ens20,(SBP>300) & !is.na(SBP))] <- NA
ens20$DBP[with(ens20,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens20 <- ens20 %>% filter(Age>=20, Age<100)

#Previous diagnosis
ens20$HX_T2D <- ifelse(ens20$H0902A.x==1, 1, 0) #Diabetes
ens20$HX_HBP <- transmute(ens20, if_all(c( #Hypertension
  H0902A.x, H0902B, H0902C), ~ is.na(.) | .!=3 ))[[1]] %>% ifelse(0, 1)
ens20$HX_CVD <- transmute(ens20, if_all(c( #Cardiovascular disease
  H0902A.x, H0902B, H0902C, H0902D), ~ is.na(.) | .!=4 ))[[1]] %>% ifelse(0,1)

#Diabetes definitions
ens20$diabetes_previous <- with(ens20, ifelse(HX_T2D==1, 1, 0))
ens20$diabetes_biochem <- with(ens20, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens20$diabetes_fin <- ((ens20$diabetes_biochem + ens20$diabetes_previous)>0) %>% ifelse(1,0)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens20_fin <- ens20
ens20_fin <- ens20_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens20_fin <- ens20_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens20_fin$categories <- (ens20_fin$SBP_cat.E*10+ens20_fin$DBP_cat.E) #Categories of measured blood pressure
ens20_fin$categories[ens20_fin$HX_HBP==1|ens20_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens20_fin$categories <- factor(ens20_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens20_fin$NEW_HTN <- (ens20_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens20_fin$PRE <- (ens20_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens20_fin$IDH <- (ens20_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens20_fin$ISH <- (ens20_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens20_fin$SDH <- (ens20_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens20_fin$DX <- (ens20_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens20_fin$HTN_fin <- (ens20_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens20_fin$HTN_goal<- ifelse((ens20_fin$DBP<80 | (ens20_fin$Age<65 & ens20_fin$SBP<130) |
                               (ens20_fin$Age>=65 & ens20_fin$SBP<140))==T,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens20_fin <- ens20_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens20_fin <- ens20_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens20_fin$categories.A <- (ens20_fin$SBP_cat.A*10+ens20_fin$DBP_cat.A) #Categories of measured blood pressure
ens20_fin$categories.A[ens20_fin$HX_HBP==1|ens20_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens20_fin$categories.A <- ordered(ens20_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens20_fin$NEW_HTN.A <- (ens20_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens20_fin$PRE.A <- (ens20_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens20_fin$IDH.A <- (ens20_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens20_fin$ISH.A <- (ens20_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens20_fin$SDH.A <- (ens20_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens20_fin$DX.A <- (ens20_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens20_fin$HTN_fin.A <- (ens20_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens20_fin$HX_treated <- as.numeric((ens20_fin$DX==1&ens20_fin$TX_HBP==1))
ens20_fin$HX_untreated <- as.numeric((ens20_fin$DX==1&ens20_fin$TX_HBP==0))

#Recode dataset
ens20_fin <- ens20_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))


#### ENSANUT 2021---- ####
setwd(WD_INN)
ens21 <- readRDS("ensanut2021_fin.rds")

#Remove n=1 missing in HX_T2D and HX_CKD
ens21 <- ens21 %>% filter(!is.na(HX_T2D))
#Clean blood pressure to avoid non-plausible or extrene values
ens21$SBP[with(ens21,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.21
ens21$DBP[with(ens21,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.21
ens21$SBP[with(ens21,(SBP>300) & !is.na(SBP))] <- NA
ens21$DBP[with(ens21,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens21 <- ens21 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens21_fin <- ens21
ens21_fin <- ens21_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens21_fin <- ens21_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens21_fin$categories <- (ens21_fin$SBP_cat.E*10+ens21_fin$DBP_cat.E) #Categories of measured blood pressure
ens21_fin$categories[ens21_fin$HX_HBP==1|ens21_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens21_fin$categories <- factor(ens21_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens21_fin$NEW_HTN <- (ens21_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens21_fin$PRE <- (ens21_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens21_fin$IDH <- (ens21_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens21_fin$ISH <- (ens21_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens21_fin$SDH <- (ens21_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens21_fin$DX <- (ens21_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens21_fin$HTN_fin <- (ens21_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens21_fin$HTN_goal<- ifelse((ens21_fin$DBP<80 | (ens21_fin$Age<65 & ens21_fin$SBP<130) |
                               (ens21_fin$Age>=65 & ens21_fin$SBP<140))==T,1,0)
ens21_fin$HTN_not_goal<- ifelse(ens21_fin$HTN_goal==0,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens21_fin <- ens21_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens21_fin <- ens21_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens21_fin$categories.A <- (ens21_fin$SBP_cat.A*10+ens21_fin$DBP_cat.A) #Categories of measured blood pressure
ens21_fin$categories.A[ens21_fin$HX_HBP==1|ens21_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens21_fin$categories.A <- ordered(ens21_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens21_fin$NEW_HTN.A <- (ens21_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens21_fin$PRE.A <- (ens21_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens21_fin$IDH.A <- (ens21_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens21_fin$ISH.A <- (ens21_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens21_fin$SDH.A <- (ens21_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens21_fin$DX.A <- (ens21_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens21_fin$HTN_fin.A <- (ens21_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Recode dataset
ens21_fin <- ens21_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))

#Untreated hypertension
ens21_fin$HX_treated <- as.numeric((ens21_fin$DX==1&ens21_fin$TX_HBP==1))
ens21_fin$HX_untreated <- as.numeric((ens21_fin$DX==1&ens21_fin$TX_HBP==0))
ens21_fin$HTN_goal_untx<- ifelse(ens21_fin$HTN_goal==1 & ens21_fin$HX_untreated==1,1,0)

##EVC
ens21_fin$evc<-ifelse(ens21_fin$a0506==1,1,0)

## Falla caríaca
ens21_fin$falla_cardiaca<-ifelse(ens21_fin$a0502c==1, 1,0)
table(ens21_fin$falla_cardiaca)

#### ENSANUT 2022---- ####
setwd(WD_INN)
ens22 <- readRDS("ensanut2022_fin.rds")

#Clean blood pressure to avoid non-plausible or extrene values
ens22$SBP[with(ens22,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.22
ens22$DBP[with(ens22,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.22
ens22$SBP[with(ens22,(SBP>300) & !is.na(SBP))] <- NA
ens22$DBP[with(ens22,(DBP>SBP) & !is.na(DBP))] <- NA
#Age filters
ens22 <- ens22 %>% filter(Age>=20, Age<100)

#----------------- European Society of Cardiology, European Society of Hypertension (ESC/ESH) -------------#
#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens22_fin <- ens22
ens22_fin <- ens22_fin %>% mutate("SBP_cat.E"=cut(SBP, c(-Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens22_fin <- ens22_fin %>% mutate("DBP_cat.E"=cut(DBP, c(-Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens22_fin$categories <- (ens22_fin$SBP_cat.E*10+ens22_fin$DBP_cat.E) #Categories of measured blood pressure
ens22_fin$categories[ens22_fin$HX_HBP==1|ens22_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens22_fin$categories <- factor(ens22_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens22_fin$NEW_HTN <- (ens22_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens22_fin$PRE <- (ens22_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens22_fin$IDH <- (ens22_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens22_fin$ISH <- (ens22_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens22_fin$SDH <- (ens22_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens22_fin$DX <- (ens22_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens22_fin$HTN_fin <- (ens22_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens22_fin$HTN_goal<- ifelse((ens22_fin$DBP<80 | (ens22_fin$Age<65 & ens22_fin$SBP<130) |
                               (ens22_fin$Age>=65 & ens22_fin$SBP<140))==T,1,0)

#------------------ American College of Cardiology, American Heart Association (ACC/AHA) ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens22_fin <- ens22_fin %>% mutate("SBP_cat.A"=cut(SBP, c(-Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens22_fin <- ens22_fin %>% mutate("DBP_cat.A"=cut(DBP, c(-Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens22_fin$categories.A <- (ens22_fin$SBP_cat.A*10+ens22_fin$DBP_cat.A) #Categories of measured blood pressure
ens22_fin$categories.A[ens22_fin$HX_HBP==1|ens22_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 22 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens22_fin$categories.A <- ordered(ens22_fin$categories.A, levels=c(11,21,12,22,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens22_fin$NEW_HTN.A <- (ens22_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens22_fin$PRE.A <- (ens22_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens22_fin$IDH.A <- (ens22_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens22_fin$ISH.A <- (ens22_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens22_fin$SDH.A <- (ens22_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens22_fin$DX.A <- (ens22_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens22_fin$HTN_fin.A <- (ens22_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens22_fin$HX_treated <- as.numeric((ens22_fin$DX==1&ens22_fin$TX_HBP==1))
ens22_fin$HX_untreated <- as.numeric((ens22_fin$DX==1&ens22_fin$TX_HBP==0))
ens22_fin$HTN_goal_untx<- ifelse(ens22_fin$HTN_goal==1 & ens22_fin$HX_untreated==1,1,0)
ens22_fin$HTN_not_goal<- ifelse(ens22_fin$HTN_goal==0,1,0)
ens22_fin$HTN_not_goal_untx<- ifelse(ens22_fin$HTN_goal==0 & ens22_fin$HX_treated==0,1,0)

#Recode dataset
ens22_fin <- ens22_fin %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,25,30,Inf), right = F) %>% ordered(labels=c("Normal weight","Overweight","Obesity")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))


#### ENSANUT 2023---- ####
setwd(WD_INN)
ens23 <- readRDS("ensanut2023_fin.rds")

# ------------------ Recode ------------------ #
ens23_fin <- ens23 %>% 
  mutate(
    "ID"=FOLIO_INT,
    "Sex"=case_when(sexo==1~0, sexo==2~1),
    "Age"=edad,
    "Year"=2023,
    "Area"=case_when(estrato.x==1~0, estrato.x==2~1,
                     estrato.x==3~1), #0:Rural // 1:Urban
    "State"=entidad.x,
    "Education"=h0317a,
    "Education_cat"=case_when(h0317a%in%c(0,1)~"None",
                              h0317a%in%c(2,6)~"Primary",
                              h0317a%in%c(3,4,5,7,8)~"Secondary",
                              h0317a%in%c(9,10,11,12)~"Tertiary") %>%  
      ordered(levels=c("None", "Primary", "Secondary", "Tertiary")),
    "Education_cat2" = case_when(h0317a %in% c(0:3, 6) ~ "Primary",
                                 h0317a %in% c(4:5, 7:8) ~ "Upper secondary",
                                 h0317a %in% c(9:12) ~ "University") %>%
      ordered(levels = c("Primary", "Upper secondary", "University")),
    "Indigenous"=case_when(h0311==2~0, h0311==1~1),
    "Smoking"=case_when(a1305==3~"Never", a1301==3~"Former",
                        a1301%in%c(1,2)~"Current") %>%
      ordered(levels=c("Never", "Former", "Current")),
    "Alcohol"=case_when(a1308%in%c(5,6)~"Not currently",
                        a1308%in%c(2,3,4)~"Less than daily",
                        a1308%in%1~"Daily") %>%
      ordered(levels=c("Not currently", "Less than daily", "Daily")),
    "HX_T2D"=case_when(a0301%in%c(2,3)~0, a0301==1~1), #Diabetes
    "HX_HBP"=case_when(a0401%in%c(2,3)~0, a0401==1~1), #Hypertension
    "HX_AMI"=case_when((a0502a==2)~0, a0502a==1~1), #Myocardial infarction
    "HX_HFA"=case_when((a0502c==2)~0, a0502c==1~1), #Heart failure
    "HX_EVC"=case_when(a0502d==2~0, a0502d==1~1), #Stroke
    "HX_CKD"=case_when(a0601c==2~0, a0601c==1~1), #CKD
    "TX_T2D"=a0307, #1:Insulin // 2:Pills // 3:Both // 4:None
    "TX_HBP"=case_when(a0404==2~0, a0404==1~1), #0:No // 1:Yes
    "Glucose"=GLU_SUERO,
    "HBA1C"=HB1AC,
    "Insulin"=INSULINA,
    "TG"=TRIG,
    "CT"=COLEST,
    "HDL"=COL_HDL,
    "LDL"=COL_LDL,
    "Fasting"=san04)

#HEIGHT (cm)
ens23_fin <- ens23_fin %>% filter(Age>=20) #Filter >20 yo
ens23_fin$Height_1 <- ens23_fin$an04_1 #1st measurement
ens23_fin$Height_1[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an15_1)] <- ens23_fin$an15_1[!is.na(ens23_fin$an15_1)]
ens23_fin$Height_2 <- ens23_fin$an04_2 #2nd measurement
ens23_fin$Height_2[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an15_2)] <- ens23_fin$an15_2[!is.na(ens23_fin$an15_2)]

#Remove "not measured" values: 222.2
ens23_fin$Height_1[ens23_fin$Height_1==222.2] <- NA
ens23_fin$Height_2[ens23_fin$Height_2==222.2] <- NA
#Final height: average of both measurements
ens23_fin$Height <- ens23_fin %>% select(
  Height_1, Height_2) %>% apply(1, mean, na.rm=T)

#WEIGHT (kg)
ens23_fin$Weight_1 <- ens23_fin$an01_1 #1st measurement
ens23_fin$Weight_1[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an12_1)] <- ens23_fin$an12_1[!is.na(ens23_fin$an12_1)]
ens23_fin$Weight_2 <- ens23_fin$an01_2 #2nd measurement
ens23_fin$Weight_2[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an12_2)] <- ens23_fin$an12_2[!is.na(ens23_fin$an12_2)]
#Remove "not measured" values: 222.22
ens23_fin$Weight_1[ens23_fin$Weight_1==222.22] <- NA
ens23_fin$Weight_2[ens23_fin$Weight_2==222.22] <- NA
#Final Weight: average of both measurements
ens23_fin$Weight <- ens23_fin %>% select(
  Weight_1, Weight_2) %>% apply(1, mean, na.rm=T)

#WAIST (cm)
ens23_fin$Waist_1 <- ens23_fin$an08_1 #1st measurement
ens23_fin$Waist_1[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an21_1)] <- ens23_fin$an21_1[!is.na(ens23_fin$an21_1)]
ens23_fin$Waist_2 <- ens23_fin$an08_2 #2nd measurement
ens23_fin$Waist_2[ens23_fin$Age>=60 & !is.na(
  ens23_fin$an21_2)] <- ens23_fin$an21_2[!is.na(ens23_fin$an21_2)]
#Remove "not measured" values: 222.22
ens23_fin$Waist_1[ens23_fin$Waist_1==222.22] <- NA
ens23_fin$Waist_2[ens23_fin$Waist_2==222.22] <- NA
#Final Waist: average of both measurements
ens23_fin$Waist <- ens23_fin %>% select(
  Waist_1, Waist_2) %>% apply(1, mean, na.rm=T)

#BMI, WHtR
ens23_fin <- ens23_fin %>% mutate(
  "BMI"=Weight/(Height/100)^2, "WHtR"=Waist/Height)

#SBP/DBP
#Remove "not measured" values: 999
ens23_fin$an27_01s[ens23_fin$an27_01s==999]<-NA; ens23_fin$an27_01d[ens23_fin$an27_01d==999]<-NA
ens23_fin$an27_02s[ens23_fin$an27_02s==999]<-NA; ens23_fin$an27_02d[ens23_fin$an27_02d==999]<-NA
ens23_fin$an27_03s[ens23_fin$an27_03s==999]<-NA; ens23_fin$an27_03d[ens23_fin$an27_03d==999]<-NA
#Final SBP/DBP: average of all three measurements
ens23_fin$SBP <- ens23_fin %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ens23_fin$DBP <- ens23_fin %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)

#CVD
ens23_fin$HX_CVD <- ((ens23_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens23_fin$HX_CVD[with(ens23_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA

#Diabetes
ens23_fin$diabetes_previous <- with(ens23_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens23_fin$diabetes_biochem <- with(ens23_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens23_fin$diabetes_undx <- with(ens23_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens23_fin$diabetes_fin <- with(ens23_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))

#Age filters
ens23_fin <- ens23_fin %>% filter(Age>=20, Age<100)

#----------------- ESC/ESH definition -------------#
#Clean blood pressure to avoid non-plausible or extrene values
ens23_fin$SBP[with(ens23_fin,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.23
ens23_fin$DBP[with(ens23_fin,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.23
ens23_fin$SBP[with(ens23_fin,(SBP>300) & !is.na(SBP))] <- NA
ens23_fin$DBP[with(ens23_fin,(DBP>SBP) & !is.na(DBP))] <- NA

#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens23_fin <- ens23_fin %>% mutate("SBP_cat.E"=cut(SBP, c(
  -Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens23_fin <- ens23_fin %>% mutate("DBP_cat.E"=cut(DBP, c(
  -Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens23_fin$categories <- (ens23_fin$SBP_cat.E*10+ens23_fin$DBP_cat.E) #Categories of measured blood pressure
ens23_fin$categories[ens23_fin$HX_HBP==1|ens23_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens23_fin$categories <- factor(ens23_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens23_fin$NEW_HTN <- (ens23_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens23_fin$PRE <- (ens23_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens23_fin$IDH <- (ens23_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens23_fin$ISH <- (ens23_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens23_fin$SDH <- (ens23_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens23_fin$DX <- (ens23_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens23_fin$HTN_fin <- (ens23_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN
ens23_fin$HTN_goal<- ifelse((ens23_fin$DBP<80 | (ens23_fin$Age<65 & ens23_fin$SBP<130) |
                               (ens23_fin$Age>=65 & ens23_fin$SBP<140))==T,1,0)

#------------------ ACC/AHA definition ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens23_fin <- ens23_fin %>% mutate("SBP_cat.A"=cut(SBP, c(
  -Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens23_fin <- ens23_fin %>% mutate("DBP_cat.A"=cut(DBP, c(
  -Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens23_fin$categories.A <- (ens23_fin$SBP_cat.A*10+ens23_fin$DBP_cat.A) #Categories of measured blood pressure
ens23_fin$categories.A[ens23_fin$HX_HBP==1|ens23_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 23 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens23_fin$categories.A <- ordered(ens23_fin$categories.A, levels=c(11,21,12,23,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens23_fin$NEW_HTN.A <- (ens23_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens23_fin$PRE.A <- (ens23_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens23_fin$IDH.A <- (ens23_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens23_fin$ISH.A <- (ens23_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens23_fin$SDH.A <- (ens23_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens23_fin$DX.A <- (ens23_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens23_fin$HTN_fin.A <- (ens23_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens23_fin$HX_treated <- as.numeric((ens23_fin$DX==1&ens23_fin$TX_HBP==1))
ens23_fin$HX_untreated <- as.numeric((ens23_fin$DX==1&ens23_fin$TX_HBP==0))
ens23_fin$HTN_goal_untx<- ifelse(ens23_fin$HTN_goal==1 & ens23_fin$HX_untreated==1,1,0)
ens23_fin$HTN_not_goal<- ifelse(ens23_fin$HTN_goal==0,1,0)
ens23_fin$HTN_not_goal_untx<- ifelse(ens23_fin$HTN_goal==0 & ens23_fin$HX_treated==0,1,0)

#Recode dataset
ens23_fin <- ens23_fin %>% mutate(
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, labels=c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))

#### ENSANUT 2024---- ####
setwd(WD_INN)
ens24 <- readRDS("ensanut2024_fin.rds")

#Recode
ens24_fin <- ens24 %>%
  mutate(
    "Education_cat" = ordered(
      Education_cat, levels=c("None", "Primary", "Secondary", "Tertiary")),
    "Education_cat2" = ordered(
      Education_cat2, levels = c("Primary", "Upper secondary", "University")),
    "Smoking" = ordered(
      Smoking, levels=c("Never", "Former", "Current")),
    "Alcohol" = ordered(
      Alcohol, levels=c("Not currently", "Less than daily", "Daily")))

#Age filters
ens24_fin <- ens24_fin %>% filter(Age>=20, Age<100)

#----------------- ESC/ESH definition -------------#
#Clean blood pressure to avoid non-plausible or extrene values
ens24_fin$SBP[with(ens24_fin,(SBP>300) & !is.na(SBP))] %>% length() -> rm.SBP.24
ens24_fin$DBP[with(ens24_fin,(DBP>SBP) & !is.na(DBP))] %>% length() -> rm.DBP.24
ens24_fin$SBP[with(ens24_fin,(SBP>300) & !is.na(SBP))] <- NA
ens24_fin$DBP[with(ens24_fin,(DBP>SBP) & !is.na(DBP))] <- NA

#Normal: <130/85
#High-Normal: 130-139/85-90
#Hypertension: ≥140/90
ens24_fin <- ens24_fin %>% mutate("SBP_cat.E"=cut(SBP, c(
  -Inf, 130, 140, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens24_fin <- ens24_fin %>% mutate("DBP_cat.E"=cut(DBP, c(
  -Inf, 85, 90, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens24_fin$categories <- (ens24_fin$SBP_cat.E*10+ens24_fin$DBP_cat.E) #Categories of measured blood pressure
ens24_fin$categories[ens24_fin$HX_HBP==1|ens24_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (High-Normal), 3 (HTN)
#Both: 11(Normal), 12 (High-Normal), 13 (IDH),
#      21 (High-Normal), 22 (High-Normal), 23 (IDH),
#      31 (ISH), 32 (ISH), 33 (SDH), 100 (Diagnosed)
ens24_fin$categories <- factor(ens24_fin$categories, levels=c(11,12,13,21,22,23,31,32,33,100), labels=c(
  "Normal", "High-Normal", "IDH", "High-Normal", "High-Normal", "IDH", "ISH", "ISH", "SDH", "Diagnosed"))
ens24_fin$NEW_HTN <- (ens24_fin$categories%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens24_fin$PRE <- (ens24_fin$categories=="High-Normal") %>% ifelse(1,0) #High-Normal
ens24_fin$IDH <- (ens24_fin$categories=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens24_fin$ISH <- (ens24_fin$categories=="ISH") %>% ifelse(1,0) #Isolated systolic
ens24_fin$SDH <- (ens24_fin$categories=="SDH") %>% ifelse(1,0) #Systodiastolic
ens24_fin$DX <- (ens24_fin$categories=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens24_fin$HTN_fin <- (ens24_fin$categories%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#------------------ ACC/AHA definition ------------------#
#Normal: <120/80
#High-Normal: 120-129/<80
#Hypertension: ≥130/80
ens24_fin <- ens24_fin %>% mutate("SBP_cat.A"=cut(SBP, c(
  -Inf, 120, 130, Inf), right = F) %>% factor(labels=c(1, 2, 3)) %>% as.numeric())
ens24_fin <- ens24_fin %>% mutate("DBP_cat.A"=cut(DBP, c(
  -Inf, 80, Inf), right = F) %>% factor(labels=c(1, 2)) %>% as.numeric())
ens24_fin$categories.A <- (ens24_fin$SBP_cat.A*10+ens24_fin$DBP_cat.A) #Categories of measured blood pressure
ens24_fin$categories.A[ens24_fin$HX_HBP==1|ens24_fin$TX_HBP==1] <- 100 #Previously diagnosed hypertension
#SBP: 10 (Normal), 20 (High-Normal), 30 (HTN) /// #DBP: 1 (Normal), 2 (HTN)
#Both: 11(Normal), 12 (IDH),
#      21 (Elevated), 23 (IDH),
#      31 (ISH), 32 (SDH), 100 (Diagnosed)
ens24_fin$categories.A <- ordered(ens24_fin$categories.A, levels=c(11,21,12,23,31,32,100), labels=c(
  "Normal", "Elevated", "IDH", "IDH", "ISH", "SDH", "Diagnosed"))
ens24_fin$NEW_HTN.A <- (ens24_fin$categories.A%in%c("IDH","ISH","SDH")) %>% ifelse(1,0) #Undiagnosed hypertension
ens24_fin$PRE.A <- (ens24_fin$categories.A=="Elevated") %>% ifelse(1,0) #Elevated
ens24_fin$IDH.A <- (ens24_fin$categories.A=="IDH") %>% ifelse(1,0) #Isolated diastolic
ens24_fin$ISH.A <- (ens24_fin$categories.A=="ISH") %>% ifelse(1,0) #Isolated systolic
ens24_fin$SDH.A <- (ens24_fin$categories.A=="SDH") %>% ifelse(1,0) #Systodiastolic
ens24_fin$DX.A <- (ens24_fin$categories.A=="Diagnosed") %>% ifelse(1,0) #Diagnosed
ens24_fin$HTN_fin.A <- (ens24_fin$categories.A%in%c("Diagnosed","IDH","ISH","SDH")) %>% ifelse(1,0) #Overall HTN

#Untreated hypertension
ens24_fin$HX_treated <- as.numeric((ens24_fin$DX==1&ens24_fin$TX_HBP==1))
ens24_fin$HX_untreated <- as.numeric((ens24_fin$DX==1&ens24_fin$TX_HBP==0))

#Goal achievement
ens24_fin$HTN_goal <- ifelse(with(
  ens24_fin, DBP<80 | ((Age<65 & SBP<130) | (Age>=65 & SBP<140))), 1, 0)
ens24_fin$HTN_goal_untx<- ifelse(ens24_fin$HTN_goal==1 & ens24_fin$HX_untreated==1,1,0)
ens24_fin$HTN_not_goal<- ifelse(ens24_fin$HTN_goal==0,1,0)
ens24_fin$HTN_not_goal_untx<- ifelse(ens24_fin$HTN_goal==0 & ens24_fin$HX_treated==0,1,0)

#Recode dataset
ens24_fin <- ens24_fin %>% mutate(
  "Sex2" = ordered(Sex, 0:1, c("Men","Women")),
  "Age_cat" = cut(Age, c(-Inf,40,60,Inf), right = F) %>% ordered(labels=c("20-39","40-59","≥60")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf), right = F) %>% ordered(labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, labels=c("Non-indigenous","Indigenous")),
  "Diabetes2" = ordered(diabetes_fin, 0:1, c("No Diabetes","Diabetes")))

#### SLI & DISLI----- ####
## Load DISLI databases
setwd(WD_INN)
Density <- readxl::read_xlsx("Bases/DISLI/Densidad.xlsx")[-1,]; names(
  Density) <- c(names(Density)[1:2], paste0("x", names(Density[-(1:2)])))

#2000
DISLI_2000 <- (readxl::read_xlsx("Bases/DISLI/IRS_2000.xlsx") %>%
                 `names<-`(c("State","St_name","SLI","SLI_cat")))[-1,]
DISLI_2000 <- DISLI_2000 %>% mutate("dens"=Density$x2000) %>%
  na.omit() %>% mutate("DISLI"=lm(SLI~dens) %>% residuals) %>% mutate(
    "DISLI_cat"=cut(DISLI, breaks=c(-Inf,quantile(DISLI,c(.25,.5,.75)),Inf),
                    labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat%in%c("Very-Low","Low","Moderate"),0,1))

#2005
DISLI_2005 <- (readxl::read_xlsx("Bases/DISLI/IRS_2005.xlsx") %>%
                 `names<-`(c("State","St_name","SLI","SLI_cat")))[-1,]
DISLI_2005 <- DISLI_2005 %>% mutate("dens"=Density$x2005) %>%
  na.omit() %>% mutate("DISLI"=lm(SLI~dens) %>% residuals) %>% mutate(
    "DISLI_cat"=cut(DISLI, breaks=c(-Inf,quantile(DISLI,c(.25,.5,.75)),Inf),
                    labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat%in%c("Very-Low","Low","Moderate"),0,1))

#2010
DISLI_2010 <- (readxl::read_xlsx("Bases/DISLI/IRS_2010.xlsx") %>%
                 `names<-`(c("State","St_name","SLI","SLI_cat")))[-1,]
DISLI_2010 <- DISLI_2010 %>% mutate("dens"=Density$x2010) %>%
  na.omit() %>% mutate("DISLI"=lm(SLI~dens) %>% residuals) %>% mutate(
    "DISLI_cat"=cut(DISLI, breaks=c(-Inf,quantile(DISLI,c(.25,.5,.75)),Inf),
                    labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat%in%c("Very-Low","Low","Moderate"),0,1))

#2015
DISLI_2015 <- (readxl::read_xlsx("Bases/DISLI/IRS_2015.xlsx") %>%
                 `names<-`(c("State","St_name","SLI","SLI_cat")))[-1,]
DISLI_2015 <- DISLI_2015 %>% mutate("dens"=Density$x2015) %>%
  na.omit() %>% mutate("DISLI"=lm(SLI~dens) %>% residuals) %>% mutate(
    "DISLI_cat"=cut(DISLI, breaks=c(-Inf,quantile(DISLI,c(.25,.5,.75)),Inf),
                    labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat%in%c("Very-Low","Low","Moderate"),0,1))

#2020
DISLI_2020 <- (readxl::read_xlsx("Bases/DISLI/IRS_2020.xlsx") %>%
                 `names<-`(c("State","St_name","SLI","SLI_cat")))[-1,]
DISLI_2020 <- DISLI_2020 %>% mutate("dens"=Density$x2020) %>%
  na.omit() %>% mutate("DISLI"=lm(SLI~dens) %>% residuals) %>% mutate(
    "DISLI_cat"=cut(DISLI, breaks=c(-Inf,quantile(DISLI,c(.25,.5,.75)),Inf),
                    labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat%in%c("Very-Low","Low","Moderate"),0,1))

#Remove previously created DISLI variables
ens00_fin$DISLI<-NULL; ens00_fin$DISLI_cat<-NULL; ens00_fin$DISLI_cat2<-NULL
ens06_fin$DISLI<-NULL; ens06_fin$DISLI_cat<-NULL; ens06_fin$DISLI_cat2<-NULL
ens12_fin$DISLI<-NULL; ens12_fin$DISLI_cat<-NULL; ens12_fin$DISLI_cat2<-NULL
ens16_fin$DISLI<-NULL; ens16_fin$DISLI_cat<-NULL; ens16_fin$DISLI_cat2<-NULL
ens18_fin$DISLI<-NULL; ens18_fin$DISLI_cat<-NULL; ens18_fin$DISLI_cat2<-NULL
ens20_fin$DISLI<-NULL; ens20_fin$DISLI_cat<-NULL; ens20_fin$DISLI_cat2<-NULL
ens21_fin$DISLI<-NULL; ens21_fin$DISLI_cat<-NULL; ens21_fin$DISLI_cat2<-NULL
ens22_fin$DISLI<-NULL; ens22_fin$DISLI_cat<-NULL; ens22_fin$DISLI_cat2<-NULL
ens23_fin$DISLI<-NULL; ens23_fin$DISLI_cat<-NULL; ens23_fin$DISLI_cat2<-NULL
ens24_fin$DISLI<-NULL; ens24_fin$DISLI_cat<-NULL; ens24_fin$DISLI_cat2<-NULL

#Make State a numeric variable in all ENSANUT and DISLI cycles
ens00_fin$State <- as.numeric(ens00_fin$State)
ens06_fin$State <- as.numeric(ens06_fin$State)
ens12_fin$State <- as.numeric(ens12_fin$State)
ens16_fin$State <- as.numeric(ens16_fin$State)
ens18_fin$State <- as.numeric(ens18_fin$State)
ens20_fin$State <- as.numeric(ens20_fin$State)
ens21_fin$State <- as.numeric(ens21_fin$State)
ens22_fin$State <- as.numeric(ens22_fin$State)
ens23_fin$State <- as.numeric(ens23_fin$State)
ens24_fin$State <- as.numeric(ens24_fin$State)

DISLI_2000$State <- as.numeric(DISLI_2000$State)
DISLI_2005$State <- as.numeric(DISLI_2005$State)
DISLI_2010$State <- as.numeric(DISLI_2010$State)
DISLI_2015$State <- as.numeric(DISLI_2015$State)
DISLI_2020$State <- as.numeric(DISLI_2020$State)

#Merge all years by State
ens00_fin <- merge(ens00_fin, DISLI_2000[c(1,6:8)], by="State", all.x = T)
ens06_fin <- merge(ens06_fin, DISLI_2005[c(1,6:8)], by="State", all.x = T)
ens12_fin <- merge(ens12_fin, DISLI_2010[c(1,6:8)], by="State", all.x = T)
ens16_fin <- merge(ens16_fin, DISLI_2015[c(1,6:8)], by="State", all.x = T)
ens18_fin <- merge(ens18_fin, DISLI_2015[c(1,6:8)], by="State", all.x = T)
ens20_fin <- merge(ens20_fin, DISLI_2020[c(1,6:8)], by="State", all.x = T)
ens21_fin <- merge(ens21_fin, DISLI_2020[c(1,6:8)], by="State", all.x = T)
ens22_fin <- merge(ens22_fin, DISLI_2020[c(1,6:8)], by="State", all.x = T)
ens23_fin <- merge(ens23_fin, DISLI_2020[c(1,6:8)], by="State", all.x = T)
ens24_fin <- merge(ens24_fin, DISLI_2020[c(1,6:8)], by="State", all.x = T)


#### Weights--------- ####
# ENSA 2000
ens00_fin.2 <- ens00_fin %>% mutate(
  "pond.ant"=pond_ant, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens00_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens00_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq); ens00_fin.2 <- merge(
    ens00_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens00_fin.2$swts <- ens00_fin.2$pond.ant*(ens00_fin.2$jn/ens00_fin.2$sumwt)

# ENSANUT 2006
ens06_fin.2 <- ens06_fin %>% mutate(
  "pond.ant"=pond_ant, "estr.ant"=estrato.x,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens06_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens06_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens06_fin.2 <- merge(ens06_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens06_fin.2$swts <- ens06_fin.2$pond.ant*(ens06_fin.2$jn/ens06_fin.2$sumwt)

# ENSANUT 2012
ens12_fin.2 <- ens12_fin %>% mutate(
  "pond.ant"=pond_hta, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens12_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens12_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens12_fin.2 <- merge(ens12_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens12_fin.2$swts <- ens12_fin.2$pond.ant*(ens12_fin.2$jn/ens12_fin.2$sumwt)

# ENSANUT 2016
ens16_fin.2 <- ens16_fin %>% mutate(
  "pond.ant"=ponde_f.x.x, "estr.ant"=est_var.x,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens16_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens16_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens16_fin.2 <- merge(ens16_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens16_fin.2$swts <- ens16_fin.2$pond.ant*(ens16_fin.2$jn/ens16_fin.2$sumwt)

# ENSANUT 2018
ens18_fin.2 <- ens18_fin %>% mutate(
  "pond.ant"=F_ANTROP_INSP, "estr.ant"=EST_DIS.x,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens18_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens18_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens18_fin.2 <- merge(ens18_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens18_fin.2$swts <- ens18_fin.2$pond.ant*(ens18_fin.2$jn/ens18_fin.2$sumwt)


# ENSANUT 2020
ens20_fin.2 <- ens20_fin %>% mutate(
  "pond.ant"=pond_ant, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens20_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens20_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens20_fin.2 <- merge(ens20_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens20_fin.2$swts <- ens20_fin.2$pond.ant*(ens20_fin.2$jn/ens20_fin.2$sumwt)

# ENSANUT 2021
ens21_fin.2 <- ens21_fin %>% mutate(
  "pond.ant"=pond_ant, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens21_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens21_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens21_fin.2 <- merge(ens21_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens21_fin.2$swts <- ens21_fin.2$pond.ant*(ens21_fin.2$jn/ens21_fin.2$sumwt)

# ENSANUT 2022
ens22_fin.2 <- ens22_fin %>% mutate(
  "pond.ant"=ponde_f.x.x, "estr.ant"=est_sel.x.x,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens22_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens22_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens22_fin.2 <- merge(ens22_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens22_fin.2$swts <- ens22_fin.2$pond.ant*(ens22_fin.2$jn/ens22_fin.2$sumwt)

# ENSANUT 2023
ens23_fin.2 <- ens23_fin %>% mutate(
  "pond.ant"=pond_hta, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens23_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens23_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens23_fin.2 <- merge(ens23_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens23_fin.2$swts <- ens23_fin.2$pond.ant*(ens23_fin.2$jn/ens23_fin.2$sumwt)

# ENSANUT 2024
ens24_fin.2 <- ens24_fin %>% mutate(
  "pond.ant"=pond_hta, "estr.ant"=est_final,
  "sampleid"=paste(State, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant)
  ); wts <- with(ens24_fin.2, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens24_fin.2$sampleid)); wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens24_fin.2 <- merge(ens24_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens24_fin.2$swts <- ens24_fin.2$pond.ant*(ens24_fin.2$jn/ens24_fin.2$sumwt)



#### Recoding-------- ####
#Alcohol intake
ens00_fin.2$Alc <- with(ens00_fin.2, case_when(a2107==0~"A", a2108==2~"B", a2109==2~"C", a2109==1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens06_fin.2$Alc <- with(ens06_fin.2, case_when(a1305==0~"A", a1305==2~"B", a1306a%in%2:4~"C", a1306a==1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens12_fin.2$Alc <- with(ens12_fin.2, case_when(a1310==0~"A", a1311==13~"B", a1311%in%5:12~"C", a1311%in%1:4~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens16_fin.2$Alc <- with(ens16_fin.2, case_when(a1312c==98~"A", a1312c%in%1:24~"C", a1312c%in%25:30~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens18_fin.2$Alc <- with(ens18_fin.2, case_when(P13_11==3~"A", P13_11==2~"B", P13_12_1%in%2:4~"C", P13_12_1%in%1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens20_fin.2$Alc <- with(ens20_fin.2, case_when(ADUL1A06==5~"B", ADUL1A06%in%2:4~"C", ADUL1A06%in%1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens21_fin.2$Alc <- with(ens21_fin.2, case_when(a1308==6~"A", a1308==5~"B", a1308%in%2:4~"C", a1308%in%1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens22_fin.2$Alc <- with(ens22_fin.2, case_when(a1308==6~"A", a1308==5~"B", a1308%in%2:4~"C", a1308%in%1~"D")) %>% 
  ordered(levels=c("A","B","C","D"), labels=c("Not currently", "Not currently", "Less than daily", "Daily"))
ens23_fin.2$Alc <- ens23_fin.2$Alcohol
ens24_fin.2$Alc <- ens24_fin.2$Alcohol

#Education level
#2000 (h211_2) - 0,1: None // 2: Primary // 3,4,5,6: Secondary // 10,11: Tertiary
#2006 (h216g.x) - 0,1: None // 2: Primary // 3,4,5,6,7,8,9: Secondary // 10,11,12,13: Tertiary
#2012 (h218a) - 0,1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#2016 (h218a) - 0,1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#2018 (NIVEL) - 0,1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#2020 (H0317A) - 0,1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#2021 (h0317a) - h0316==10, 1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#2022 (h0317a) - h0316==10, 1: None // 2,6: Primary // 3,4,5,7,8: Secondary // 9,10,11,12: Tertiary
#0: None // 1: Primary // 2: Secondary // 3: Tertiary
ens00_fin.2$Edu <- with(ens00_fin.2, case_when(h211_2%in%0:1~"None", h211_2==2~"Primary", h211_2%in%3:6~"Secondary", h211_2%in%10:11~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens06_fin.2$Edu <- with(ens06_fin.2, case_when(h216g.x%in%0:1~"None", h216g.x==2~"Primary", h216g.x%in%3:9~"Secondary", h216g.x%in%10:13~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens12_fin.2$Edu <- with(ens12_fin.2, case_when(h218a%in%0:1~"None", h218a%in%c(2,6)~"Primary", h218a%in%c(3:5,7:8)~"Secondary", h218a%in%9:12~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens16_fin.2$Edu <- with(ens16_fin.2, case_when(h218a%in%0:1~"None", h218a%in%c(2,6)~"Primary", h218a%in%c(3:5,7:8)~"Secondary", h218a%in%9:12~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens18_fin.2$NIVEL <- as.numeric(ens18_fin.2$NIVEL)
ens18_fin.2$Edu <- with(ens18_fin.2, case_when(NIVEL%in%0:1~"None", NIVEL%in%c(2,6)~"Primary", NIVEL%in%c(3:5,7:8)~"Secondary", NIVEL%in%9:12~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens20_fin.2$Edu <- with(ens20_fin.2, case_when(H0317A%in%0:1~"None", H0317A%in%c(2,6)~"Primary", H0317A%in%c(3:5,7:8)~"Secondary", H0317A%in%9:12~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens21_fin.2$Edu <- with(ens21_fin.2, case_when(h0317a%in%0:1~"None", h0317a%in%c(2,6)~"Primary", h0317a%in%c(3:5,7:8)~"Secondary", h0317a%in%9:12~"Tertiary")) %>%  
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens22_fin.2$Edu <- with(ens22_fin.2, case_when(h0317a%in%0:1~"None", h0317a%in%c(2,6)~"Primary", h0317a%in%c(3:5,7:8)~"Secondary", h0317a%in%9:12~"Tertiary")) %>% 
  ordered(levels=c("None", "Primary", "Secondary", "Tertiary"))
ens23_fin.2$Edu <- ens23_fin.2$Education_cat
ens24_fin.2$Edu <- ens24_fin.2$Education_cat

#Unemployment: 0: employed // 1: unemployed
#ens06_fin.2$Job <- with(ens06_fin.2, ifelse(h218.x==1, "Employed", "Unemployed"))
#ens12_fin.2$Job <- with(ens12_fin.2, ifelse(h221==1, "Employed", "Unemployed"))
#ens16_fin.2$Job <- with(ens16_fin.2, ifelse(h221==1, "Employed", "Unemployed"))
#ens20_fin.2$Job <- with(ens20_fin.2, ifelse(H0321==1, "Employed", "Unemployed"))
#ens21_fin.2$Job <- with(ens21_fin.2, ifelse(h0321==1, "Employed", "Unemployed"))
#ens22_fin.2$Job <- with(ens22_fin.2, ifelse(h0321==1, "Employed", "Unemployed"))
#ens23_fin.2$Job <- with(ens23_fin.2, ifelse(h0321==1, "Employed", "Unemployed"))

#Affiliation to health services:
#2000 (h302_1) - 1:IMSS | 2:ISSSTE | 3:Privado trabajo | 4:Privado personal | 5:PEMEX | 6:SEDENA | 10:SEMAR | 11:Estatal | 12: No tiene | 77: Otra | 88/99:NA
ens00_fin.2$SocSec <- with(ens00_fin.2, case_when(h302_1%in%3:4~"Private", h302_1%in%c(1:2,5:11,77)~"Public", h302_1==12~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2006 (h208a.x) - 1:IMSS | 2:SSA (Popular) | 3:ISSSTE estatal | 4:ISSSTE | 5:SEDENA/SEMAR | 6:PEMEX | 7: Privado | 11:No | 77: Otra | 88/99:NA
ens06_fin.2$SocSec <- with(ens06_fin.2, case_when(h208a.x%in%7~"Private", h208a.x%in%c(1:6,77)~"Public", h208a.x==11~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2012 (h211a) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5:SEDENA/SEMAR | 6:SSA (Popular, nueva generación) | 7: Privado | 8:Otra | 9:No | 99:NA
ens12_fin.2$SocSec <- with(ens12_fin.2, case_when(h211a%in%7~"Private", h211a%in%c(1:6,8)~"Public", h211a==9~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2016 (h211a) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5:SEDENA/SEMAR | 6:SSA (Popular, nueva generación) | 7: Privado | 8:Otra | 9:No | 99:NA
ens16_fin.2$SocSec <- with(ens16_fin.2, case_when(h211a%in%7~"Private", h211a%in%c(1:6,8)~"Public", h211a==9~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2018 (P3_10_[02-11,99]) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 7-8:SSA (Popular, PROSPERA) | 9: Privado | 10:Otra | 11:No | 99:NA
ens18_fin.2$SS_public <- with(ens18_fin.2, ifelse((
  P3_10_01==1|P3_10_02==1|P3_10_03==1|P3_10_04==1|P3_10_05==1|P3_10_06==1|P3_10_07==1|P3_10_08==1|P3_10_10==1),1,0))
ens18_fin.2$SocSec <- with(ens18_fin.2, case_when(P3_10_09==1~"Private", SS_public==1~"Public", P3_10_11==1~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2020 (H0310A) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 7:SSA (IMSS Bienestar) | 8: Privado | 9:Otra | 10-11:No | 99:NA
ens20_fin.2$H0310A <- as.numeric(ens20_fin.2$H0310A); ens21_fin.2$H0310A <- as.numeric(ens21_fin.2$H0310A); ens22_fin.2$H0310A <- as.numeric(ens22_fin.2$H0310A)
ens20_fin.2$SocSec <- with(ens20_fin.2, case_when(H0310A%in%8~"Private", H0310A%in%c(1:7,9)~"Public", H0310A%in%10:11~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2021 (H0310A) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 7:SSA (IMSS Bienestar) | 8: Privado | 9:Otra | 10-11:No | 99:NA
ens21_fin.2$SocSec <- with(ens21_fin.2, case_when(H0310A%in%8~"Private", H0310A%in%c(1:7,9)~"Public", H0310A%in%10:11~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2022 (H0310A) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 7:SSA (IMSS Bienestar) | 8: Privado | 9:Otra | 10-11:No | 99:NA
ens22_fin.2$SocSec <- with(ens22_fin.2, case_when(H0310A%in%8~"Private", H0310A%in%c(1:7,9)~"Public", H0310A%in%10:11~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2023 (H0310A) - 1:IMSS | 2:ISSSTE | 3:ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 7:SSA (IMSS Bienestar) | 8: Privado | 9:Otra | 10-11:No | 99:NA
ens23_fin.2$SocSec <- with(ens23_fin.2, case_when(H0310A%in%8~"Private", H0310A%in%c(1:7,9)~"Public", H0310A%in%10:11~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))
#2024 (H0310A) - 1:IMSS | 2:ISSSTE/ISSSTE estatal | 4:PEMEX | 5-6:SEDENA/SEMAR | 8: Privado | 9:Otra | 10:Ninguno | 11:IMSS-Bienestar | 99:NA
ens24_fin.2$SocSec <- with(ens24_fin.2, case_when(H0310A%in%8~"Private", H0310A%in%c(1:6,9,11)~"Public", H0310A%in%10~"Not affiliated")) %>%
  ordered(levels=c("Not affiliated","Public", "Private"))

#0: Not affiliated // 1: Affiliated
ens00_fin.2$SocSec2 <- with(ens00_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens06_fin.2$SocSec2 <- with(ens06_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens12_fin.2$SocSec2 <- with(ens12_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens16_fin.2$SocSec2 <- with(ens16_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens18_fin.2$SocSec2 <- with(ens18_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens20_fin.2$SocSec2 <- with(ens20_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens21_fin.2$SocSec2 <- with(ens21_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens22_fin.2$SocSec2 <- with(ens22_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens23_fin.2$SocSec2 <- with(ens23_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))
ens24_fin.2$SocSec2 <- with(ens24_fin.2, ifelse(SocSec=="Not affiliated","Not affiliated","Affiliated")) %>% ordered(levels=c("Affiliated", "Not affiliated"))

#Previous DM diagnosis:
ens00_fin.2$Diabetes3 <- with(ens00_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens06_fin.2$Diabetes3 <- with(ens06_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens12_fin.2$Diabetes3 <- with(ens12_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens16_fin.2$Diabetes3 <- with(ens16_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens18_fin.2$Diabetes3 <- with(ens18_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens20_fin.2$Diabetes3 <- with(ens20_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens21_fin.2$Diabetes3 <- with(ens21_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens22_fin.2$Diabetes3 <- with(ens22_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens23_fin.2$Diabetes3 <- with(ens23_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))
ens24_fin.2$Diabetes3 <- with(ens24_fin.2, ifelse(HX_T2D==1, "Diagnosed", "Non-diagosed")) %>% ordered(levels=c("Non-diagosed", "Diagnosed"))

#Area:
ens00_fin.2$Area2 <- with(ens00_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens06_fin.2$Area2 <- with(ens06_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens12_fin.2$Area2 <- with(ens12_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens16_fin.2$Area2 <- with(ens16_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens18_fin.2$Area2 <- with(ens18_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens20_fin.2$Area2 <- with(ens20_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens21_fin.2$Area2 <- with(ens21_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens22_fin.2$Area2 <- with(ens22_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens23_fin.2$Area2 <- with(ens23_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))
ens24_fin.2$Area2 <- with(ens24_fin.2, ifelse(Area==1, 1,2)) %>% ordered(levels=1:2, labels = c("Urban", "Rural"))

#DISLI:
ens00_fin.2$DISLI_cat3 <- ens00_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens06_fin.2$DISLI_cat3 <- ens06_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens12_fin.2$DISLI_cat3 <- ens12_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens16_fin.2$DISLI_cat3 <- ens16_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens18_fin.2$DISLI_cat3 <- ens18_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens20_fin.2$DISLI_cat3 <- ens20_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens21_fin.2$DISLI_cat3 <- ens21_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens22_fin.2$DISLI_cat3 <- ens22_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens23_fin.2$DISLI_cat3 <- ens23_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))
ens24_fin.2$DISLI_cat3 <- ens24_fin.2$DISLI_cat %>% ordered(levels=c("Very-Low","Low","Moderate","High"))


#### Svey designs---- ####
#Adjust options for single PSU
options(survey.adjust.domain.lonely = TRUE, survey.lonely.psu="adjust")

#---------------------------------- 2000 ----------------------------------#
ens00_fin.2$psu_fin <- paste0(ens00_fin.2$ent.x, ens00_fin.2$mun.x)
ens00_fin.2$wgt_fin <- ens00_fin.2$pond.ant
ens00_fin.2$str_fin <- paste0(ens00_fin.2$ent.x, ens00_fin.2$estr.ant)

#Sampling design
ens00_fin_survey0 <- svydesign(data=ens00_fin.2, ids=~psu_fin, #PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata
#Survey (main): all participants with complete BP data
ens00_fin_survey <- subset(ens00_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens00_fin.2 <- ens00_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens00_fin_survey2 <- subset(ens00_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens00_fin_survey3 <- subset(ens00_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens00_fin_survey4 <- subset(ens00_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens00_fin_survey5 <- subset(ens00_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens00_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2006 ----------------------------------#
ens06_fin.2$psu_fin <- ens06_fin.2$code_upm.x
ens06_fin.2$wgt_fin <- ens06_fin.2$pond.ant
ens06_fin.2$str_fin <- ens06_fin.2$estr.ant

#Sampling design
ens06_fin_survey0 <- svydesign(data=ens06_fin.2, ids=~code_upm.x, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens06_fin_survey <- subset(ens06_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens06_fin.2 <- ens06_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens06_fin_survey2 <- subset(ens06_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens06_fin_survey3 <- subset(ens06_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens06_fin_survey4 <- subset(ens06_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens06_fin_survey5 <- subset(ens06_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens06_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2012 ----------------------------------#
ens12_fin.2$psu_fin <- ens12_fin.2$code_upm.x
ens12_fin.2$wgt_fin <- ens12_fin.2$pond.ant
ens12_fin.2$str_fin <- ens12_fin.2$estr.ant

#Sampling design
ens12_fin_survey0 <- svydesign(data=ens12_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens12_fin_survey <- subset(ens12_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens12_fin.2 <- ens12_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens12_fin_survey2 <- subset(ens12_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens12_fin_survey3 <- subset(ens12_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens12_fin_survey4 <- subset(ens12_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens12_fin_survey5 <- subset(ens12_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens12_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2016 ----------------------------------#
ens16_fin.2$psu_fin <- ens16_fin.2$code_upm.x
ens16_fin.2$wgt_fin <- ens16_fin.2$pond.ant
ens16_fin.2$str_fin <- ens16_fin.2$estr.ant

#Sampling design
ens16_fin_survey0 <- svydesign(data=ens16_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens16_fin_survey <- subset(ens16_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens16_fin.2 <- ens16_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens16_fin_survey2 <- subset(ens16_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens16_fin_survey3 <- subset(ens16_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens16_fin_survey4 <- subset(ens16_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens16_fin_survey5 <- subset(ens16_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens16_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2018 ----------------------------------#
ens18_fin.2$psu_fin <- ens18_fin.2$UPM_ant
ens18_fin.2$wgt_fin <- ens18_fin.2$pond.ant
ens18_fin.2$str_fin <- ens18_fin.2$estr.ant

#Sampling design
ens18_fin_survey0 <- svydesign(data=ens18_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens18_fin_survey <- subset(ens18_fin_survey0, !is.na(SBP)&!is.na(DBP)&wgt_fin>0)
ens18_fin.2 <- ens18_fin.2 %>% filter(!is.na(SBP), !is.na(DBP), wgt_fin>0)
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens18_fin_survey2 <- subset(ens18_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens18_fin_survey3 <- subset(ens18_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens18_fin_survey4 <- subset(ens18_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens18_fin_survey5 <- subset(ens18_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens18_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2020 ----------------------------------#
ens20_fin.2$psu_fin <- ens20_fin.2$Upm.x
ens20_fin.2$wgt_fin <- ens20_fin.2$pond.ant
ens20_fin.2$str_fin <- ens20_fin.2$estr.ant

#Sampling design
ens20_fin_survey0 <- svydesign(data=ens20_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens20_fin_survey <- subset(ens20_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens20_fin.2 <- ens20_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens20_fin_survey2 <- subset(ens20_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens20_fin_survey3 <- subset(ens20_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens20_fin_survey4 <- subset(ens20_fin_survey2, DX==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens20_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2021 ----------------------------------#
ens21_fin.2$psu_fin <- ens21_fin.2$upm.x
ens21_fin.2$wgt_fin <- ens21_fin.2$pond.ant
ens21_fin.2$str_fin <- ens21_fin.2$estr.ant

#Sampling design
ens21_fin_survey0 <- svydesign(data=ens21_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens21_fin_survey <- subset(ens21_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens21_fin.2 <- ens21_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens21_fin_survey2 <- subset(ens21_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens21_fin_survey3 <- subset(ens21_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens21_fin_survey4 <- subset(ens21_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens21_fin_survey5 <- subset(ens21_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens21_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2022 ----------------------------------#
ens22_fin.2$psu_fin <- ens22_fin.2$upm.x
ens22_fin.2$wgt_fin <- ens22_fin.2$pond.ant
ens22_fin.2$str_fin <- ens22_fin.2$estr.ant

#Sampling design
ens22_fin_survey0 <- svydesign(data=ens22_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens22_fin_survey <- subset(ens22_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens22_fin.2 <- ens22_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens22_fin_survey2 <- subset(ens22_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens22_fin_survey3 <- subset(ens22_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens22_fin_survey4 <- subset(ens22_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens22_fin_survey5 <- subset(ens22_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens22_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2023 ----------------------------------#
ens23_fin.2$psu_fin <- ens23_fin.2$upm.x.x
ens23_fin.2$wgt_fin <- ens23_fin.2$pond.ant
ens23_fin.2$str_fin <- ens23_fin.2$estr.ant

#Sampling design
ens23_fin_survey0 <- svydesign(data=ens23_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens23_fin_survey <- subset(ens23_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens23_fin.2 <- ens23_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens23_fin_survey2 <- subset(ens23_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens23_fin_survey3 <- subset(ens23_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens23_fin_survey4 <- subset(ens23_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens23_fin_survey5 <- subset(ens23_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens23_fin_survey, na.rm = T)*100 #Undiagnosed


#---------------------------------- 2024 ----------------------------------#
ens24_fin.2$psu_fin <- ens24_fin.2$upm.x.x
ens24_fin.2$wgt_fin <- ens24_fin.2$pond.ant
ens24_fin.2$str_fin <- ens24_fin.2$estr.ant

#Sampling design
ens24_fin_survey0 <- svydesign(data=ens24_fin.2, ids=~psu_fin, #Data & PSU
                               weights=~wgt_fin, #Anthropometry & BP weights
                               strata=~str_fin, nest=TRUE) #Strata 1st stage
#Survey (main): all participants with complete BP data
ens24_fin_survey <- subset(ens24_fin_survey0, !is.na(SBP)&!is.na(DBP))
ens24_fin.2 <- ens24_fin.2 %>% filter(!is.na(SBP), !is.na(DBP))
#Survey 2: participants with hypertension (compare diagnosed/UDH)
ens24_fin_survey2 <- subset(ens24_fin_survey, HTN_fin==1)
#Survey 3: participants with undiagnosed hypertension (compare ISH/IDH/SDH)
ens24_fin_survey3 <- subset(ens24_fin_survey2, NEW_HTN==1)
#Survey 4: participants with diagnosed hypertension (compare treated/UHT)
ens24_fin_survey4 <- subset(ens24_fin_survey2, DX==1)
#Survey 5: participants with treated hypertension (compare goals/not goals)
ens24_fin_survey5 <- subset(ens24_fin_survey4, HX_treated==1)
#HBP prevalence
svymean(~HTN_fin+NEW_HTN+DX, ens24_fin_survey, na.rm = T)*100 #Undiagnosed
svymean(~HX_treated+HTN_goal, ens24_fin_survey4, na.rm = T)*100 #Undiagnosed


#### Flowchart (SF1)- ####
setwd(WD_INN)

# Overall ENSANUT
ens_in00 <- read_spss("Bases/ENSA_2000/Hogar_integrantes.sav")
ens_in06 <- read_spss("Bases/ENSANUT_2006/Hogar_integrantes.sav")
ens_in12 <- read_spss("Bases/ENSANUT_2012/Hogar_integrantes.sav")
ens_in16 <- read_spss("Bases/ENSANUT_2016/Hogar_Integrantes_procesada.sav")
ens_in18 <- read_spss("Bases/ENSANUT_2018/CS_RESIDENTES.sav")
ens_in20 <- read_spss("Bases/ENSANUT_2020/integrantes_ensanut2020_w.sav")
ens_in21 <- read_spss("Bases/ENSANUT_2021/integrantes_ensanut2021_w_12_01_2022.sav")
ens_in22 <- read_spss("Bases/ENSANUT_2022/integrantes_ensanut2022_w.sav")
ens_in23 <- read_spss("Bases/ENSANUT_2023/integrantes_ensanut2023_w_n.sav")
ens_in24 <- read_spss("Bases/ENSANUT_2024/integrantes_ensanut2024_w_ICB.sav")

psu00 <- paste0(ens_in00$ent, ens_in00$mun) %>% table %>% length()
psu06 <- ens_in06$code_upm %>% table %>% length()
psu12 <- ens_in12$code_upm %>% table %>% length()
psu16 <- ens_in16$code_upm %>% table %>% length()
psu18 <- ens_in18$UPM_DIS %>% table %>% length()
psu20 <- ens_in20$Upm %>% table %>% length()
psu21 <- ens_in21$upm %>% table %>% length()
psu22 <- ens_in22$upm %>% table %>% length()
psu23 <- ens_in23$upm %>% table %>% length()
psu24 <- ens_in24$upm %>% table %>% length()

n0.00 <- nrow(ens_in00); remove(ens_in00)
n0.06 <- nrow(ens_in06); remove(ens_in06)
n0.12 <- nrow(ens_in12); remove(ens_in12)
n0.16 <- nrow(ens_in16); remove(ens_in16)
n0.18 <- nrow(ens_in18); remove(ens_in18)
n0.20 <- nrow(ens_in20); remove(ens_in20)
n0.21 <- nrow(ens_in21); remove(ens_in21)
n0.22 <- nrow(ens_in22); remove(ens_in22)
n0.23 <- nrow(ens_in23); remove(ens_in23)
n0.24 <- nrow(ens_in24); remove(ens_in24)

# Adults ≥20 years old
n1.00 <- read_spss("Bases/ENSA_2000/Adultos.sav") %>% nrow
n1.06 <- read_spss("Bases/ENSANUT_2006/Adultos.sav") %>% nrow
n1.12 <- read_spss("Bases/ENSANUT_2012/Adultos.sav") %>% nrow
n1.16 <- read_spss("Bases/ENSANUT_2016/adultos_cronicas.sav") %>% nrow
n1.18 <- read_spss("Bases/ENSANUT_2018/CS_ADULTOS.sav") %>% nrow
n1.20 <- read_spss("Bases/ENSANUT_2020/integrantes_ensanut2020_w.sav") %>%
  filter(H0303>=20 & H0303!=999) %>% nrow
n1.21 <- read_spss("Bases/ENSANUT_2021/ensadul2021_entrega_w_15_12_2021.sav") %>% nrow
n1.22 <- read_spss("Bases/ENSANUT_2022/ensadul2022_entrega_w.sav") %>% nrow
n1.23 <- read_spss("Bases/ENSANUT_2023/adultos_ensanut2023_w_n.sav") %>% nrow
n1.24 <- read_spss("Bases/ENSANUT_2024/adultos_ensanut2024_w.sav") %>% nrow

# ANALYSIS 1: Complete BP + survey weights
n2.00 <- ens00_fin_survey %>% nrow
n2.06 <- ens06_fin_survey %>% nrow
n2.12 <- ens12_fin_survey %>% nrow
n2.16 <- ens16_fin_survey %>% nrow
n2.18 <- ens18_fin_survey %>% nrow
n2.20 <- ens20_fin_survey %>% nrow
n2.21 <- ens21_fin_survey %>% nrow
n2.22 <- ens22_fin_survey %>% nrow
n2.23 <- ens23_fin_survey %>% nrow
n2.24 <- ens24_fin_survey %>% nrow

# Expands to
t.00 <- ens00_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.06 <- ens06_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.12 <- ens12_fin_survey[["variables"]][["pond.ant"]] %>% sum
ens12_fin$pond_ant %>% sum(na.rm=T) %>% format(big.mark=",")
t.16 <- ens16_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.18 <- ens18_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.20 <- ens20_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.21 <- ens21_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.22 <- ens22_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.23 <- ens23_fin_survey[["variables"]][["pond.ant"]] %>% sum
t.24 <- ens24_fin_survey[["variables"]][["pond.ant"]] %>% sum

# ANALYSIS 2: Complete covariates + total/diagnosed HTA
v1 <- c("HTN_fin","DX","Age","BMI_cat","Diabetes3","Smoking","Sex2","Area2",
  "SocSec2","Edu","DISLI_cat2","Indigenous2","Alc"); v2 <- c(v1,"HX_untreated")

na.nrow <- function(x){x %>% na.omit %>% nrow}
n3.00 <- ens00_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.06 <- ens06_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.12 <- ens12_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.16 <- ens16_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.18 <- ens18_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.20 <- ens20_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.21 <- ens21_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.22 <- ens22_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.23 <- ens23_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow
n3.24 <- ens24_fin_survey[["variables"]][v1] %>% filter(HTN_fin==1) %>% na.nrow

n4.00 <- ens00_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.06 <- ens06_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.12 <- ens12_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.16 <- ens16_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.18 <- ens18_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.20 <- ens20_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.21 <- ens21_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.22 <- ens22_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.23 <- ens23_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow
n4.24 <- ens24_fin_survey[["variables"]][v2] %>% filter(DX==1) %>% na.nrow

# Flowchart
rbind(
  paste0("n0.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma),
  paste0("n1.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma),
  paste0("n2.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma),
  paste0("t.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma),
  paste0("n3.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma),
  paste0("n4.", c("00","06","12","16","18","20","21","22","23","24")) %>%
    as.list %>% lapply(get) %>% sapply(scales::comma)) %>% as.data.frame() %>%
  `colnames<-`(c("00","06","12","16","18","20","21","22","23","24")) %>%
  `rownames<-`(c("Overall", "≥ 20 years", "Complete BP", "Expansion",
                 "UDX OR", "UTX OR"))

## Combined sample sizes ##
# Overall ENSANUT
paste0("n0.", c("00","06","12","16","18","20","21","22","23","24")) %>%
  as.list %>% sapply(get) %>% sum %>% scales::comma()
# Adults ≥20 years old
paste0("n1.", c("00","06","12","16","18","20","21","22","23","24")) %>%
  as.list %>% sapply(get) %>% sum %>% scales::comma()
# Complete BP data
paste0("n2.", c("00","06","12","16","18","20","21","22","23","24")) %>%
  as.list %>% sapply(get) %>% sum %>% scales::comma()
# Hypertension + complete data
paste0("n3.", c("00","06","12","16","18","20","21","22","23","24")) %>%
  as.list %>% sapply(get) %>% sum %>% scales::comma()
# Diagnosed hypertension + complete data
paste0("n4.", c("00","06","12","16","18","21","22","23","24")) %>%
  as.list %>% sapply(get) %>% sum %>% scales::comma()

## PSUs ##
psu00
psu06
psu12
psu16
psu18
psu20
psu21
psu22
psu23
psu24


## Strata ##
options(max.print=20000)
#2000: 4 estratos (2 por estado [rural/urbano], 4 en EDOMEX)
ens00_fin.2$estr.ant %>% table %>% length
ens00_fin.2$str_fin %>% table(ens00_fin.2$State) #Nueva variable por estado
#2006: 172 estratos (State[32] x Urb[3] x Oportunidades[2])
ens06_fin.2$estr.ant %>% table %>% length
ens06_fin.2$estr.ant %>% table(ens06_fin.2$State)
#2012: 155 estratos (State[32] x ]Urb[4] x SLI[2])
ens12_fin.2$estr.ant %>% table %>% length
ens12_fin.2$estr.ant %>% table(ens12_fin.2$State)
#2016: 39 estratos (Reg[4] x Urb[2] x SLI[3] + estados individuales)
ens16_fin.2$estr.ant %>% table %>% length
#2018: 325 estratos (State[32] x Urb[3] x SES[4])
ens18_fin.2$estr.ant %>% table %>% length
ens18_fin.2$estr.ant %>% table(ens18_fin.2$State)
#2021-2024 (State[32] x Urb[3])
ens21_fin.2$estr.ant %>% table %>% length #91
ens22_fin.2$estr.ant %>% table %>% length #94
ens23_fin.2$estr.ant %>% table %>% length #61
ens24_fin.2$estr.ant %>% table %>% length #91


####----------------------- Hypertension prevalence --------------------#### ----####
#### Custom functions------------- ####
get.prev <- function(data, by){
  #1:Year // 2:Sex // 3:Age // 4:BMI // 5:T2D 
  svyby(formula=out, by=strats[[by]], design=data,
        svymean, na.rm=T) %>% prev.mutate(by)}

prev.mutate <- function(x, by){
  if(by==1){
    x %>% `colnames<-`(c("Year", "out", "se")) %>% mutate(
      "prop"=round(100*out, 1), "IC95"=round(100*(se*qnorm(0.975)), 2),
      "lIC95"=prop-IC95, "uIC95"=prop+IC95, "cluster"=label) %>%
      select(prop, IC95, lIC95, uIC95, cluster, Year)}
  else{
    x %>% `colnames<-`(c("Year", "group", "out", "se")) %>% mutate(
      "prop"=round(100*out, 1), "IC95"=round(100*(se*qnorm(0.975)), 2),
      "lIC95"=prop-IC95, "uIC95"=prop+IC95, "cluster"=label) %>%
      select(prop, IC95, lIC95, uIC95, cluster, group, Year)}}

strats <- as.list(paste0("~Year", c(
  "","+Sex2","+Age_cat","+BMI_cat","+Diabetes2"))) %>%
  lapply(as.formula, env=parent.frame())

# Function: direct age-sex standardization across multiple survey cycles
std_prev <- function(svy_list, 
                     outcome, 
                     age_var, 
                     sex_var, 
                     cycle_names = NULL,
                     ref_dist) {
  library(survey)
  library(dplyr)
  options(survey.adjust.domain.lonely = TRUE, survey.lonely.psu="adjust")
  # svy_list: list of survey design objects (one per cycle)
  # outcome: character, binary variable name (e.g., "hypertension")
  # age_var: character, categorical age group variable name (e.g., "age_cat")
  # sex_var: character, sex variable name (e.g., "sex")
  # cycle_names: optional vector of names (length = length(svy_list))
  # ref_dist: data.frame with columns: strata_as, prop_ref (sum = 1)
  
  # Ensure cycle names
  if (is.null(cycle_names)) {
    cycle_names <- paste0("Cycle_", seq_along(svy_list))
  }
  
  results <- list()
  
  for (i in seq_along(svy_list)) {
    svy <- svy_list[[i]]
    cycle <- cycle_names[i]
    
    # Create stratum variable (age-sex)
    svy <- update(svy, strata_as = interaction(
      get(age_var), get(sex_var), sep = "_"))
    
    # Stratum-specific prevalence
    stratum_prev <- svyby(as.formula(paste0("~", outcome)),
                          ~strata_as,
                          svy,
                          svymean,
                          vartype = "se",
                          na.rm = TRUE)
    
    # Merge with reference distribution
    merged <- left_join(stratum_prev, ref_dist, by = c("strata_as"))
    
    # Standardized prevalence and SE (delta method)
    std_prev <- 100*sum(merged[[outcome]] * merged$prop_ref, na.rm = TRUE)
    var_std  <- sum((merged$se^2) * (merged$prop_ref^2), na.rm = TRUE)
    se_std   <- 100*sqrt(var_std)
    
    results[[cycle]] <- data.frame(
      cycle = cycle,
      std_prev = std_prev,
      se_std = se_std,
      lower = std_prev - 1.96 * se_std,
      upper = std_prev + 1.96 * se_std
    )
  }
  
  bind_rows(results)
}


#### Age and sex standardization-- ####

#Load CONAPO dataset
setwd(WD_OUT)
conapo <-read_csv("Hypertension phenotypes/pobproy_quinq1.csv")

ref_df<- conapo %>%
  group_by(ANO, SEXO) %>% summarise(
    #Total population by age categories
    "POB_20_39"=sum(POB_20_24)+sum(POB_25_29)+sum(POB_30_34)+sum(POB_35_39),
    "POB_40_59"=sum(POB_40_44)+sum(POB_45_49)+sum(POB_50_54)+sum(POB_55_59),
    "POB_60"=sum(POB_60_64)+sum(POB_65_69)+sum(POB_70_74)+sum(POB_75_79)+
      sum(POB_80_84)+sum(POB_85_mm), .groups = "drop") %>%
  #Rearrange dataset
  pivot_longer(
    cols = starts_with("POB_"), names_to = "AGE_GROUP",
    values_to = "population_total") %>%
  rename("Age_cat"=AGE_GROUP, "Sex2"=SEXO) %>%
  #Create new variables
  mutate(
    "Age_cat"=factor(Age_cat, labels=c("20-39", "40-59", "60+"), ordered=T),
    "Sex2"=factor(Sex2, labels =c("Men", "Women"), ordered = T)) %>%
  filter(ANO %in% c(2018)) %>%
  mutate(
    "strata_as" = levels(interaction(
      c("20-39","40-59","60+"), c("Men","Women"), sep = "_")),
    "prop_ref"=population_total/sum(population_total)) %>%
  select(strata_as, prop_ref)

ref_df$strata_as[ref_df$strata_as=="60+_Men"] <- "≥60_Men"
ref_df$strata_as[ref_df$strata_as=="60+_Women"] <- "≥60_Women"

svy_list <- list(ens00_fin_survey, ens06_fin_survey, ens12_fin_survey,
                 ens16_fin_survey, ens18_fin_survey, ens20_fin_survey,
                 ens21_fin_survey, ens22_fin_survey, ens23_fin_survey,
                 ens24_fin_survey)

cycle_names <- as.character(2000 + c(0, 6, 12, 16, 18, 20:24))

std_HTN_fin <- std_prev(
  svy_list, outcome = "HTN_fin",
                        age_var = "Age_cat",
                        sex_var = "Sex2",
                        cycle_names = cycle_names,
                        ref_dist = ref_df) %>%
  mutate(cluster="Total Hypertension Adjusted")

std_NEW_HTN <- std_prev(svy_list,
                        outcome = "NEW_HTN",
                        age_var = "Age_cat",
                        sex_var = "Sex2",
                        cycle_names = cycle_names,
                        ref_dist = ref_df)%>%
  mutate(cluster="Undiagnosed Hypertension Adjusted")

std_DX <- std_prev(svy_list,
                        outcome = "DX",
                        age_var = "Age_cat",
                        sex_var = "Sex2",
                        cycle_names = cycle_names,
                        ref_dist = ref_df)%>%
  mutate(cluster="Diagnosed Hypertension Adjusted")

std_ISH <- std_prev(svy_list,
                   outcome = "ISH",
                   age_var = "Age_cat",
                   sex_var = "Sex2",
                   cycle_names = cycle_names,
                   ref_dist = ref_df)%>%
  mutate(cluster="Isolated Systolic Hypertension Adjusted")

std_IDH <- std_prev(svy_list,
                    outcome = "IDH",
                    age_var = "Age_cat",
                    sex_var = "Sex2",
                    cycle_names = cycle_names,
                    ref_dist = ref_df)%>%
  mutate(cluster="Isolated Diastolic Hypertension Adjusted")

std_SDH <- std_prev(svy_list,
                    outcome = "SDH",
                    age_var = "Age_cat",
                    sex_var = "Sex2",
                    cycle_names = cycle_names,
                    ref_dist = ref_df)%>%
  mutate(cluster="Systodiastolic Hypertension Adjusted")

svy_list2 <- list(
  ens00_fin_survey4, ens06_fin_survey4, ens12_fin_survey4,
  ens16_fin_survey4, ens18_fin_survey4, ens21_fin_survey4,
  ens22_fin_survey4, ens23_fin_survey4, ens24_fin_survey4)

cycle_names2 <- as.character(2000 + c(0, 6, 12, 16, 18, 21:24))

std_Treated<- std_prev(svy_list2,
                    outcome = "HX_treated",
                    age_var = "Age_cat",
                    sex_var = "Sex2",
                    cycle_names = cycle_names2,
                    ref_dist = ref_df) %>%
  mutate(cluster="Treated Hypertension Adjusted")

std_Untreated<- std_prev(svy_list2,
                       outcome = "HX_untreated",
                       age_var = "Age_cat",
                       sex_var = "Sex2",
                       cycle_names = cycle_names2,
                       ref_dist = ref_df)%>%
  mutate(cluster="Untreated Hypertension Adjusted")

std_Controlled<- std_prev(svy_list2,
                         outcome = "HTN_goal",
                         age_var = "Age_cat",
                         sex_var = "Sex2",
                         cycle_names = cycle_names2,
                         ref_dist = ref_df)%>%
  mutate(cluster="Controlled Hypertension Adjusted")

std_Uncontrolled<- std_prev(svy_list2,
                         outcome = "HTN_not_goal",
                         age_var = "Age_cat",
                         sex_var = "Sex2",
                         cycle_names = cycle_names2,
                         ref_dist = ref_df)%>%
  mutate(cluster="Uncontrolled Hypertension Adjusted")


#### Total -------------(HTN_fin)- ####
### 1) Prevalence by Year (ESC/ESH)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin"); st <- 1; HTN_fin_prev_year <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin.A"); st <- 1; HTN_fin_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin"); st <- 2; HTN_fin_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin"); st <- 3; HTN_fin_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin"); st <- 4; HTN_fin_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Total Hypertension"; out <- as.formula(
  "~HTN_fin"); st <- 5; HTN_fin_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)



#### Undiagnosed -------(NEW_HTN)- ####
### 1) Prevalence by Year (ESC/ESH)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 1; UHT_prev_year <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN.A"); st <- 1; UHT_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 2; UHT_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 3; UHT_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 4; UHT_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 5; UHT_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### Diagnosed --------------(DX)- ####
### 1) Prevalence by Year (ESC/ESH)
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 1; DX_prev_year <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 2; DX_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 3; DX_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 4; DX_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 5; DX_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### High-Normal BP --------(PRE)- ####
### 1) Prevalence by Year (ESC/ESH)
label <- "High-Normal BP"; out <- as.formula(
  "~PRE"); st <- 1; PRE_prev_year <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Elevated BP"; out <- as.formula(
  "~PRE.A"); st <- 1; PRE_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "High-Normal BP"; out <- as.formula(
  "~PRE"); st <- 2; PRE_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "High-Normal BP"; out <- as.formula(
  "~PRE"); st <- 3; PRE_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "High-Normal BP"; out <- as.formula(
  "~PRE"); st <- 4; PRE_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "High-Normal BP"; out <- as.formula(
  "~PRE"); st <- 5; PRE_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### Isolated Systolic -----(ISH)- ####
### 0) Prevalence within UDH (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 1; ISH_prev_year <- rbind(
    get.prev(ens00_fin_survey3, st), get.prev(ens06_fin_survey3, st),
    get.prev(ens12_fin_survey3, st), get.prev(ens16_fin_survey3, st),
    get.prev(ens18_fin_survey3, st), get.prev(ens20_fin_survey3, st),
    get.prev(ens21_fin_survey3, st), get.prev(ens22_fin_survey3, st),
    get.prev(ens23_fin_survey3, st), get.prev(ens24_fin_survey3, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 1; ISH_prev_year_overall <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH.A"); st <- 1; ISH_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 2; ISH_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 3; ISH_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 4; ISH_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Isolated Systolic Hypertension"; out <- as.formula(
  "~ISH"); st <- 5; ISH_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### Isolated Diastolic ----(IDH)- ####
### 0) Prevalence within UDH (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 1; IDH_prev_year <- rbind(
    get.prev(ens00_fin_survey3, st), get.prev(ens06_fin_survey3, st),
    get.prev(ens12_fin_survey3, st), get.prev(ens16_fin_survey3, st),
    get.prev(ens18_fin_survey3, st), get.prev(ens20_fin_survey3, st),
    get.prev(ens21_fin_survey3, st), get.prev(ens22_fin_survey3, st),
    get.prev(ens23_fin_survey3, st), get.prev(ens24_fin_survey3, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 1; IDH_prev_year_overall <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH.A"); st <- 1; IDH_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 2; IDH_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 3; IDH_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 4; IDH_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Isolated Diastolic Hypertension"; out <- as.formula(
  "~IDH"); st <- 5; IDH_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### Systodiastolic --------(SDH)- ####
### 0) Prevalence within UDH (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 1; SDH_prev_year <- rbind(
    get.prev(ens00_fin_survey3, st), get.prev(ens06_fin_survey3, st),
    get.prev(ens12_fin_survey3, st), get.prev(ens16_fin_survey3, st),
    get.prev(ens18_fin_survey3, st), get.prev(ens20_fin_survey3, st),
    get.prev(ens21_fin_survey3, st), get.prev(ens22_fin_survey3, st),
    get.prev(ens23_fin_survey3, st), get.prev(ens24_fin_survey3, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 1; SDH_prev_year_overall <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>% 
  as.data.frame %>% `rownames<-`(NULL)

### 1) Prevalence by Year (ACC/AHA)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH.A"); st <- 1; SDH_prev_year.A <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>% 
  as.data.frame %>% `rownames<-`(NULL)

### 2) Prevalence by Sex (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 2; SDH_prev_sex <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### 3) Prevalence by Age (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 3; SDH_prev_age <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>% 
  as.data.frame %>% `rownames<-`(NULL)

### 4) Prevalence by BMI (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 4; SDH_prev_bmi <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>% 
  as.data.frame %>% `rownames<-`(NULL)

### 5) Prevalence by T2D (ESC/ESH)
label <- "Systodiastolic Hypertension"; out <- as.formula(
  "~SDH"); st <- 5; SDH_prev_t2d <- rbind(
    get.prev(ens00_fin_survey, st), get.prev(ens06_fin_survey, st),
    get.prev(ens12_fin_survey, st), get.prev(ens16_fin_survey, st),
    get.prev(ens18_fin_survey, st), get.prev(ens20_fin_survey, st),
    get.prev(ens21_fin_survey, st), get.prev(ens22_fin_survey, st),
    get.prev(ens23_fin_survey, st), get.prev(ens24_fin_survey, st)) %>% 
  as.data.frame %>% `rownames<-`(NULL)


#### Treatment and control in DX-- ####
### Undiagnosed Hypertension Proportion
label <- "Undiagnosed Hypertension"; out <- as.formula(
  "~NEW_HTN"); st <- 1; UTH_prop_year <- rbind(
    get.prev(ens00_fin_survey2, st), get.prev(ens06_fin_survey2, st),
    get.prev(ens12_fin_survey2, st), get.prev(ens16_fin_survey2, st),
    get.prev(ens18_fin_survey2, st), get.prev(ens20_fin_survey2, st),
    get.prev(ens21_fin_survey2, st), get.prev(ens22_fin_survey2, st),
    get.prev(ens23_fin_survey2, st), get.prev(ens24_fin_survey2, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### Diagnosed Hypertension Proportion
label <- "Diagnosed Hypertension"; out <- as.formula(
  "~DX"); st <- 1; DX_prop_year <- rbind(
    get.prev(ens00_fin_survey2, st), get.prev(ens06_fin_survey2, st),
    get.prev(ens12_fin_survey2, st), get.prev(ens16_fin_survey2, st),
    get.prev(ens18_fin_survey2, st), get.prev(ens20_fin_survey2, st),
    get.prev(ens21_fin_survey2, st), get.prev(ens22_fin_survey2, st),
    get.prev(ens23_fin_survey2, st), get.prev(ens24_fin_survey2, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### Treated Hypertension Proportion
label <- "Treated Hypertension"; out <- as.formula(
  "~HX_treated"); st <- 1; HX_treated_prev_year <- rbind(
    get.prev(ens00_fin_survey4, st), get.prev(ens06_fin_survey4, st),
    get.prev(ens12_fin_survey4, st), get.prev(ens16_fin_survey4, st),
    get.prev(ens18_fin_survey4, st), #Data not available for 2020
    get.prev(ens21_fin_survey4, st), get.prev(ens22_fin_survey4, st),
    get.prev(ens23_fin_survey4, st), get.prev(ens24_fin_survey4, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### Untreated Hypertension Proportion
label <- "Untreated Hypertension"; out <- as.formula(
  "~HX_untreated"); st <- 1; HX_untreated_prev_year <- rbind(
    get.prev(ens00_fin_survey4, st), get.prev(ens06_fin_survey4, st),
    get.prev(ens12_fin_survey4, st), get.prev(ens16_fin_survey4, st),
    get.prev(ens18_fin_survey4, st), #Data not available for 2020
    get.prev(ens21_fin_survey4, st), get.prev(ens22_fin_survey4, st),
    get.prev(ens23_fin_survey4, st), get.prev(ens24_fin_survey4, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### Goal Achievement Proportion
label <- "Controlled"; out <- as.formula(
  "~HTN_goal"); st <- 1; HX_goal_prev_year <- rbind(
    get.prev(ens00_fin_survey4, st), get.prev(ens06_fin_survey4, st),
    get.prev(ens12_fin_survey4, st), get.prev(ens16_fin_survey4, st),
    get.prev(ens18_fin_survey4, st), #Data not available for 2020
    get.prev(ens21_fin_survey4, st), get.prev(ens22_fin_survey4, st),
    get.prev(ens23_fin_survey4, st), get.prev(ens24_fin_survey4, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)

### Goal Non-Achievement Proportion
label <- "Uncontrolled"; out <- as.formula(
  "~HTN_not_goal"); st <- 1; HX_not_goal_prev_year <- rbind(
    get.prev(ens00_fin_survey4, st), get.prev(ens06_fin_survey4, st),
    get.prev(ens12_fin_survey4, st), get.prev(ens16_fin_survey4, st),
    get.prev(ens18_fin_survey4, st), #Data not available for 2020
    get.prev(ens21_fin_survey4, st), get.prev(ens22_fin_survey4, st),
    get.prev(ens23_fin_survey4, st), get.prev(ens24_fin_survey4, st)) %>%
  as.data.frame %>% `rownames<-`(NULL)


#### Total estimates for 2024----- ####

### Total adults with undiagnosed hypertension
svyby(~HTN_fin, by=~Year, design=ens24_fin_survey, svytotal, na.rm=T) %>%  
  mutate(prop=round(HTN_fin, digits=1), IC95=round((se*1.96), digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Undiagnosed") %>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.HTN_fin

### Total adults with undiagnosed hypertension
svyby(~NEW_HTN, by=~Year, design=ens24_fin_survey, svytotal, na.rm=T) %>%  
  mutate(prop=round(NEW_HTN, digits=1), IC95=round((se*1.96), digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Undiagnosed") %>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.NEW_HTN

### Total adults with diagnosed hypertension
svyby(~DX, by=~Year, design=ens24_fin_survey, svytotal, na.rm=T) %>%  
  mutate(prop=round(DX, digits=1), IC95=round((se*1.96), digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Diagnosed") %>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.DX

### Total adults without treatment
svyby(formula=~HX_untreated, by=~Year,
      design=ens24_fin_survey4, svytotal, na.rm=T) %>%  
  mutate(prop=round(HX_untreated, 1), IC95=round((se*1.96), 2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Untreated") %>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.HX_untreated

### Total adults not reaching goals
svyby(formula=~HTN_not_goal, by=~Year,
      design=ens24_fin_survey4, svytotal, na.rm=T) %>%  
  mutate(prop=round(HTN_not_goal, 1), IC95=round((se*1.96), 2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Uncontrolled") %>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.HTN_not_goal

### Total adults untreated and not reaching goals
svyby(formula=~HTN_not_goal_untx, by=~Year,
      design=ens24_fin_survey4, svytotal, na.rm=T) %>%  
  mutate(prop=round(HTN_not_goal_untx, 1), IC95=round((se*1.96), 2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Untreated, uncontrolled")%>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.HTN_not_goal_untx

### Total adults untreated and reaching goals
svyby(formula=~HTN_goal_untx, by=~Year,
      design=ens24_fin_survey4, svytotal, na.rm=T) %>%  
  mutate(prop=round(HTN_goal_untx, 1), IC95=round((se*1.96), 2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Untreated, uncontrolled")%>%
  select(prop,IC95,lIC95,uIC95,cluster,Year) -> total.HTN_goal_untx



DX_prop_year[10,c(1,3:4)] %>% round(1) #Percentage aware of their condition
total.DX[c(1,3:4)] %>% format(big.mark=",") ### Total with diagnosed hypertension

UTH_prop_year[10,c(1,3:4)] %>% round(1) #Percentage of undiagnosed hypertension
total.NEW_HTN[c(1,3:4)] %>% format(big.mark=",") ### Undiagnosed hypertension

total.HX_untreated[c(1,3:4)] %>% format(big.mark=",") ## Without treatment
total.HTN_not_goal[c(1,3:4)] %>% format(big.mark=",") ## Not reaching goals
total.HTN_not_goal_untx[c(1,3:4)] %>% floor %>% format(big.mark=",") ## UTX+UCT

total.HX_untreated[c(1)]/total.DX[c(1)] #% diagnosed w/o tx
total.HTN_not_goal[c(1)]/total.DX[c(1)] #% diagnosed w/o control

total.HTN_fin[1] %>% format(big.mark=",") #Total living with hypertension


####----------------------- Trends (Poisson) + predictors (logistic) ---#### ----####
#### Pois: Pooled data + design--- ####
## Variables that will be used in Poisson models
old_n <- c("PRE", "PRE.A", "HTN_fin", "HTN_fin.A", "NEW_HTN", "NEW_HTN.A",
           "DX", "ISH", "IDH", "SDH", "Age", "Sex", "Year",
           "swts", "psu_fin", "str_fin", "wgt_fin")

## New (simpler) names for those variables
new_n <- c("PRE", "PRE.A", "HTN_fin", "HTN_fin.A", "NEW_HTN", "NEW_HTN.A",
           "DX", "ISH", "IDH", "SDH", "age", "sex", "year",
           "swts", "psu", "str", "wgt")

## List of surveys to apply loops more easily
surv_list <- list(
  ens00_fin.2, ens06_fin.2, ens12_fin.2, ens18_fin.2, ens16_fin.2,
  ens20_fin.2, ens21_fin.2, ens22_fin.2, ens23_fin.2, ens24_fin.2)

## For each survey, select the desired variables and stack in a single data set
pool.poisson0 <- list(); for(i in seq_along(surv_list)){
  pool <- surv_list[[i]] %>% select(all_of(old_n))
  pool.poisson0 <- append(pool.poisson0, list(pool))
  }; pool.poisson0 <- do.call(rbind, pool.poisson0) # Pooled data set

## Take the pooled dataset and rename each variable with the new names
for(i in seq_along(old_n)){pool.poisson0[new_n[i]] <- pool.poisson0[old_n[i]]}
pool.poisson <- pool.poisson0 %>% select(all_of(new_n)) %>%
  mutate("yrs" = year - 2000, "wgt2"=wgt/10, "eff" = paste0(year, str))
pool.poisson.c <- filter(pool.poisson, yrs>=16) #Concise data set

## Pooled survey design
#Main survey
#survpoiss.f1 <- svydesign(data=pool.poisson, ids=~psu, weights=~wgt2,
#                          strata=~interaction(year,str), nest=T)
#Concise models (2016-2024)
#survpoiss.c1 <- subset(survpoiss.f1, yrs>=16)
#
#Model example
#svyglm(formula = HTN_fin ~ yrs + scale(age) + sex,
#       design = survpoiss.f1, family = quasipoisson())



#### Pois: Mixed effects models--- ####

##-------------------------- 01 PRE ------------------------------------------##
#Full model (2000-2024)
poiss.01f <- glmer(PRE~yrs+scale(age)+sex+(1|eff), family = poisson(),
                   data = filter(pool.poisson, HTN_fin==0), weights = swts)
#Concise model (2016-2024)
poiss.01c <- glmer(PRE~yrs+scale(age)+sex+(1|eff), family = poisson(),
                   data = filter(pool.poisson.c, HTN_fin==0), weights = swts)


##-------------------------- 02 PRE.A ----------------------------------------##
#Full model (2000-2024)
poiss.02f <- glmer(PRE.A~yrs+scale(age)+sex+(1|eff), family = poisson(),
                   data = filter(pool.poisson, HTN_fin.A==0), weights = swts)
#Concise model (2016-2024)
poiss.02c <- glmer(PRE.A~yrs+scale(age)+sex+(1|eff), family = poisson(),
                   data = filter(pool.poisson.c, HTN_fin.A==0), weights = swts)


##-------------------------- 03 HTN_fin --------------------------------------##
#Full model (2000-2024)
poiss.03f <- glmer(HTN_fin~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.03c <- glmer(HTN_fin~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 04 HTN_fin.A ------------------------------------##
#Full model (2000-2024)
poiss.04f <- glmer(HTN_fin.A~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.04c <- glmer(HTN_fin.A~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 05 NEW_HTN --------------------------------------##
#Full model (2000-2024)
poiss.05f <- glmer(NEW_HTN~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.05c <- glmer(NEW_HTN~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 06 NEW_HTN.A ------------------------------------##
#Full model (2000-2024)
poiss.06f <- glmer(NEW_HTN.A~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.06c <- glmer(NEW_HTN.A~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 07 DX -------------------------------------------##
#Full model (2000-2024)
poiss.07f <- glmer(DX~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.07c <- glmer(DX~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 08 ISH ------------------------------------------##
#Full model (2000-2024)
poiss.08f <- glmer(ISH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.08c <- glmer(ISH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 09 IDH ------------------------------------------##
#Full model (2000-2024)
poiss.09f <- glmer(IDH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.09c <- glmer(IDH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)


##-------------------------- 10 SDH ------------------------------------------##
#Full model (2000-2024)
poiss.10f <- glmer(SDH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson, weights = swts)
#Concise model (2016-2024)
poiss.10c <- glmer(SDH~yrs+scale(age)+sex+(1|eff),
                   family = poisson(), data = pool.poisson.c, weights = swts)



#### Pois: PRs + overdispersion--- ####

# Function: summarize PR for the 1st non-intercept coefficient in a glmer
PR_fin <- function(x){
  a <- summary(x)$coef[2,1:2]
  b <- c(a[1], a[1] + c(-1,1)*a[2]*qnorm(.975))
  c <- sprintf(exp(b), fmt="%#.3f")
  d <- paste0(c[1], " (", c[2], "-", c[3], ")")
  return(d)}


# PRs ready to insert into manuscript text
#ESC/ESH #2000-2024
PR_fin(poiss.03f) #HTN_fin
PR_fin(poiss.05f) #NEW_HTN
PR_fin(poiss.07f) #DX
PR_fin(poiss.01f) #PRE

#ESC/ESH #2016-2024
PR_fin(poiss.03c) #HTN_fin
PR_fin(poiss.05c) #NEW_HTN
PR_fin(poiss.07c) #DX
PR_fin(poiss.01c) #PRE

#ACC/AHA #2000-2024
PR_fin(poiss.04f) #HTN_fin.A
PR_fin(poiss.06f) #NEW_HTN.A
PR_fin(poiss.07f) #DX
PR_fin(poiss.02f) #PRE.A

#ACC/AHA #2016-2024
PR_fin(poiss.04c) #HTN_fin.A
PR_fin(poiss.06c) #NEW_HTN.A
PR_fin(poiss.07c) #DX
PR_fin(poiss.02c) #PRE.A

#UDH Phenotypes
PR_fin(poiss.08f) #ISH 2000-2024
PR_fin(poiss.09f) #IDH 2000-2024
PR_fin(poiss.10f) #SDH 2000-2024

PR_fin(poiss.08c) #ISH 2016-2024
PR_fin(poiss.09c) #IDH 2016-2024
PR_fin(poiss.10c) #SDH 2016-2024


# Check all Poisson models for overdispersion
check_overdispersion(poiss.01f) #PRE
check_overdispersion(poiss.02f) #PRE.A
check_overdispersion(poiss.03f) #HTN_fin
check_overdispersion(poiss.04f) #HTN_fin.A
check_overdispersion(poiss.05f) #NEW_HTN
check_overdispersion(poiss.06f) #NEW_HTN.A
check_overdispersion(poiss.07f) #DX
check_overdispersion(poiss.08f) #ISH
check_overdispersion(poiss.09f) #IDH
check_overdispersion(poiss.10f) #SDH

#Note: Consider using COM-Poisson models to deal with underdispersion
#pacman::p_load(DHARMa)
#testDispersion(simulateResiduals(fittedModel = poiss.01f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.02f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.03f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.04f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.05f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.06f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.07f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.08f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.09f, plot = F))
#testDispersion(simulateResiduals(fittedModel = poiss.10f, plot = F))



#### Pois: Seq. adjusted trends--- ####

## Variables that will be used in Poisson models
old_n_ <- c("HTN_fin", "NEW_HTN", "DX", "PRE",
           "Age", "Sex", "Year", "Alc", "Smoking",
           "Edu", "DISLI_cat2", "Area2", "SocSec2", "Indigenous2",
           "BMI_cat", "Diabetes3",
           "swts", "psu_fin", "str_fin", "wgt_fin")

## New (simpler) names for those variables
new_n_ <- c("HTN_fin", "NEW_HTN", "DX", "PRE",
           "age", "sex", "year",
           "alc", "smo", "edu", "sli", "urb", "sec", "ind", "obe", "dm2",
           "swts", "psu", "str", "wgt")

## List of surveys to apply loops more easily
surv_list_ <- list(
  ens00_fin.2, ens06_fin.2, ens12_fin.2, ens18_fin.2, ens16_fin.2,
  ens20_fin.2, ens21_fin.2, ens22_fin.2, ens23_fin.2, ens24_fin.2)

## For each survey, select the desired variables and stack in a single data set
pool.poisson0_ <- list(); for(i in seq_along(surv_list_)){
  pool <- surv_list_[[i]] %>% select(all_of(old_n_))
  pool.poisson0_ <- append(pool.poisson0_, list(pool))
}; pool.poisson0_ <- do.call(rbind, pool.poisson0_) # Pooled data set

## Take the pooled dataset and rename each variable with the new names
for(i in seq_along(old_n_)){
  pool.poisson0_[new_n_[i]] <- pool.poisson0_[old_n_[i]]}
pool.poisson_ <- pool.poisson0_ %>% select(all_of(new_n_)) %>%
  mutate("yrs" = year - 2000, "wgt2"=wgt/10, "eff" = paste0(year, str))
pool.poisson.c_ <- filter(pool.poisson_, yrs>=16) 



#Fully adjusted
#Overall
poiss.adj.1 <- glm(
  HTN_fin~yrs+scale(age)+sex+edu+sli+urb+sec+ind+smo+alc+dm2+obe,
  family = poisson(), data = pool.poisson_, weights = swts)
#New
poiss.adj.2 <- glm(
  NEW_HTN~yrs+scale(age)+sex+edu+sli+urb+sec+ind+smo+alc+dm2+obe,
  family = poisson(), data = pool.poisson_, weights = swts)
#DX
poiss.adj.3 <- glm(
  DX~yrs+scale(age)+sex+edu+sli+urb+sec+ind+smo+alc+dm2+obe,
  family = poisson(), data = pool.poisson_, weights = swts)
#PRE
poiss.adj.4 <- glm(
  PRE~yrs+scale(age)+sex+edu+sli+urb+sec+ind+smo+alc+dm2+obe,
  family = poisson(), data = filter(pool.poisson_, HTN_fin==0), weights=swts)


PR_fin(poiss.adj.1)
PR_fin(poiss.adj.2)
PR_fin(poiss.adj.3)
PR_fin(poiss.adj.4)



#### Logistic: Pooled data-------- ####
## Variables that will be used in logistic models
old_n <- c("HTN_fin", "NEW_HTN", "DX", "HX_untreated", "HTN_not_goal",
           "Year", "Age", "Sex2", "Area2", "BMI_cat", "Diabetes3",
           "Alc", "Smoking", "SocSec2", "Edu", "DISLI_cat2",
           "Indigenous2", "swts", "psu_fin", "str_fin", "wgt_fin")

## New (simpler) names for those variables
new_n <- c("HTN_fin", "NEW_HTN", "DX", "HX_untreated", "HTN_not_goal",
           "year", "age", "sex", "urb", "obe", "dm2",
           "alc", "tab", "afi", "edu", "sli",
           "ind", "swts", "psu", "str", "wgt")

## List of surveys to apply loops more easily
ens20_fin.2$HTN_not_goal <- ifelse(ens20_fin.2$HTN_goal==1, 0, 1)
surv_list <- list(
  ens00_fin.2, ens06_fin.2, ens12_fin.2, ens18_fin.2, ens16_fin.2,
  ens20_fin.2, ens21_fin.2, ens22_fin.2, ens23_fin.2, ens24_fin.2)

## For each survey, select the desired variables and stack in a single data set
pool.binomial0 <- list(); for(i in seq_along(surv_list)){
  pool <- surv_list[[i]] %>% select(all_of(old_n))
  pool.binomial0 <- append(pool.binomial0, list(pool))
}; pool.binomial0 <- do.call(rbind, pool.binomial0) # Pooled data set


## Take the pooled dataset and rename each variable with the new names
for(i in seq_along(old_n)){pool.binomial0[new_n[i]] <- pool.binomial0[old_n[i]]}
pool.binomial <- pool.binomial0 %>% select(all_of(new_n)) %>%
  
  ## Create new variables, set predictors as factors and adjust levels
  mutate(
    "yrs"=year - 2000, "wgt2"=wgt/10, "eff"=paste0(year, str), "age"=age/10,
    "sex"=factor(sex, ordered = F, levels=c("Women", "Men")),
    "tab"=factor(tab, ordered = F, levels=c("Never", "Former", "Current")),
    "alc"=factor(alc, ordered = F, levels=c(
      "Not currently", "Less than daily", "Daily")),
    "edu"=factor(edu, ordered = F, levels=c(
      "None", "Primary", "Secondary", "Tertiary"), labels=c(
        "None/Primary", "None/Primary", "Secondary", "Tertiary")),
    "ind"=factor(ind, ordered = F, levels=c("Non-indigenous", "Indigenous")),
    "afi"=factor(afi, ordered = F, levels=c("Affiliated", "Not affiliated")),
    "sli"=factor(sli, levels = 0:1, labels=c("Low/Medium", "High")),
    "urb"=factor(urb, ordered = F, levels=c("Urban", "Rural")),
    "obe"=factor(obe, ordered = F, levels=c("No obesity", "Obesity")),
    "dm2"=factor(dm2, ordered = F, levels=c("Non-diagosed", "Diagnosed"),
                 labels=c("No diabetes", "Diabetes")))

## Set labels for predictors
lr.names <- c(
  "Age", "Sex", "Smoking", "Alcohol\nintake", "Edu.\nlevel", "Ind.\nidentity",
  "Social\nsecurity", "DISLI", "Area", "BMI-\nobesity", "Prior\nTD2", "Years")

attr(pool.binomial$age, "label") <- lr.names[1]
attr(pool.binomial$sex, "label") <- lr.names[2]
attr(pool.binomial$tab, "label") <- lr.names[3]
attr(pool.binomial$alc, "label") <- lr.names[4]
attr(pool.binomial$edu, "label") <- lr.names[5]
attr(pool.binomial$ind, "label") <- lr.names[6]
attr(pool.binomial$afi, "label") <- lr.names[7]
attr(pool.binomial$sli, "label") <- lr.names[8]
attr(pool.binomial$urb, "label") <- lr.names[9]
attr(pool.binomial$obe, "label") <- lr.names[10]
attr(pool.binomial$dm2, "label") <- lr.names[11]
attr(pool.binomial$yrs, "label") <- lr.names[12]



#### Logistic: Models------------- ####

##--------------------- BOBYQA controller ----------------------##
#Bound Optimization BY Quadratic Approximation
#Stable random effect estimate
#Stable fixed effects estimates
#Computation time: 3.023 min
#start2 <- Sys.time()
##-------------------------- 01 NEW_HTN --------------------------------------##
pool.binomial %>% select(-HX_untreated) %>% na.omit %>% filter(HTN_fin==1) %>%
  glmer(formula = NEW_HTN~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
          (1|year), family = binomial(), weights = swts,
        control = glmerControl(optimizer = "bobyqa")) -> binom.01
##-------------------------- 02 HX_untreated ---------------------------------##
pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) %>%
  glmer(formula = HX_untreated~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
          (1|year), family = binomial(), weights = swts,
        control = glmerControl(optimizer = "bobyqa")) -> binom.02
##-------------------------- 03 HTN_not_goal ---------------------------------##
pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) %>%
  glmer(formula = HTN_not_goal~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
          (1|year), family = binomial(), weights = swts,
        control = glmerControl(optimizer = "bobyqa")) -> binom.03
#end2 <- Sys.time()


##--------------------- Standard Nelder-Mead controller ----------------------##
#Fails to converge for random effect
#Stable fixed effects estimates
#Computation time: 3.348 min
#start1 <- Sys.time()
##-------------------------- 01 HTN_fin --------------------------------------##
#pool.binomial %>% select(-HX_untreated) %>% na.omit %>% filter(HTN_fin==1) %>%
#  glmer(formula = NEW_HTN~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
#          (1|year), family = binomial(), weights = swts) -> binom.01
##-------------------------- 02 HX_untreated ---------------------------------##
#pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) %>%
#  glmer(formula = HX_untreated~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
#          (1|year), family = binomial(), weights = swts) -> binom.02
##-------------------------- 03 HTN_not_goal ---------------------------------##
#pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) %>%
#  glmer(formula = HTN_not_goal~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs+
#          (1|year), family = binomial(), weights = swts) -> binom.03
#end1 <- Sys.time()


#### Logistic: GOF---------------- ####
pool.binomial %>% select(-HX_untreated) %>% na.omit %>%
  filter(HTN_fin==1) -> pooled.DX
pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) -> pooled.TX
pool.binomial %>% na.omit %>% filter(DX==1, year!=2020) -> pooled.GOAL

#Undiagnosed
glm(formula = out~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
      family=binomial(), weights=swts, data=pooled.DX) -> m_test1
glmtoolbox::hltest(m_test1)
#(ns, p=0.36206)

pooled.DX$age_cat <- cut(pooled.DX$age*10, c(-Inf,40,60,Inf),
                         right = F, labels=c("20-39","40-59","≥60"))
glm(formula = out~age_cat+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
    family=binomial(), weights=swts, data=pooled.DX) -> m_test1B
glmtoolbox::hltest(m_test1B)
#(ns, p=0.11199)


#Untreated
glm(formula = out~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
      family=binomial(), weights=swts, data=pooled.TX) -> m_test2
glmtoolbox::hltest(m_test2)
#(***, 2.3216e-06)
#Could be fixed by using categorical age instead of age
pooled.TX$age_cat <- cut(pooled.TX$age*10, c(-Inf,40,60,Inf),
                         right = F, labels=c("20-39","40-59","≥60"))
glm(formula = out~age_cat+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
    family=binomial(), weights=swts, data=pooled.TX) -> m_test2B
glmtoolbox::hltest(m_test2B)
#(ns, 0.44869)


#Uncontrolled
glm(formula = (out==0)~age+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
    family=binomial(), weights=swts, data=pooled.GOAL) -> m_test3
glmtoolbox::hltest(m_test3)
#(***, 0.0002988)
#Could be fixed by using categorical age instead of age
pooled.GOAL$age_cat <- cut(pooled.GOAL$age*10, c(-Inf,40,60,Inf),
                         right = F, labels=c("20-39","40-59","≥60"))
glm(formula = (out==0)~age_cat+sex+tab+alc+edu+ind+afi+sli+urb+obe+dm2+yrs,
    family=binomial(), weights=swts, data=pooled.GOAL) -> m_test3B
glmtoolbox::hltest(m_test3B)
#(*, 0.022)

#exp(m_test1B$coefficients[2:3]) %>% round(3)
#exp(m_test2B$coefficients[2:3]) %>% round(3)
#exp(m_test3B$coefficients[2:3]) %>% round(3)


####----------------------- Figures ------------------------------------#### ----####
#### Prevalence --------------------(F1/FS4)- #####
#------------------------------------ ESC/ESH ---------------------------------#
setwd(WD_OUT)

prev <- rbind(
  UHT_prev_year, DX_prev_year, PRE_prev_year, HTN_fin_prev_year) %>%
  mutate("prop2"=round(prop,1),"IC952"=round(IC95,1),
         "lIC952"=round(lIC95,1),"uIC952"=round(uIC95,1),
         "lab"=paste0(prop2,"%","\n","(",lIC952,"-",uIC952,")"),
         "clust"=ordered(cluster, levels=c(
           "Total Hypertension", "Diagnosed Hypertension",
           "Undiagnosed Hypertension", "High-Normal BP")))

f1 <- prev %>% ggplot(aes(
    x=as.numeric(Year), y=prop,group=clust, colour=clust, linetype=clust)) +
  geom_line(linewidth=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2,
                position=position_dodge(0.01), linetype=1) + theme_bw() +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle", colour="", linetype="",
       title="Hypertension prevalence (ESC/ESH definitions)") +
  scale_x_continuous(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_color_manual(values=c("#000000","#AB2F26","#804FB3","#EE7474")) +
  scale_linetype_manual(values=c(1,4,1,4)) + 
  geom_text(data = filter(prev, clust=="Total Hypertension") %>%
              mutate(sp=ifelse(Year%in%c(2021,2024), lIC952-1.3, uIC952+1.5)),
            fontface = "bold.italic", size=4, lineheight=0.75,
            aes(label = lab, y = sp)) + ylim(5,40) +
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13, hjust=0.5),
        panel.grid.minor.x = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=14, hjust=0.5),
        legend.key.size = unit(12, "mm")); f1

ggsave(f1, file="Hypertension phenotypes/Figures/Figure1.pdf", bg="transparent",
       width=37.5, height=22.5, units=c("cm"), dpi=600, limitsize = FALSE)


#------------------------------------ ACC/AHA ---------------------------------#
prev.A <- rbind(
  UHT_prev_year.A, DX_prev_year, PRE_prev_year.A, HTN_fin_prev_year.A) %>%
  mutate("prop2"=round(prop,1),"IC952"=round(IC95,1),
         "lIC952"=round(lIC95,1),"uIC952"=round(uIC95,1),
         "lab"=paste0(prop2,"%","\n","(",lIC952,"-",uIC952,")"),
         "clust"=ordered(cluster, levels=c(
           "Total Hypertension", "Diagnosed Hypertension",
           "Undiagnosed Hypertension", "Elevated BP")))

sf2 <- prev.A %>% ggplot(aes(
  x=as.numeric(Year), y=prop,group=clust, colour=clust, linetype=clust)) +
  geom_line(linewidth=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2,
                position=position_dodge(0.01), linetype=1) + theme_bw() +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle", colour="", linetype="",
       title="Hypertension prevalence (ACC/AHA definitions)") +
  scale_x_continuous(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_color_manual(values=c("#000000","#AB2F26","#804FB3","#EE7474")) +
  scale_linetype_manual(values=c(1,4,1,4)) + 
  geom_text(data = filter(prev.A, clust=="Total Hypertension") %>%
              mutate(sp=ifelse(Year%in%c(2021), lIC952-3.5, uIC952+3.5)),
            fontface = "bold.italic", size=4, lineheight=1,
            aes(label = lab, y = sp)) + ylim(6.5,75) +
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13, hjust=0.5),
        panel.grid.minor.x = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=14, hjust=0.5),
        legend.key.size = unit(12, "mm")); sf2

ggsave(sf2, file="Hypertension phenotypes/Figures/FigureS4.png", dpi=600,
       bg="transparent", width=37.5, height=22.5, units="cm", limitsize = F)


#### Adjusted Prevalence ----------(FS2/FS3)- #####
#------------------------------------ ESC/ESH ---------------------------------#
setwd(WD_OUT)

adj_prev <- rbind(std_NEW_HTN, std_DX, std_HTN_fin) %>%
  rename(Year=cycle, prop=std_prev) %>%
  mutate(
    "prop2"=round(prop,1), "lIC952"=round(lower,1),"uIC952"=round(upper,1),
    "lab"=paste0(prop2,"%","\n","(",lIC952,"-",uIC952,")"),
    "clust"=ordered(cluster, levels=c(
      "Total Hypertension Adjusted","Undiagnosed Hypertension Adjusted",
      "Diagnosed Hypertension Adjusted"))) %>%
  mutate(
    "spy" = case_when(
    clust=="Total Hypertension Adjusted"&Year%in%c(2021) ~ lIC952-1.25,
    clust=="Undiagnosed Hypertension Adjusted"&Year%in%2006:2023 ~ lIC952-1.25,
    clust=="Diagnosed Hypertension Adjusted"&Year%in%2006:2023 ~ uIC952+1.25,
    Year%in%2024 ~ prop, T ~ uIC952+1.25),
    "spx" = case_when(Year==2024 ~ 2024+1, T ~ as.numeric(Year))) %>%
  mutate("sz" = ifelse(Year==2024, 1, 0))

f1_a <- adj_prev %>% 
  ggplot(aes(x=as.numeric(Year), y=prop,group=clust,
             colour=clust, linetype=clust)) +
  geom_line(linewidth=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.01), linetype=1) + theme_bw() +
  scale_x_continuous(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24))+ ggtitle(
    "Age- and sex-adjusted hypertension prevalence (ESC/ESH definitions)") +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle", colour="",linetype="") +
  scale_color_manual(values=c("#000000","#804FB3","#AB2F26")) +
  scale_linetype_manual(values=c(1,1,4)) + 
  geom_text(data = adj_prev, fontface = "bold.italic", show.legend = F,
            lineheight=0.7, aes(label = lab, y = spy, x = spx, size = sz)) +
  scale_size_continuous(range=c(3.15, 3.75)) +
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13, hjust=0.5),
        panel.grid.minor.x = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=14, hjust=0.5),
        legend.key.size = unit(12, "mm")); f1_a

ggsave(f1_a, file="Hypertension phenotypes/Figures/FigureS2.png", limitsize=F,
       bg="transparent", width=37.5, height=22.5, units="cm", dpi=600)


#ESC/ESH phenotypes
toplab <- with(std_NEW_HTN, paste0(
  round(std_prev,1), "% (", round(lower,1),
  "-", round(upper,1), ")")); f2_a <- rbind(
  mutate(std_ISH, x="C"),
  mutate(std_IDH, x="B"),
  mutate(std_SDH, x="A")) %>% 
  rename(Year=cycle) %>% transmute(
    Year, "x"=ordered(x,labels=c("SDH","IDH","ISH")), "p"=std_prev,
    "ic"=(upper-lower)/2) %>% group_by(Year) %>% arrange(desc(x)) %>%
  mutate("up"=cumsum(p)+ic, "low"=cumsum(p)-ic) %>% ungroup() %>%
  mutate("Year"=as.factor(Year), "lab"=paste0(
    round(p,1), "%", " (", round(p-ic,1),"-", round(p+ic,1),")")) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) + theme_bw() +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = low, ymax = up), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="") + ggtitle(
    "Age- and sex-adjusted prevalence of UDH phenotypes (ESC/ESH definitions)")+
  scale_fill_manual(values=c("#552586", "#804FB3", "#B589D6")) + ylim(0,27) +
  geom_text(aes(label = lab), size=4, col="white",
            position = position_stack(vjust = .55), fontface="bold.italic") +
  annotate("label", y=std_NEW_HTN$upper+0.5, x=1:10, lineheight=1,
           label=toplab, size=4, fontface="bold.italic")+
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="inside",
        legend.position.inside = c(0.9, 0.9), legend.direction ="horizontal",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2_a

f2_b <- rbind(
  mutate(std_Untreated,x=1), mutate(std_Treated,x=2)) %>%
  rename(Year=cycle)%>%
  transmute(Year, x, p=std_prev, ic=(upper-lower)/2) %>% 
  mutate(x=ordered(x,labels=c("Untreated", "Treated"))) %>% group_by(Year) %>%
  arrange(desc(x)) %>% mutate(pos = cumsum(p), u=pos+ic, l=pos-ic) %>%
  ungroup() %>% mutate(Year=as.factor(Year), lab = paste0(
    round(p,1), "%\n(", round(p-ic,1), "-", round(p+ic,1), ")"),
    l=ifelse(x=="Treated", l, NA), u=ifelse(x=="Treated", u, NA)) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = l, ymax = u), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_y_continuous(breaks = seq(0,100,25)) +  theme_bw() +
  scale_fill_manual(values=c("#EE7474", "#AB2F26")) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Treatment status (among diagnosed adults)") +
  geom_text(aes(label = lab), size=4, col="white", fontface="bold.italic",
            position = position_stack(vjust = .5), lineheight=1) +
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2_b

#Goal achievement
f2_c <- rbind(
  mutate(std_Uncontrolled,x=1), mutate(std_Controlled,x=2)) %>%
  rename(Year=cycle)%>%
  transmute(Year, x, p=std_prev, ic=(upper-lower)/2) %>% mutate(
    x=ordered(x,labels=c("Uncontrolled","Controlled"))) %>% group_by(Year) %>%
  arrange(desc(x)) %>% mutate(pos = cumsum(p), u=pos+ic, l=pos-ic) %>%
  ungroup() %>% mutate(Year=as.factor(Year), lab = paste0(
    round(p,1), "%\n(", round(p-ic,1), "-", round(p+ic,1), ")"),
    l=ifelse(x=="Controlled", l, NA), u=ifelse(x=="Controlled", u, NA)) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = l, ymax = u), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_y_continuous(breaks = seq(0,100,25)) +  theme_bw() +
  scale_fill_manual(values=c("#4682B4","#104E8B")) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Goal achievement (among diagnosed adults)") +
  geom_text(aes(label = lab), size=4, col="white", fontface="bold.italic",
            position = position_stack(vjust = .5), lineheight=1) +
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2_c

f_2 <- ggpubr::ggarrange(f2_a, ggpubr::ggarrange(
  f2_b, f2_c, ncol=2, labels=c("B", "C"), font.label = list(size=18)),
  nrow=2, labels = "A", common.legend = F, font.label = list(size=18))

ggsave(f_2, file="Hypertension phenotypes/Figures/FigureS3.png", dpi=600,
       bg="transparent", width=42.5*18/16, height=25.5*18/16,
       units="cm", limitsize = F)


#### Phenotypes + TX + control----------(F2)- #####
setwd(WD_OUT)

#ESC/ESH phenotypes
toplab <- with(UHT_prev_year, paste0(
  prop, "% (", round(lIC95,1), "-", round(uIC95,1), ")")); f2a <- rbind(
    mutate(ISH_prev_year_overall, x="C"),
    mutate(IDH_prev_year_overall, x="B"),
    mutate(SDH_prev_year_overall, x="A")) %>% transmute(
      Year, "x"=ordered(x,labels=c("SDH","IDH","ISH")), "p"=prop, "ic"=IC95) %>%
  group_by(Year) %>% arrange(desc(x)) %>%
  mutate("up"=cumsum(p)+ic, "low"=cumsum(p)-ic) %>% ungroup() %>%
  mutate("Year"=as.factor(Year), "lab"=paste0(
    round(p,1), "%", " (", round(p-ic,1),"-", round(p+ic,1),")")) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) + theme_bw() +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = low, ymax = up), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Prevalence of UDH phenotypes (ESC/ESH definitions)") +
  scale_fill_manual(values=c("#552586", "#804FB3", "#B589D6")) +  ylim(0,23) +
  geom_text(aes(label = lab), size=4, col="white",
            position = position_stack(vjust = .55), fontface="bold.italic") +
  annotate("label", y=UHT_prev_year$uIC95+.5, x=1:10, lineheight=1,
           label=toplab, size=4, fontface="bold.italic")+
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="inside",
        legend.position.inside = c(0.9, 0.9), legend.direction ="horizontal",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2a

#Treated vs Untreated
f2b <- rbind(
  mutate(HX_untreated_prev_year,x=1), mutate(HX_treated_prev_year,x=2)) %>%
  transmute(Year, x, p=prop, ic=IC95) %>% filter(Year!=2020) %>%
  mutate(x=ordered(x,labels=c("Untreated", "Treated"))) %>% group_by(Year) %>%
  arrange(desc(x)) %>% mutate(pos = cumsum(p), u=pos+ic, l=pos-ic) %>%
  ungroup() %>% mutate(Year=as.factor(Year), lab = paste0(
    round(p,1), "%\n(", round(p-ic,1), "-", round(p+ic,1), ")"),
    l=ifelse(x=="Treated", l, NA), u=ifelse(x=="Treated", u, NA)) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = l, ymax = u), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_y_continuous(breaks = seq(0,100,25)) +  theme_bw() +
  scale_fill_manual(values=c("#EE7474", "#AB2F26")) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Treatment status (among diagnosed adults)") +
  geom_text(aes(label = lab), size=4, col="white", fontface="bold.italic",
            position = position_stack(vjust = .5), lineheight=1) +
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2b

#Goal achievement
f2c <- rbind(
  mutate(HX_not_goal_prev_year,x=1), mutate(HX_goal_prev_year,x=2)) %>%
  transmute(Year, x, p=prop, ic=IC95) %>% mutate(
    x=ordered(x,labels=c("Uncontrolled","Controlled"))) %>% group_by(Year) %>%
  arrange(desc(x)) %>% mutate(pos = cumsum(p), u=pos+ic, l=pos-ic) %>%
  ungroup() %>% mutate(Year=as.factor(Year), lab = paste0(
    round(p,1), "%\n(", round(p-ic,1), "-", round(p+ic,1), ")"),
    l=ifelse(x=="Controlled", l, NA), u=ifelse(x=="Controlled", u, NA)) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = l, ymax = u), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_y_continuous(breaks = seq(0,100,25)) +  theme_bw() +
  scale_fill_manual(values=c("#4682B4","#104E8B")) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Goal achievement (among diagnosed adults)") +
  geom_text(aes(label = lab), size=4, col="white", fontface="bold.italic",
            position = position_stack(vjust = .5), lineheight=1) +
  theme(plot.title = element_text(size=19, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); f2c

f2 <- ggpubr::ggarrange(
  f2a, ggpubr::ggarrange(f2b, f2c, ncol=2, labels=c("B", "C"),
                 font.label = list(size=18)), nrow=2, labels = "A",
  common.legend = F, font.label = list(size=18))

ggsave(f2, file="Hypertension phenotypes/Figures/Figure2.pdf", dpi=600,
       bg="transparent", width=42.5*18/16, height=25.5*18/16,
       units="cm", limitsize = F)


#### Phenotypes supp. -------------(FS5/FS6)- #####
setwd(WD_OUT)

#------------------------------------ ACC/AHA ---------------------------------#
toplab <- with(UHT_prev_year.A, paste0(
  prop, "% (", round(lIC95,1), "-", round(uIC95,1), ")")); sf3 <- rbind(
    mutate(ISH_prev_year.A, x="C"),
    mutate(IDH_prev_year.A, x="B"),
    mutate(SDH_prev_year.A, x="A")) %>% transmute(
      Year, "x"=ordered(x,labels=c("SDH","IDH","ISH")), "p"=prop, "ic"=IC95) %>%
  group_by(Year) %>% arrange(desc(x)) %>%
  mutate("up"=cumsum(p)+ic, "low"=cumsum(p)-ic) %>% ungroup() %>%
  mutate("Year"=as.factor(Year), "lab"=paste0(
    round(p,1), "%", " (", round(p-ic,1),"-", round(p+ic,1),")")) %>%
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) + theme_bw() +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = low, ymax = up), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Prevalence of UDH phenotypes (ACC/AHA definitions)") +
  scale_fill_manual(values=c("#552586", "#804FB3", "#B589D6")) + ylim(0,60) +
  geom_text(aes(label = lab), size=4, col="white",
            position = position_stack(vjust = .55), fontface="bold.italic") +
  annotate("label", y=UHT_prev_year.A$uIC95+2.5, x=1:10, lineheight=1,
           label=toplab, size=4, fontface="bold.italic")+
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); sf3

ggsave(sf3, file="Hypertension phenotypes/Figures/FigureS5.png", dpi=600,
       bg="transparent", width=37.5*10/9, height=22.5*10/9,
       units="cm", limitsize = F)

#-------------------------------- Proportions ECS/ESH -------------------------#
#Phenotypes of UDH
sf4 <- rbind(
  ISH_prev_year %>% mutate(x="C"), IDH_prev_year %>% mutate(x="B"),
  SDH_prev_year %>% mutate(x="A")) %>% transmute(Year, x, p=prop, ic=IC95) %>%
  mutate(x=ordered(x, labels=c("SDH","IDH","ISH"))) %>% group_by(Year) %>%
  arrange(desc(x)) %>% mutate(pos = cumsum(p), u=pos+ic, l=pos-ic) %>%
  ungroup() %>% mutate(Year=as.factor(Year), lab = paste0(
    round(p,1), "%\n(", round(p-ic,1), "-", round(p+ic,1), ")"),
    l=ifelse(x=="SDH", NA, l), u=ifelse(x=="SDH", NA, u)) %>% 
  ggplot(aes(x=Year, y=p, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=1) +
  geom_errorbar(aes(ymin = l, ymax = u), width = .1, col = "black") +
  scale_x_discrete(breaks = 2000 + c(0,6,12,16,18,20,21,22,23,24)) +
  scale_y_continuous(breaks = seq(0,100,25)) +  theme_bw() +
  scale_fill_manual(values=c("#552586", "#804FB3", "#B589D6")) +
  labs(x="ENSANUT cycle", y="Weighted prevalence (%)", fill="",
       title="Phenotype proportions among UDH (ESC/ESH definitions)") +
  geom_text(aes(label = lab), size=4, col="white", fontface="bold.italic",
            position = position_stack(vjust = .5), lineheight=1) +
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5)); sf4

ggsave(sf4, file="Hypertension phenotypes/Figures/FigureS6.png", dpi=600,
       bg="transparent", width=37.5*10/9, height=22.5*10/9,
       units="cm", limitsize = F)


#### Modifiers of UHT-------------------(F3)- ####
setwd(WD_OUT)

lower <- c(UHT_prev_age$lIC95,UHT_prev_sex$lIC95,
           UHT_prev_bmi$lIC95,UHT_prev_t2d$lIC95) %>% min %>% floor()
upper <- c(UHT_prev_age$uIC95,UHT_prev_sex$uIC95,
           UHT_prev_bmi$uIC95,UHT_prev_t2d$uIC95) %>% max %>% ceiling()

##-- Age --##
plt <- "Age groups"; f3a <- UHT_prev_age %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- Sex --##
plt <- "Sex"; f3b <- UHT_prev_sex %>%
  mutate(group=forcats::fct_rev(group)) %>% 
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- BMI --##
plt <- "BMI obesity"; f3c <- UHT_prev_bmi %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- T2D --##
plt <- "Diabetes status"; f3d <- UHT_prev_t2d %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))


#### Modifiers of ISH------------------(FS7)- ####
setwd(WD_OUT)

lower <- c(ISH_prev_age$lIC95,ISH_prev_sex$lIC95,
           ISH_prev_bmi$lIC95,ISH_prev_t2d$lIC95) %>% min %>% floor()
upper <- c(ISH_prev_age$uIC95,ISH_prev_sex$uIC95,
           ISH_prev_bmi$uIC95,ISH_prev_t2d$uIC95) %>% max %>% ceiling()

##-- Age --##
plt <- "Age groups"; sf5a <- ISH_prev_age %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- Sex --##
plt <- "Sex"; sf5b <- ISH_prev_sex %>%
  mutate(group=forcats::fct_rev(group)) %>% 
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- BMI --##
plt <- "BMI obesity"; sf5c <- ISH_prev_bmi %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- T2D --##
plt <- "Diabetes status"; sf5d <- ISH_prev_t2d %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))


#### Modifiers of IDH------------------(FS8)- ####
setwd(WD_OUT)

lower <- c(IDH_prev_age$lIC95,IDH_prev_sex$lIC95,
           IDH_prev_bmi$lIC95,IDH_prev_t2d$lIC95) %>% min %>% floor()
upper <- c(IDH_prev_age$uIC95,IDH_prev_sex$uIC95,
           IDH_prev_bmi$uIC95,IDH_prev_t2d$uIC95) %>% max %>% ceiling()

##-- Age --##
plt <- "Age groups"; sf6a <- IDH_prev_age %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- Sex --##
plt <- "Sex"; sf6b <- IDH_prev_sex %>%
  mutate(group=forcats::fct_rev(group)) %>% 
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- BMI --##
plt <- "BMI obesity"; sf6c <- IDH_prev_bmi %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- T2D --##
plt <- "Diabetes status"; sf6d <- IDH_prev_t2d %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))


#### Modifiers of SDH------------------(FS9)- ####
setwd(WD_OUT)

lower <- c(SDH_prev_age$lIC95,SDH_prev_sex$lIC95,
           SDH_prev_bmi$lIC95,SDH_prev_t2d$lIC95) %>% min %>% floor()
upper <- c(SDH_prev_age$uIC95,SDH_prev_sex$uIC95,
           SDH_prev_bmi$uIC95,SDH_prev_t2d$uIC95) %>% max %>% ceiling()

##-- Age --##
plt <- "Age groups"; sf7a <- SDH_prev_age %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- Sex --##
plt <- "Sex"; sf7b <- SDH_prev_sex %>%
  mutate(group=forcats::fct_rev(group)) %>% 
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- BMI --##
plt <- "BMI obesity"; sf7c <- SDH_prev_bmi %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))

##-- T2D --##
plt <- "Diabetes status"; sf7d <- SDH_prev_t2d %>%
  ggplot(aes(x=Year, y=prop, group=group, colour=group)) +
  geom_line(linewidth=1.5) + geom_point(size=2) + geom_errorbar(aes(
    ymin=lIC95, ymax=uIC95), width=.2, position=position_dodge(0.01)) +
  scale_x_continuous(breaks = 2000+c(0,6,12,16,18,20,21,22,23,24), labels = c(
    paste0("20",c("00","06","12","16","18")), 20:24)) +
  labs(y="Weighted prevalence (%)", x="ENSANUT cycle",
       colour=NULL, title=plt) + theme_bw() + ylim(lower,upper) +
  scale_color_manual(values=c("#B498BF", "#353134")) + 
  theme(plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5), legend.position="bottom",
        panel.grid.minor.x = element_blank(), legend.key.size = unit(10,"mm"),
        legend.text = element_text(size=14, hjust=0.5))


#### Modifiers: arrange & export------------- ####
setwd(WD_OUT)

f3 <- ggpubr::ggarrange(
  f3a, f3b, f3c, f3d, labels = LETTERS[1:4], nrow=2, ncol=2,
  font.label = list(size=18)) %>% annotate_figure(top = text_grob(
    "Modifiers of Undiagnosed Hypertension Prevalence\n",
    face = "bold", size = 19, lineheight = 1))
ggsave(f3, file="Hypertension phenotypes/Figures/Figure3.pdf",  dpi=600,
       bg="transparent", width=40, height=24, units="cm",
       limitsize = F, device = cairo_pdf)


sf5 <- ggpubr::ggarrange(
  sf5a, sf5b, sf5c, sf5d, labels = LETTERS[1:4], nrow=2, ncol=2,
  font.label = list(size=18)) %>% annotate_figure(top = text_grob(
    "Modifiers of ISH Prevalence\n", face = "bold", size = 19, lineheight = 1))
ggsave(sf5, file="Hypertension phenotypes/Figures/FigureS7.png",  dpi=600,
       bg="transparent", width=40, height=24, units="cm", limitsize = F)

sf6 <- ggpubr::ggarrange(
  sf6a, sf6b, sf6c, sf6d, labels = LETTERS[1:4], nrow=2, ncol=2,
  font.label = list(size=18)) %>% annotate_figure(top = text_grob(
    "Modifiers of IDH Prevalence\n", face = "bold", size = 19, lineheight = 1))
ggsave(sf6, file="Hypertension phenotypes/Figures/FigureS8.png",  dpi=600,
       bg="transparent", width=40, height=24, units="cm", limitsize = F)

sf7 <- ggpubr::ggarrange(
  sf7a, sf7b, sf7c, sf7d, labels = LETTERS[1:4], nrow=2, ncol=2,
  font.label = list(size=18)) %>% annotate_figure(top = text_grob(
    "Modifiers of SDH Prevalence\n", face = "bold", size = 19, lineheight = 1))
ggsave(sf7, file="Hypertension phenotypes/Figures/FigureS9.png",  dpi=600,
       bg="transparent", width=40, height=24, units="cm", limitsize = F)



#### Logistic models--------------------(F4)- ####
setwd(WD_OUT)

#------------------------------------------------------------------------------#
#------------------------------- UNDIAGNOSED ----------------------------------#
#------------------------------------------------------------------------------#
#Get coefficients from logistic model
summary(binom.01)$coefficients[-1,1:2] %>% as.data.frame() %>%
  `colnames<-`(c("est", "se")) %>% rownames_to_column("term") %>% transmute(
    "est"=est, "se"=se, "var"=substr(term, 1, 3), "term"=ifelse(
      nchar(term)<=3, term, substr(term, 4, 100))) %>%
  mutate("term"=case_when(term=="age"~"(10 year\nincrements)",
                          term=="yrs"~"(Years\nsince 2000)", TRUE~term)) %>%
  group_by(var) %>% mutate("term_order"=row_number()) %>% ungroup() -> coef.01

#Get reference levels for categorical variables
model.frame(binom.01) %>% lapply(\(x) if (is.factor(x)) levels(x)[1]) %>%  
  .[!sapply(., is.null)] %>% unlist %>% as.data.frame() %>%
  `colnames<-`("term") %>% rownames_to_column("var") %>%
  transmute(est=0, se=0, var=var, term=term, term_order=0) %>%
  as_tibble() -> term.01

#Merge, arrange, and get odds ratios
rbind(coef.01, term.01) %>% mutate(
  "var_order" = match(var, unique(var)),
  "lower"=est-qnorm(0.975)*se, "upper"=est+qnorm(0.975)*se,
  "OR"=exp(est), "OR.l"=exp(lower), "OR.u"=exp(upper),
  "label" = ifelse(term == "", var, paste0("  ", term)),
  "label" = factor(label, levels = rev(label))) -> data.01

#Get proper variable names
model.names <- model.frame(binom.01) %>% (\(x) x[-c(1, length(x))]) %>%
  names; select(pool.binomial, all_of(model.names)) %>% as.list %>%
  lapply(attr, "label") %>% Filter(Negate(is.null), .) %>%
  unlist %>% data.frame %>% `colnames<-`("var_name") %>% 
  rownames_to_column("var") -> name.01


#------------------------------------------------------------------------------#
#------------------------------- UNTREATED ------------------------------------#
#------------------------------------------------------------------------------#
#Get coefficients from logistic model
summary(binom.02)$coefficients[-1,1:2] %>% as.data.frame() %>%
  `colnames<-`(c("est", "se")) %>% rownames_to_column("term") %>% transmute(
    "est"=est, "se"=se, "var"=substr(term, 1, 3), "term"=ifelse(
      nchar(term)<=3, term, substr(term, 4, 100))) %>%
  mutate("term"=case_when(term=="age"~"(10 year\nincrements)",
                          term=="yrs"~"(Years\nsince 2000)", TRUE~term)) %>%
  group_by(var) %>% mutate("term_order"=row_number()) %>% ungroup() -> coef.02

#Get reference levels for categorical variables
model.frame(binom.02) %>% lapply(\(x) if (is.factor(x)) levels(x)[1]) %>%  
  .[!sapply(., is.null)] %>% unlist %>% as.data.frame() %>%
  `colnames<-`("term") %>% rownames_to_column("var") %>%
  transmute(est=0, se=0, var=var, term=term, term_order=0) %>%
  as_tibble() -> term.02

#Merge, arrange, and get odds ratios
rbind(coef.02, term.02) %>% mutate(
  "var_order" = match(var, unique(var)),
  "lower"=est-qnorm(0.975)*se, "upper"=est+qnorm(0.975)*se,
  "OR"=exp(est), "OR.l"=exp(lower), "OR.u"=exp(upper),
  "label" = ifelse(term == "", var, paste0("  ", term)),
  "label" = factor(label, levels = rev(label))) -> data.02

#Get proper variable names
model.names <- model.frame(binom.02) %>% (\(x) x[-c(1, length(x))]) %>%
  names; select(pool.binomial, all_of(model.names)) %>% as.list %>%
  lapply(attr, "label") %>% Filter(Negate(is.null), .) %>%
  unlist %>% data.frame %>% `colnames<-`("var_name") %>% 
  rownames_to_column("var") -> name.02


#------------------------------------------------------------------------------#
#------------------------------- UNCONTROLLED ---------------------------------#
#------------------------------------------------------------------------------#
#Get coefficients from logistic model
summary(binom.03)$coefficients[-1,1:2] %>% as.data.frame() %>%
  `colnames<-`(c("est", "se")) %>% rownames_to_column("term") %>% transmute(
    "est"=est, "se"=se, "var"=substr(term, 1, 3), "term"=ifelse(
      nchar(term)<=3, term, substr(term, 4, 100))) %>%
  mutate("term"=case_when(term=="age"~"(10 year\nincrements)",
                          term=="yrs"~"(Years\nsince 2000)", TRUE~term)) %>%
  group_by(var) %>% mutate("term_order"=row_number()) %>% ungroup() -> coef.03

#Get reference levels for categorical variables
model.frame(binom.03) %>% lapply(\(x) if (is.factor(x)) levels(x)[1]) %>%  
  .[!sapply(., is.null)] %>% unlist %>% as.data.frame() %>%
  `colnames<-`("term") %>% rownames_to_column("var") %>%
  transmute(est=0, se=0, var=var, term=term, term_order=0) %>%
  as_tibble() -> term.03

#Merge, arrange, and get odds ratios
rbind(coef.03, term.03) %>% mutate(
  "var_order" = match(var, unique(var)),
  "lower"=est-qnorm(0.975)*se, "upper"=est+qnorm(0.975)*se,
  "OR"=exp(est), "OR.l"=exp(lower), "OR.u"=exp(upper),
  "label" = ifelse(term == "", var, paste0("  ", term)),
  "label" = factor(label, levels = rev(label))) -> data.03

#Get proper variable names
model.names <- model.frame(binom.03) %>% (\(x) x[-c(1, length(x))]) %>%
  names; select(pool.binomial, all_of(model.names)) %>% as.list %>%
  lapply(attr, "label") %>% Filter(Negate(is.null), .) %>%
  unlist %>% data.frame %>% `colnames<-`("var_name") %>% 
  rownames_to_column("var") -> name.03


#------------------------------------------------------------------------------#
#------------------------------- OR PLOT --------------------------------------#
#------------------------------------------------------------------------------#
fmt4 <- function(x){sprintf(x, fmt="%#.2f")}

#Set labels and positions
merge(data.01, name.01, by="var") %>% mutate(
  "term" = reorder(term, -term_order),
  "var_name" = reorder(var_name, var_order),
  "lab0"=ifelse(term_order==0, NA, paste0(
    fmt4(OR), " (",fmt4(OR.l)," - ",fmt4(OR.u),")")),
  "lab1"=ifelse(lower*upper<0, NA, lab0),
  "lab2"=ifelse(lower*upper>0, NA, lab0),
  "pos"=0.225) -> df.01

#Set labels and positions
merge(data.02, name.02, by="var") %>% mutate(
  "term" = reorder(term, -term_order),
  "var_name" = reorder(var_name, var_order),
  "lab0"=ifelse(term_order==0, NA, paste0(
    fmt4(OR), " (",fmt4(OR.l)," - ",fmt4(OR.u),")")),
  "lab1"=ifelse(lower*upper<0, NA, lab0),
  "lab2"=ifelse(lower*upper>0, NA, lab0),
  "pos"=0.225)-> df.02

#Set labels and positions
merge(data.03, name.03, by="var") %>% mutate(
  "term" = reorder(term, -term_order),
  "var_name" = reorder(var_name, var_order),
  "lab0"=ifelse(term_order==0, NA, paste0(
    fmt4(OR), " (",fmt4(OR.l)," - ",fmt4(OR.u),")")),
  "lab1"=ifelse(lower*upper<0, NA, lab0),
  "lab2"=ifelse(lower*upper>0, NA, lab0),
  "pos"=0.225)-> df.03

#Plot
rbind(mutate(df.01, out=1), mutate(df.02, out=2), mutate(df.03, out=3)) %>%
  mutate("out" = factor(out, 1:3, paste0(c(
    "Undiagnosed", "Untreated", "Uncontrolled"), " hypertension") )) %>%
  ggplot(aes(x = OR, y = term, color = out)) +
  facet_grid(rows=vars(var_name), cols=vars(out), scales="free_y",switch="y") +
  geom_vline(xintercept=1, linetype=2) + theme_bw() +
  scale_x_log10(breaks=c(0.50, 0.75, 1.00, 1.3, 2)) +
  geom_pointrange(aes(xmin = OR.l, xmax = OR.u),
                  position = position_dodge(width = 0.5), size = 0.5) +
  labs(x="Odds Ratio (95% CI)", y=NULL) +
  scale_color_manual(values=c("#552586", "#E25E5E", "#4682B4")) +
  geom_text(aes(label = lab1, x = pos, y = term), size=4.5,
             fontface="bold.italic", hjust=0) +
  geom_text(aes(label = lab2, x = pos, y = term), size=4.5,
            fontface="plain", hjust=0, color="black") +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        strip.text = element_text(face="bold.italic", size=14),
        plot.title = element_text(size=17, face="bold", hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=2),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=-1),
        axis.text = element_text(size=13,hjust=0.5)) -> f4; f4

ggsave(f4, file="Hypertension phenotypes/Figures/Figure4.pdf",  dpi=600,
       bg="transparent", width=51, height=30, units="cm", limitsize = F)





#### Missing data analysis------------(FS10)- ####
setwd(WD_OUT)

vars <- c("HTN_fin", "NEW_HTN", "Year", "Age", "Sex2", "Area2", "SBP",
          "BMI_cat", "Diabetes3", "Alc", "Smoking", "SocSec2", "Edu",
          "DISLI_cat2", "Indigenous2", "DX", "HX_untreated", "HTN_goal")

names <- c("Total_HTN", "Undiagnosed_HTN", "Cycle",  "Age", "Sex", "Area",
           "SBP", "BMI-\nobesity", "Prior\nTD2", "Alcohol\nintake", "Smoking",
           "Social\nsecurity", "Edu.\nlevel", "DISLI", "Ind.\nidentity",
           "Diagnosed_HTN", "Untreated_HTN", "Uncontrolled_HTN")

miss.data <- rbind(
  ens00_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens06_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens12_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens16_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens18_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens20_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens21_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens22_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens23_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names),
  ens24_fin.2 %>% select(all_of(vars)) %>% `colnames<-`(names)) %>% 
  mutate("Uncontrolled_HTN"=ifelse(Uncontrolled_HTN==0, 1, 0)) 

n <- format(nrow(miss.data), big.mark=","); miss.data %>%
  mutate("fct"=quantcut(Age, q = 10)) %>% select(
    -Sex, -Age, -Undiagnosed_HTN, -Uncontrolled_HTN, -Total_HTN, -Diagnosed_HTN,
    -Untreated_HTN, -Cycle, -SBP, -DISLI) %>% gg_miss_fct(fct) +
  labs(y="Variable", x="Age deciles (years)", fill="% Missing", title=paste0(
    "Adults 20-99 years old\nn = ", n)) +
  scale_fill_gradient(low = "gray80", high="black", limits=c(0,10)) + theme(
    legend.position = "bottom", axis.text.x = element_text(angle=0, hjust=0.5),
    plot.title = element_text(hjust=0.5, face="bold")) -> miss.1A

n <- format(nrow(miss.data), big.mark=","); miss.data %>%
  mutate("fct"=quantcut(SBP, q = 10)) %>% select(
    -Sex, -Age, -Undiagnosed_HTN, -Uncontrolled_HTN, -Total_HTN, -Diagnosed_HTN,
    -Untreated_HTN, -Cycle, -SBP, -DISLI) %>% gg_miss_fct(fct) +
  labs(y="Variable", x="SBP deciles (years)", fill="% Missing", title=paste0(
    "Adults 20-99 years old\nn = ", n)) +
  scale_fill_gradient(low = "gray80", high="black", limits=c(0,10)) + theme(
    legend.position = "bottom", axis.text.x = element_text(angle=0, hjust=0.5),
    plot.title = element_text(hjust=0.5, face="bold")) -> miss.1B

n <- format(nrow(miss.data), big.mark=","); col_custom<-"black"; miss.data %>%
  select(-Age, -Undiagnosed_HTN, -Uncontrolled_HTN, -Total_HTN,
         -Diagnosed_HTN, -Untreated_HTN, -Cycle, -SBP, -DISLI) %>% 
  gg_miss_var(show_pct = T, facet = Sex) +
  geom_bar(aes(y = pct_miss), stat = "identity", position = "dodge",
           width = 0.025, colour = col_custom, fill = col_custom) +
  geom_point(aes(y = pct_miss), colour = col_custom, fill = col_custom) +
  theme_bw() + scale_y_continuous(lim=c(0,10)) + labs(title=paste0(
    "Adults 20-99 years old\nn = ", n)) +
  geom_label(aes(label=round(as.numeric(pct_miss),1), y=pct_miss+1)) +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(hjust=0.5, face="bold"),
        plot.title = element_text(hjust=0.5, face="bold")) -> miss.1C

n <- format(nrow(miss.data %>% filter(Total_HTN==1)),
            big.mark=","); col_custom <- "#552586"; miss.data %>%
  mutate(Undiagnosed_HTN=ifelse(
    Undiagnosed_HTN==1, "Undiagnosed HTN", "Diagnosed HTN")) %>% select(
      -Age, -Sex, -Uncontrolled_HTN, -Total_HTN, -Diagnosed_HTN,
      -Untreated_HTN, -Cycle, -SBP, -DISLI) %>% 
  gg_miss_var(show_pct = T, facet = Undiagnosed_HTN) + 
  geom_bar(aes(y = pct_miss), stat = "identity", position = "dodge",
           width = 0.025, colour = col_custom, fill = col_custom) +
  geom_point(aes(y = pct_miss), colour = col_custom, fill = col_custom) +
  theme_bw() + scale_y_continuous(lim=c(0,10)) + labs(title=paste0(
    "Adults 20-99 years old with hypertension (HTN)\nn = ", n)) +
  geom_label(aes(label=round(as.numeric(pct_miss),1), y=pct_miss+1)) +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(hjust=0.5, face="bold"),
        plot.title = element_text(hjust=0.5, face="bold")) -> miss.1D

n <- format(nrow(miss.data %>% filter(Diagnosed_HTN==1)),
            big.mark=","); col_custom <- "#E25E5E"; miss.data %>%
  mutate(Untreated_HTN=ifelse(
    Untreated_HTN==1, "Untreated HTN", "Treated HTN")) %>% filter(
      !is.na(Untreated_HTN), Diagnosed_HTN==1) %>% select(
        -Age, -Sex, -Uncontrolled_HTN, -Total_HTN, -Diagnosed_HTN,
        -Undiagnosed_HTN, -Cycle, -SBP, -DISLI) %>% 
  gg_miss_var(show_pct = T, facet = Untreated_HTN) +
  geom_bar(aes(y = pct_miss), stat = "identity", position = "dodge",
           width = 0.025, colour = col_custom, fill = col_custom) +
  geom_point(aes(y = pct_miss), colour = col_custom, fill = col_custom) +
  theme_bw() + scale_y_continuous(lim=c(0,10)) + labs(title=paste0(
    "Adults 20-99 years old with diagnosed HTN\nn = ", n)) +
  geom_label(aes(label=round(as.numeric(pct_miss),1), y=pct_miss+1)) +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(hjust=0.5, face="bold"),
        plot.title = element_text(hjust=0.5, face="bold")) -> miss.1E

n <- format(nrow(miss.data %>% filter(Diagnosed_HTN==1)),
            big.mark=","); col_custom <- "#4682B4"; miss.data %>%
  mutate(Uncontrolled_HTN=ifelse(
    Uncontrolled_HTN==1, "Uncontrolled HTN", "Controlled HTN")) %>% filter(
      !is.na(Untreated_HTN), Diagnosed_HTN==1) %>% select(
        -Age, -Sex, -Total_HTN, -Diagnosed_HTN, -Untreated_HTN,
        -Undiagnosed_HTN, -Cycle, -SBP, -DISLI) %>% 
  gg_miss_var(show_pct = T, facet = Uncontrolled_HTN) +
  geom_bar(aes(y = pct_miss), stat = "identity", position = "dodge",
           width = 0.025, colour = col_custom, fill = col_custom) +
  geom_point(aes(y = pct_miss), colour = col_custom, fill = col_custom) +
  theme_bw() + scale_y_continuous(lim=c(0,10)) + labs(title=paste0(
    "Adults 20-99 years old with diagnosed HTN\nn = ", n)) +
  geom_label(aes(label=round(as.numeric(pct_miss),1), y=pct_miss+1)) +
  theme(axis.title.y = element_blank(),
        strip.text = element_text(hjust=0.5, face="bold"),
        plot.title = element_text(hjust=0.5, face="bold")) -> miss.1F

ggarrange(miss.1A, miss.1B, miss.1C, "", "", "", miss.1D, miss.1E, miss.1F,
          nrow=3, ncol=3, heights = c(1, 0.05, 1),
          labels=c(LETTERS[1:3], rep("",3), LETTERS[4:6])) -> Fig.Miss

ggsave(plot = Fig.Miss, width = 65*.8, height=40*.75, units="cm", limitsize = F,
       filename = "Hypertension phenotypes/Figures/FigureS10.png", dpi = 600)


####----------------------- Tables -------------------------------------#### ----####
#### Population tables---------(T1/TS2-TS11)- ####

#------------------------------------ Labels ----------------------------------#
ens00_fin.2$BMI_cat2 <- ifelse(ens00_fin.2$BMI_cat=="Obesity", 1, 0)
ens06_fin.2$BMI_cat2 <- ifelse(ens06_fin.2$BMI_cat=="Obesity", 1, 0)
ens12_fin.2$BMI_cat2 <- ifelse(ens12_fin.2$BMI_cat=="Obesity", 1, 0)
ens16_fin.2$BMI_cat2 <- ifelse(ens16_fin.2$BMI_cat=="Obesity", 1, 0)
ens18_fin.2$BMI_cat2 <- ifelse(ens18_fin.2$BMI_cat=="Obesity", 1, 0)
ens20_fin.2$BMI_cat2 <- ifelse(ens20_fin.2$BMI_cat=="Obesity", 1, 0)
ens21_fin.2$BMI_cat2 <- ifelse(ens21_fin.2$BMI_cat=="Obesity", 1, 0)
ens22_fin.2$BMI_cat2 <- ifelse(ens22_fin.2$BMI_cat=="Obesity", 1, 0)
ens23_fin.2$BMI_cat2 <- ifelse(ens23_fin.2$BMI_cat=="Obesity", 1, 0)
ens24_fin.2$BMI_cat2 <- ifelse(ens24_fin.2$BMI_cat=="Obesity", 1, 0)

vars <- c(
  "Age", "Sex", "Smoking", "Alc", "Edu", "Indigenous", "SocSec2",
  "DISLI_cat2", "SBP", "DBP", "WHtR", "BMI", "BMI_cat2", "Glucose",
  "HX_T2D", "HX_CKD", "HX_CVD")

labs <- c(
  "Age (years)", "Female (%)", "Smoking (%)", "Alcohol intake (%)",
  "Education level (%)", "Indigenous identity (%)", "Social security (%)",
  "High DISLI (%)", "SBP (mmHg)", "DBP (mmHg)", "WHtR", "BMI (kg/m^2)",
  "BMI-defined obesity (%)", "Glucose (mg/dL)", "Diabetes (%)",
  "CKD (%)", "CVD (%)")

for(i in 1:length(vars)){setattr(ens00_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens06_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens12_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens16_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens18_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens20_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens21_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens22_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens23_fin.2[[vars[i]]], "label", labs[i])}
for(i in 1:length(vars)){setattr(ens24_fin.2[[vars[i]]], "label", labs[i])}

#------------------------------------ Table 1 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2000_1 <- ens00_fin.2%>%dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2006_1 <- ens06_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2012_1 <- ens12_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2016_1 <- ens16_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2018_1 <- ens18_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2020_1 <- ens20_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels();tab_2021_1 <- ens21_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels() ;tab_2022_1 <- ens22_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels() ;tab_2023_1 <- ens23_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels(); tab_2024_1 <- ens24_fin.2%>%
  dplyr::select(dplyr::all_of(vars)) %>%
  tbl_summary(missing = "no", missing_text = "Missing") %>%
  bold_labels()

#Save table
tab2 <-tbl_merge(tbls = list(
  tab_2000_1, tab_2006_1, tab_2012_1, tab_2016_1, tab_2018_1, tab_2020_1,
  tab_2021_1, tab_2022_1, tab_2023_1, tab_2024_1), tab_spanner = c(
    "ENSA 2000", "ENSANUT 2006", "ENSANUT 2012",
    "ENSANUT 2016","ENSANUT 2018","ENSANUT 2020",
    "ENSANUT 2021","ENSANUT 2022","ENSANUT 2023", "ENSANUT 2024")) %>%
  as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%
  body_end_section_landscape() %>% print(
    target = "Hypertension phenotypes/Tables/Table1.docx"); remove(
      tab_2000_1, tab_2006_1, tab_2012_1, tab_2016_1, tab_2018_1, 
      tab_2020_1, tab_2021_1, tab_2022_1, tab_2023_1, tab2)

#------------------------------------ 2000 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2000 <- ens00_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSA 2000**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2000, path="Hypertension phenotypes/Tables/TableS2.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2006 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2006 <- ens06_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2006**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2006, path="Hypertension phenotypes/Tables/TableS3.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2012 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2012 <- ens12_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2012**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2012, path="Hypertension phenotypes/Tables/TableS4.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2016 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2016 <- ens16_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2016**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2016, path="Hypertension phenotypes/Tables/TableS5.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2018 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2018 <- ens18_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2018**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2018, path="Hypertension phenotypes/Tables/TableS6.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2020 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2020 <- ens20_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2020**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2020, path="Hypertension phenotypes/Tables/TableS7.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2021 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2021 <- ens21_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2021**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2021, path="Hypertension phenotypes/Tables/TableS8.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2022 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2022 <- ens22_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2022**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2022, path="Hypertension phenotypes/Tables/TableS9.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2023 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2023 <- ens23_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2023**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2023, path="Hypertension phenotypes/Tables/TableS10.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ 2024 ------------------------------------#
set_flextable_defaults(
  font.size = 9, font.family = "Arial",
  table.layout="autofit", split=F, table_align="center"
); tab_2024 <- ens24_fin.2%>%dplyr::select(categories, dplyr::all_of(vars)) %>%
  tbl_summary(by = categories, missing = "ifany", missing_text = "Missing") %>%
  bold_labels() %>% add_p() %>% add_overall() %>% modify_spanning_header(
    all_stat_cols() ~ "**ENSANUT 2024**") %>% as_flex_table() %>% 
  align(align = "center", part = "all"); save_as_docx(
    tab_2024, path="Hypertension phenotypes/Tables/TableS11.docx",
    pr_section = prop_section(page_size = page_size(orient = "landscape")))

#------------------------------------ TEXT ------------------------------------#
paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list %>% lapply(get) %>% lapply(select, Sex, categories) %>%#lapply(table)
  sapply(function(x){round(prop.table(table(x),2), 3)[2,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, Smoking, categories) %>%
  sapply(function(x){round(prop.table(table(x),2), 3)[1,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, Alc, categories) %>%
  sapply(function(x){round(prop.table(table(x),2), 3)[1,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, Edu, categories) %>%
  sapply(function(x){round(prop.table(table(x),2), 3)[4,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, Indigenous, categories) %>%
  sapply(function(x){round(prop.table(table(x),2), 3)[2,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, SocSec2, categories) %>% 
  sapply(function(x){round(prop.table(table(x),2), 3)[2,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))

paste0("ens", c("00","06","12","16","18", 20:24), "_fin.2") %>%
  as.list() %>% lapply(get) %>% lapply(select, DISLI_cat2, categories) %>%
  sapply(function(x){round(prop.table(table(x),2), 3)[2,] * 100}) %>%
  `colnames<-`(c("00","06","12","16","18", 20:24))


#### Prevalence: year (E/A)------(TS12/TS13)- ####

displayP.CI <- function(x){
  a<-sprintf(fmt="%#.1f",x$prop)
  b<-sprintf(fmt="%#.1f",x$lIC95)
  c<-sprintf(fmt="%#.1f",x$uIC95)
  paste0(a, " (", b, " - ", c, ")")}


#------------------------------------ESC/ESH-----------------------------------#
C1 <- HTN_fin_prev_year$Year %>% as.character()
C2 <- displayP.CI(HTN_fin_prev_year)
C3 <- displayP.CI(UHT_prev_year)
C4 <- displayP.CI(DX_prev_year)
C5 <- displayP.CI(PRE_prev_year)
tabla1 <- data.frame(C1, C2, C3, C4, C5)
names(tabla1) <- c("Year", "Total Hypertension", "Undiagnosed Hypertension",
                   "Diagnosed Hypertension", "High-normal")
tabla1 <- rbind(
  tabla1,
  c("2000-2024 PR (95% CI)", paste0("poiss.0",c(3,5,7,1),"f") %>%
      lapply(get) %>% lapply(PR_fin) %>% as.character),
  c("Adjusted PR (95% CI)", paste0("poiss.adj.",1:4) %>%
      lapply(get) %>% lapply(PR_fin) %>% as.character),
  c("2016-2024 PR (95% CI)", paste0("poiss.0",c(3,5,7,1),"c") %>%
      lapply(get) %>% lapply(PR_fin) %>% as.character))

tabla_ESC <- flextable(tabla1, cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_ESC, path="Hypertension phenotypes/Tables/TableS12.docx")


#------------------------------------ACC/AHA-----------------------------------#
C1 <- HTN_fin_prev_year.A$Year %>% as.character()
C2 <- displayP.CI(HTN_fin_prev_year.A)
C3 <- displayP.CI(UHT_prev_year.A)
C4 <- displayP.CI(DX_prev_year)
C5 <- displayP.CI(PRE_prev_year.A)
tabla0 <- data.frame(C1, C2, C3, C4, C5)
names(tabla0) <- c("Year", "Total Hypertension", "Undiagnosed Hypertension",
                   "Diagnosed Hypertension", "Elevated BP")

tabla0 <- rbind(
  tabla0,
  c("2000-2024 PR (95% CI)", paste0("poiss.0",c(4,6,7,2),"f") %>%
      lapply(get) %>% lapply(PR_fin) %>% as.character),
  c("2016-2024 PR (95% CI)", paste0("poiss.0",c(4,6,7,2),"c") %>%
      lapply(get) %>% lapply(PR_fin) %>% as.character))

tabla_AHA <- flextable(tabla0, cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_AHA, path="Hypertension phenotypes/Tables/TableS13.docx")



#### Prevalence: year+sex ESC/ESH-----(TS14)- ####
C1 <- HTN_fin_prev_sex$Year %>% as.character()
C2 <- HTN_fin_prev_sex$group %>% as.character()
C3 <- displayP.CI(HTN_fin_prev_sex)
C4 <- displayP.CI(UHT_prev_sex)
C5 <- displayP.CI(DX_prev_sex)
C6 <- displayP.CI(PRE_prev_sex)
tabla2 <- data.frame(C1, C2, C3, C4, C5, C6)
names(tabla2) <- c(
  "Year", "Sex", "Total Hypertension", "Undiagnosed Hypertension",
  "Diagnosed Hypertension", "High-normal")
tabla_sex <- flextable(tabla2[1:6], cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_sex, path="Hypertension phenotypes/Tables/TableS14.docx")

#### Prevalence: year+age ESC/ESH-----(TS15)- ####
C1 <- HTN_fin_prev_age$Year %>% as.character()
C2 <- HTN_fin_prev_age$group %>% as.character()
C3 <- displayP.CI(HTN_fin_prev_age)
C4 <- displayP.CI(UHT_prev_age)
C5 <- displayP.CI(DX_prev_age)
C6 <- displayP.CI(PRE_prev_age)
tabla3 <- data.frame(C1, C2, C3, C4, C5, C6)
names(tabla3) <- c(
  "Year", "Age", "Total Hypertension", "Undiagnosed Hypertension",
  "Diagnosed Hypertension", "High-normal")
tabla_age <- flextable(tabla3, cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_age, path="Hypertension phenotypes/Tables/TableS15.docx")

#### Prevalence: year+BMI ESC/ESH-----(TS16)- ####
C1 <- HTN_fin_prev_bmi$Year %>% as.character()
C2 <- HTN_fin_prev_bmi$group %>% as.character()
C3 <- displayP.CI(HTN_fin_prev_bmi)
C4 <- displayP.CI(UHT_prev_bmi)
C5 <- displayP.CI(DX_prev_bmi)
C6 <- displayP.CI(PRE_prev_bmi)
tabla4 <- data.frame(C1, C2, C3, C4, C5, C6)
names(tabla4) <- c(
  "Year", "BMI", "Total Hypertension", "Undiagnosed Hypertension",
  "Diagnosed Hypertension", "High-normal")

tabla_bmi <- flextable(tabla4, cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_bmi, path="Hypertension phenotypes/Tables/TableS16.docx")


#### Prevalence: year+T2D ESC/ESH-----(TS17)- ####
C1 <- HTN_fin_prev_t2d$Year %>% as.character()
C2 <- HTN_fin_prev_t2d$group %>% as.character()
C3 <- displayP.CI(HTN_fin_prev_t2d)
C4 <- displayP.CI(UHT_prev_t2d)
C5 <- displayP.CI(DX_prev_t2d)
C6 <- displayP.CI(PRE_prev_t2d)
tabla5 <- data.frame(C1, C2, C3, C4, C5, C6)
names(tabla5) <- c(
  "Year", "Diabetes", "Total Hypertension", "Undiagnosed Hypertension",
  "Diagnosed Hypertension", "High-normal")
tabla_dm <- flextable(tabla5[1:6], cwidth = c(2,1.5,1.5)) %>% align(
  align = "center", part = "all") %>% autofit(); flextable::save_as_docx(
    tabla_dm, path="Hypertension phenotypes/Tables/TableS17.docx")

####----------------------- Additional analyses ------------------------#### ----####
#### 2000-2012 Sen/Spe adjustment----------------------- ####

#------------------------------------------------------ Undiagnosed/Total -----#
#Sensitivity and specificity of digital device (vs. manual device)
se <- 0.749
sp <- 0.922
D <- se + sp - 1

#Undiagnosed, diagnosed and total hypertension data sets (all years)
uht <- UHT_prev_year
dht <- DX_prev_year
tht <- HTN_fin_prev_year

#Prevalence estimates for 2000, 2006, 2012
P.u <- uht[1:3,1]/100
P.d <- dht[1:3,1]/100
P.t <- tht[1:3,1]/100

#Variance estimates for 2000, 2006, 2012
V.u <- (uht[1:3,2]/qnorm(.975)/100)**2
V.d <- (dht[1:3,2]/qnorm(.975)/100)**2

#Undiagnosed hypertension: adjusted prevalence and variance (delta method)
P.u_ <- (P.u + sp - 1)/(se + sp - 1)
V.u_ <- V.u/D**2

#Total hypertension: adjusted prevalence and variance (delta method)
P.t_ <- P.u_ + P.d
V.t_ <- V.d + V.u_ - (2*V.d/D)

#New hypertension estimates (with 95% CI)
#Undiagnosed
uht_ <- data.frame("prop" = (100*P.u_) %>% round(1), "var" = V.u_) %>%
  mutate("IC95" = (100*qnorm(0.975)*sqrt(var)) %>% round(2),
         "lIC95" = prop - IC95, "uIC95" = prop + IC95,
         "cluster" = "Undiagnosed Hypertension",
         "Year" = c(2000, 2006, 2012)) %>%
  select(-var) %>% rbind(uht[4:10,])

#Total
tht_ <- data.frame("prop" = (100*P.t_) %>% round(1), "var" = V.t_) %>%
  mutate("IC95" = (100*qnorm(0.975)*sqrt(var)) %>% round(2),
         "lIC95" = prop - IC95, "uIC95" = prop + IC95,
         "cluster" = "Total Hypertension",
         "Year" = c(2000, 2006, 2012)) %>%
  select(-var) %>% rbind(tht[4:10,])


#------------------------------------------------------------ ISH/IDH/SDH -----#
# ISH, IDH and SDH hypertension data sets (all years)
ish_w <- ISH_prev_year; ish <- ISH_prev_year_overall
idh_w <- IDH_prev_year; idh <- IDH_prev_year_overall
sdh_w <- SDH_prev_year; sdh <- SDH_prev_year_overall

# Prevalence estimates for 2000, 2006, 2012
P.ish_w <- ish_w[1:3,1]/100; P.ish <- ish[1:3,1]/100
P.idh_w <- idh_w[1:3,1]/100; P.idh <- idh[1:3,1]/100
P.sdh_w <- sdh_w[1:3,1]/100; P.sdh <- sdh[1:3,1]/100

#Variance estimates for 2000, 2006, 2012
V.ish_w <- (ish_w[1:3,2]/qnorm(.975)/100)**2
V.idh_w <- (idh_w[1:3,2]/qnorm(.975)/100)**2
V.sdh_w <- (sdh_w[1:3,2]/qnorm(.975)/100)**2

# New overall prevalences (adjusted UDH * group proportion)
P.ish_ <- P.u_ * P.ish_w
P.idh_ <- P.u_ * P.idh_w
P.sdh_ <- P.u_ * P.sdh_w

# Adjusted variance
# V(XY) = Y^2*Var(X) + X^2*Var(Y) + 2XY*Cov(X,Y)
# Where:
# - X = Adjusted UDH prevalence (P.u_)
# - Y = Proportion of each subgroup (P.ish_w, P.idh_w, P.sdh_w)
# - Var(X) = V.u_
# - Var(Y) = V.ish_w, V.idh_w, V.sdh_w
# - Assuming Cov(X,Y) = 0
V.ish_ <- ((P.ish_w^2)*V.u_) + ((P.u_^2)*V.ish_w) + 0
V.idh_ <- ((P.idh_w^2)*V.u_) + ((P.u_^2)*V.idh_w) + 0
V.sdh_ <- ((P.sdh_w^2)*V.u_) + ((P.u_^2)*V.sdh_w) + 0

#New hypertension estimates (with 95% CI)
#ISH
ish_ <- data.frame("prop" = (100*P.ish_) %>% round(1), "var" = V.ish_) %>%
  mutate("IC95" = (100*qnorm(0.975)*sqrt(var)) %>% round(2),
         "lIC95" = prop - IC95, "uIC95" = prop + IC95,
         "cluster" = "Isolated Systolic Hypertension",
         "Year" = c(2000, 2006, 2012)) %>%
  select(-var) %>% rbind(ish[4:10,])

#IDH
idh_ <- data.frame("prop" = (100*P.idh_) %>% round(1), "var" = V.idh_) %>%
  mutate("IC95" = (100*qnorm(0.975)*sqrt(var)) %>% round(2),
         "lIC95" = prop - IC95, "uIC95" = prop + IC95,
         "cluster" = "Isolated Diastolic Hypertension",
         "Year" = c(2000, 2006, 2012)) %>%
  select(-var) %>% rbind(idh[4:10,])

#SDH
sdh_ <- data.frame("prop" = (100*P.sdh_) %>% round(1), "var" = V.sdh_) %>%
  mutate("IC95" = (100*qnorm(0.975)*sqrt(var)) %>% round(2),
         "lIC95" = prop - IC95, "uIC95" = prop + IC95,
         "cluster" = "Systodiastolic Hypertension",
         "Year" = c(2000, 2006, 2012)) %>%
  select(-var) %>% rbind(sdh[4:10,])


#### WHO age standardization---------------------------- ####

#cdn.who.int/media/docs/default-source/gho-documents/global-health-estimates/
#gpe_discussion_paper_series_paper31_2001_age_standardization_rates.pdf
who_ref <- data.frame(
  strata_as = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
    "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
    "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100+"),
  prop_ref = c(
    8.86, 8.69, 8.60, 8.47, 8.22, 7.93, 7.61, 7.15, 6.59, 6.04,
    5.37, 4.55, 3.72, 2.96, 2.21, 1.52, 0.91, 0.44, 0.15, 0.04, 0.005)) %>%
  filter(!strata_as%in%c( "0-4", "5-9", "10-14", "15-19")) %>%
  transmute(
    prop_ref=prop_ref/100,
    strata_as = case_when(
      strata_as %in% c("20-24", "25-29", "30-34", "35-39") ~ "20-39",
      strata_as %in% c("40-44", "45-49", "50-54", "55-59") ~ "40-59",
      strata_as %in% c("60-64", "65-69", "70-74", "75-79", "80-84",
                       "85-89", "90-94", "95-99", "100+") ~ "60+")) %>%
  group_by(strata_as) %>%
  summarise(prop_ref = sum(prop_ref, na.rm = TRUE), .groups = "drop")

who_ref$strata_as[who_ref$strata_as=="60+"] <- "≥60"


who_prev_HTN_fin <- data.frame(); for(i in 1:length(svy_list)){
  prev <- svyby(~HTN_fin, ~Age_cat+Sex2, svy_list[[i]],
                svymean, vartype = "se", na.rm = TRUE) %>%
    `rownames<-`(NULL) %>% rename("strata_as"=Age_cat)
  prev_m <- prev %>% filter(Sex2=="Men")
  prev_w <- prev %>% filter(Sex2=="Women")
  
  merged <- merge(prev_m, who_ref, by="strata_as")
  prev_std <- 100*sum(merged[3] * merged$prop_ref, na.rm = TRUE)
  var_std  <- sum((merged$se^2) * (merged$prop_ref^2), na.rm = TRUE)
  se_std   <- 100*sqrt(var_std)
  results_m <- data.frame(
    cycle=cycle_names[i], sex="Men", prev_std=prev_std, se_std=se_std,
    lower = prev_std - 1.96 * se_std, upper = prev_std + 1.96 * se_std)
  
  merged <- merge(prev_w, who_ref, by="strata_as")
  prev_std <- 100*sum(merged[3] * merged$prop_ref, na.rm = TRUE)
  var_std  <- sum((merged$se^2) * (merged$prop_ref^2), na.rm = TRUE)
  se_std   <- 100*sqrt(var_std)
  results_w <- data.frame(
    cycle=cycle_names[i], sex="Women", prev_std=prev_std, se_std=se_std,
    lower = prev_std - 1.96 * se_std, upper = prev_std + 1.96 * se_std)
  
  who_prev_HTN_fin <- rbind(who_prev_HTN_fin, results_m, results_w)
}

who_prev_HTN_fin %>% filter(sex=="Men")
who_prev_HTN_fin %>% filter(sex=="Women")



#### Abstract ACC--------------------------------------- ####
#fig.ACC <- ggarrange(
#  labels=c("","","C"), nrow=3, heights = c(1,0.03, 0.75),
#  ggarrange(fig1A, fig1B, ncol=2, labels=LETTERS[1:2],
#            widths=c(1,.75)),"",fig3)
#ggsave(fig.ACC, file="Hypertension phenotypes/Figures/Figure_ACC.png",
#       bg="white", width=55, height=40, units="cm", dpi=300, limitsize=F)

#### Descriptive---------------------------------------- ####

#Year 2000
#xt_2000<-ens00_fin_survey %>% tbl_svysummary(label = list(
#  Sex ~ "Sex", Age_cat ~ "Age group", BMI_cat ~ "Obesity diagnosis",
#  Edu ~ "Education level", SocSec2 ~ "Social security",
#  Indigenous2 ~ "Indigenous identity"), statistic = list(
#    all_categorical() ~ "{p}", ll_continuous() ~ "{mean}"), include = c(
#      Sex, Age_cat, BMI_cat, SBP, DBP, Diabetes2, Edu, SocSec2,
#      Indigenous2), missing="no") %>% add_ci(pattern="{stat} ({ci})") 

#Save table
#tab2 <-tbl_merge(tbls = list(t_2016, t_2018, t_2021, t_2022, t_2023))%>%
#  as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
#doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%
#  body_end_section_landscape() %>% print(target = "Tables/tablaS1.docx")
#remove(t_2016, t_2018, t_2021, t_2022, t_2023, tab2)

