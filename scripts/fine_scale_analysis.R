###### BICCO-Net fine scale analysis ######
###### Gary Powney & Nick Isaac, July 2013 #######

##############################################################################
#
# This is the fine-scale analysis of population dynamics in response to climate change
#
##############################################################################

rm(list=ls())    # clear R

############################################################################## FIND DATA

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- 'W:/PYWELL_SHARED/Pywell Projects/BRC/BICCO-Net/Phase 2/Analysis/Data'

############################################################################## LOAD DATA

### load the species/site data ###
spp_data <- read.csv(paste0(datadir,"/UKBMS_spp_data_for_bicco_net_analysis.csv"),header=T)    # add species data
site_data <- read.table(paste0(datadir,"/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt"),header=T,sep="\t")   # add site level habitat data
clim_ano <- read.csv(paste0(datadir,"/climate_anomaly.csv"),header=T)

### Drop species ###
spp_data<-spp_data[!spp_data$COMMON_NAME=="Large Blue",] # Large Blue


############# 
###  Mixed effects models 
############# 

# perhaps replace year effect with 1st order autoregressive. lme? gls? can we do gls with poisson?
# fixed effects to add:
#1. local temperature anomaly and it's interaction with TEMP

### Loop through each species, extract the various model parameters of interest ###
spp_list <- as.character(unique(spp_data$COMMON_NAME))
require(MuMIn)
require(lme4)

# rename temp column
names(spp_data)[8]<-"temper"


### merge climate anomaly data with species data ###
# rename for merge
names(clim_ano)[1] <- "SITE"
names(clim_ano)[2] <- "YEAR"

spp_data <- merge(spp_data,clim_ano)




##########################################
# find-scale models 1: Climate + anomaly #
##########################################

# sindex ~ sindex of previous year + temper + precipitation + temp_ano + prec_ano + temper:temp_ano + temper:prec_ano + precipitation:temp_ano + precipitation:prec_ano + (1|SITE) + (1|YEAR) #  NEED TO ADD TEMP ANOMALY

intercepts<-NULL
b_PRE_SIN<-NULL
b_TEMP<-NULL
b_PREC<-NULL
b_TEMP_ANO<-NULL
b_PREC_ANO<-NULL
b_TEMP_TEMP_ANO_INT<-NULL
b_TEMP_PREC_ANO_INT<-NULL
b_PREC_TEMP_ANO_INT<-NULL
b_PREC_PREC_ANO_INT<-NULL
cond_r_squ<-NULL
intercepts_p<-NULL
b_PRE_SIN_p<-NULL
b_TEMP_p<-NULL
b_PREC_p<-NULL
b_TEMP_ANO_p<-NULL
b_PREC_ANO_p<-NULL
b_TEMP_TEMP_ANO_INT_p<-NULL
b_TEMP_PREC_ANO_INT_p<-NULL
b_PREC_TEMP_ANO_INT_p<-NULL
b_PREC_PREC_ANO_INT_p<-NULL
num_sites<-NULL
num_years<-NULL
mod_AIC <- NULL

# quick code for extracting parameters from the GLMM model outputs #

for (i in spp_list){
  cat(i,"\n")
  temp <- spp_data[spp_data$COMMON_NAME==i,]
  mod_1 <- glmer(SINDEX ~ pre_sindex + temper + precipitation + temp_ano + prec_ano + temper:temp_ano + temper:prec_ano + precipitation:temp_ano + precipitation:prec_ano + (1|SITE) + (1|YEAR), data = temp, family = poisson)  # Species model
  intercepts <- c(intercepts,summary(mod_1)$coefficients[1,1])
  b_PRE_SIN <- c(b_PRE_SIN,summary(mod_1)$coefficients[2,1])
  b_TEMP <- c(b_TEMP,summary(mod_1)$coefficients[3,1])
  b_PREC <- c(b_PREC,summary(mod_1)$coefficients[4,1])
  b_TEMP_ANO <- c(b_TEMP_ANO,summary(mod_1)$coefficients[5,1])
  b_PREC_ANO <- c(b_PREC_ANO,summary(mod_1)$coefficients[6,1])
  b_TEMP_TEMP_ANO_INT <- c(b_TEMP_TEMP_ANO_INT,summary(mod_1)$coefficients[7,1])
  b_TEMP_PREC_ANO_INT <- c(b_TEMP_PREC_ANO_INT,summary(mod_1)$coefficients[8,1])
  b_PREC_TEMP_ANO_INT <- c(b_PREC_TEMP_ANO_INT,summary(mod_1)$coefficients[9,1])
  b_PREC_PREC_ANO_INT <- c(b_PREC_PREC_ANO_INT,summary(mod_1)$coefficients[10,1])
  cond_r_squ <- c(cond_r_squ,r.squaredGLMM(mod_1)[1])
  intercepts_p <- c(intercepts_p,summary(mod_1)$coefficients[1,4])
  b_PRE_SIN_p <- c(b_PRE_SIN_p,summary(mod_1)$coefficients[2,4])
  b_TEMP_p <- c(b_TEMP_p,summary(mod_1)$coefficients[3,4])
  b_PREC_p <- c(b_PREC_p,summary(mod_1)$coefficients[4,4])
  b_TEMP_ANO_p <- c(b_TEMP_ANO_p,summary(mod_1)$coefficients[5,4])
  b_PREC_ANO_p <- c(b_PREC_ANO_p,summary(mod_1)$coefficients[6,4])
  b_TEMP_TEMP_ANO_INT_p <- c(b_TEMP_TEMP_ANO_INT_p,summary(mod_1)$coefficients[7,4])
  b_TEMP_PREC_ANO_INT_p <- c(b_TEMP_PREC_ANO_INT_p,summary(mod_1)$coefficients[8,4])
  b_PREC_TEMP_ANO_INT_p <- c(b_PREC_TEMP_ANO_INT_p,summary(mod_1)$coefficients[9,4])
  b_PREC_PREC_ANO_INT_p <- c(b_PREC_PREC_ANO_INT_p,summary(mod_1)$coefficients[10,4])
  num_sites <- c(num_sites,length(unique(temp$SITE)))
  num_years <- c(num_years,length(unique(temp$YEAR)))
  mod_AIC <- c(mod_AIC,AIC(mod_1))
}

results <- data.frame(COMMON_NAME=spp_list,intercept=intercepts,intercepts_p=intercepts_p,b_PRE_SIN=b_PRE_SIN,b_PRE_SIN_p=b_PRE_SIN_p,b_TEMP=b_TEMP,b_TEMP_p=b_TEMP_p,b_PREC=b_PREC,b_PREC_p=b_PREC_p,b_PREC_ANO=b_PREC_ANO,b_PREC_ANO_p=b_PREC_ANO_p,b_TEMP_ANO=b_TEMP_ANO,b_TEMP_ANO_p=b_TEMP_ANO_p,b_TEMP_TEMP_ANO_INT=b_TEMP_TEMP_ANO_INT,b_TEMP_TEMP_ANO_INT_p=b_TEMP_TEMP_ANO_INT_p,b_TEMP_PREC_ANO_INT=b_TEMP_PREC_ANO_INT,b_TEMP_PREC_ANO_INT_p=b_TEMP_PREC_ANO_INT_p,b_PREC_TEMP_ANO_INT=b_PREC_TEMP_ANO_INT,b_PREC_TEMP_ANO_INT_p=b_PREC_TEMP_ANO_INT_p,b_PREC_PREC_ANO_INT=b_PREC_PREC_ANO_INT,b_PREC_PREC_ANO_INT_p=b_PREC_PREC_ANO_INT_p,fe_r_squ=cond_r_squ,num_sites=num_sites,num_years=num_years,mod_AIC=mod_AIC)



### MAKE THIS PASTE USING THE DIRECTORY AND ALSO ADD THE DATE AND GROUP (BTO/UKBMS) ###
write.csv(results,file="W:/PYWELL_SHARED/Pywell Projects/BRC/Gary/BICCO-NET/Results/GLMM outputs/model_1_results.csv",row.names=FALSE)




############## NEED TO DECIDE WHETHER TO USE MUMIN HERE? #############
# Discuss with Nick:
#   - Do we want to just have the smae full model across all species or simplify the model based on MuMIn?


############## WHAT IS THE BEST WAY TO PRESENT THESE MODELS ACROSS MULTIPLE SPECIES #############
# Use histograms of the coefs for each group and model parameter #
# OR...
# Box plots of the coefs split for each taxa and parameter #







########################################
# model 3 - adding habitat interaction #
########################################

###### THIS NEEDS TO BE UPDATED FOLLOWING THE DECISION ON WHICH HYPOTHESIS WE TEST  ######

#require("reshape2")

## look at habitat in a 2 km buffer ##
#site_2 <- site_data[site_data$buffer==2000&site_data$Surv=="UKBMS",c("siteno.gref","A","BgRo","Br","BW","C","CW","F","G","H","M","S","R","UG")]

## id dominant habitat type for the site ##
#dominant_habitat <- NULL
#for (i in 1:nrow(site_2)){
#  temp_table<-melt(site_2[i,2:14])
#  dominant_habitat <- c(dominant_habitat,as.character(temp_table[temp_table$value==max(temp_table$value),1]))
#}
#site_2$dominant_habitat<-dominant_habitat
# 
# # join dominant habitat info onto the spp_data
# names(site_2)[1] <- "SITE"
# site_2_merge <- site_2[,c(1,15)]
# final_data <- merge(spp_data,site_2_merge)
# 
# intercepts<-NULL
# b_PRE_SIN<-NULL
# b_TEMP<-NULL
# b_PREC<-NULL
# cond_r_squ<-NULL
# intercepts_p<-NULL
# b_PRE_SIN_p<-NULL
# b_TEMP_p<-NULL
# b_PREC_p<-NULL
# num_sites<-NULL
# num_years<-NULL
# 
# # need to update this part #
# 
# for (i in spp_list){
#   cat(i,"\n")
#   temp <- final_data[final_data$COMMON_NAME==i,]
#   mod_1 <- glmer(SINDEX ~ pre_sindex + temper + temper:dominant_habitat + precipitation + precipitation:dominant_habitat + dominant_habitat + (1|SITE), data = temp, family = poisson)  # Species model
#   intercepts <- c(intercepts,summary(mod_1)$coefficients[1,1])
#   b_PRE_SIN <- c(b_PRE_SIN,summary(mod_1)$coefficients[2,1])
#   b_TEMP <- c(b_TEMP,summary(mod_1)$coefficients[3,1])
#   b_PREC <- c(b_PREC,summary(mod_1)$coefficients[4,1])
#   cond_r_squ <- c(cond_r_squ,r.squaredGLMM(mod_1)[1])
#   intercepts_p <- c(intercepts_p,summary(mod_1)$coefficients[1,4])
#   b_PRE_SIN_p <- c(b_PRE_SIN_p,summary(mod_1)$coefficients[2,4])
#   b_TEMP_p <- c(b_TEMP_p,summary(mod_1)$coefficients[3,4])
#   b_PREC_p <- c(b_PREC_p,summary(mod_1)$coefficients[4,4])
#   num_sites <- c(num_sites,length(unique(temp$SITE)))
#   num_years <- c(num_years,length(unique(temp$YEAR))) 
# }
# 
# results <- data.frame(COMMON_NAME=spp_list,intercept=intercepts,intercepts_p=intercepts_p,b_PRE_SIN=b_PRE_SIN,b_PRE_SIN_p=b_PRE_SIN_p,b_TEMP=b_TEMP,b_TEMP_p=b_TEMP_p,b_PREC=b_PREC,b_PREC_p=b_PREC_p,fe_r_squ=cond_r_squ,num_sites=num_sites,num_years=num_years)
# 
