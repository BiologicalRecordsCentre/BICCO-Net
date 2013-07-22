###### BICCO-Net fine scale analysis ######
###### Gary Powney & Nick Isaac, July 2013 #######

##############################################################################
#
# This is a template for the fine-scale analysis of population dynamics in response to climate change
#
##############################################################################

rm(list=ls())    # clear R

############################################################################## FIND DATA

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- '../Data'

# Option B: point to the Dropbox folder. The user should uncomment the appropriate line
#datadir <- 'C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net' #Gary
#datadir <- 'C:/Documents and Settings/NJBI/My Documents/Dropbox/work/BICCO-Net data' # Nick work
#datadir <- NULL # Blaise to add

############################################################################## LOAD DATA

sp8 <- read.csv(paste0(datadir,"/species_8.csv"),header=T) 					# add ringlet example
sp58 <- read.csv(paste0(datadir,"/species_58.csv"),header=T)					# add Silver-spotter Skipper example
fake_fourier <- read.csv(paste0(datadir,"/fourier_example.csv"),header=T)		# add in fake fourier data
#the initial table will be coefficients, not net effect. To calculate (see below)


### Join fourier scores onto the species index data 
sp8 <- merge(sp8,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))			# for species 8
sp58 <- merge(sp58,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))		# for species 58

############# missing step
#### caulculate LOCAL Temp and Precipitation effect from:
# 1. species-specific fourier covariates 
# 2. local P & T
############# 

###  Mixed effects models 
library(lme4)

# perhaps replace year effect with 1st order autoregressive. lme? gls? can we do gls with poisson?
# fixed effects to add:
#1. local temperature anomoly and it;'s interaction with TEMP
#2. precip, as above
#3. habitat term(s) 
#4. habitat:temp and habitat:precip terms

sp8_mod <- glmer(SINDEX ~ PREV_SINDEX + TEMP + PREC + (1|YEAR) + (1|SITE), data = sp8, family = poisson)	# Species 8 model
summary(sp8_mod)

sp58_mod <- glmer(SINDEX ~ PREV_SINDEX + TEMP + PREC + (1|YEAR) + (1|SITE), data = sp58, family = poisson)	# Species 58 model
summary(sp58_mod)
