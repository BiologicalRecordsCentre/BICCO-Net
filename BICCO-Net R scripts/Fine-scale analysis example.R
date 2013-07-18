###### BICCO-Net fine scale analysis ######
###### G Powney, July 2013 ######

rm(list=ls())  	# clear R

### add data - Direct these to the BICCO-Net dropbox folder 		
sp8 <- read.csv("C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net\\species_8.csv",header=T) 					# add ringlet example
sp58 <- read.csv("C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net\\species_58.csv",header=T)					# add Silver-spotter Skipper example
fake_fourier <- read.csv("C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net\\fourier_example.csv",header=T)		# add in fake fourier data

### Join fourier scores onto the species index data 
sp8 <- merge(sp8,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))			# for species 8
sp58 <- merge(sp58,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))		# for species 58

###  Mixed effects models 
library(lme4)

sp8_mod <- glmer(SINDEX ~ PREV_SINDEX + TEMP + PREC + (1|YEAR) + (1|SITE), data = sp8, family = poisson)	# Species 8 model
summary(sp8_mod)

sp58_mod <- glmer(SINDEX ~ PREV_SINDEX + TEMP + PREC + (1|YEAR) + (1|SITE), data = sp58, family = poisson)	# Species 58 model
summary(sp58_mod)
