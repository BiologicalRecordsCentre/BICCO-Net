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
datadir <- "C:/Users/blaise/Dropbox/BICCO-Net/Climate transformation" # Blaise

############################################################################## LOAD DATA

sp8 <- read.csv(paste0(datadir,"/species_8.csv"),header=T) 					# add ringlet example
sp58 <- read.csv(paste0(datadir,"/species_58.csv"),header=T)					# add Silver-spotter Skipper example
fake_fourier <- read.csv(paste0(datadir,"/fourier_example.csv"),header=T)		# add in fake fourier data
#the initial table will be coefficients, not net effect. To calculate (see below)


### Join fourier scores onto the species index data 
sp8 <- merge(sp8,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))			# for species 8
sp58 <- merge(sp58,fake_fourier,by.x=c("SITE","YEAR",by.y=c("SITE","YEAR")))		# for species 58


###########################################################################

###### converting local climate into difference from national averages ######

climate<-read.table(paste0(datadir,"/fake_weather.txt"),header=T)
transdata<-read.table(paste0(datadir,"/transdata.txt"),header=T)
spCovariates<-read.table(paste0(datadir,"/spCovariates.txt"),header=T)

climate1<-climate[,c(13:2,25:14)]
attach(climate)

# setting up site-specific climate data in required format:
# each row contains temperature and precipitation from year of survey to ten years before the survey
# (it looks like 11 years are included because the 10 included years do not correspond to calendar years)
ystart<-c(1:length(Year))[Year==1966]
yend<-c(1:length(Year))[Year==2011]

# i = types of weather i.e. mean temp and precipitation
# j = years back i.e. we look at weather for up to 10 years before the survey period
climate2<-matrix(0,46,264)
for(i in 1:2){for(j in 1:11){
  weathertype=cbind(1:12,13:24)
climate2[,((12*(j-1)+1):(12*j))+(132*(i-1))]<-as.matrix(climate1[(ystart+1-j):(yend+1-j),weathertype[,i]])
}}

climate3<-matrix(0,46,264)

# transforming site-specific climate to difference from UK mean (scaled  by variance)
for(i in 1:46){ for(j in 1:264){
  climate3[i,j]<-(climate2[i,j]-transdata[1,j+1])/transdata[2,j+1]}}


### getting species response to each month's weather then averaging for annual effect:
#sp<-1 # here select which species we are interested in
weatherCovs<-matrix(0,46,264)
for(i in 1:46){ for(j in 1:264){
  weatherCovs[i,j]<-as.matrix(climate3[i,j]*spCovariates[sp,j+1])}}

temp<-rowSums(weatherCovs[,1:132])/120
precipitation<-rowSums(weatherCovs[,133:264])/120

# annualweather variables (temp & precipitation) are an annual estimation of 
# population growth in species[sp] due to temperature and rainfall.
annualweather<-as.data.frame(cbind(c(1966:2011),temp,precipitation))
names(annualweather)<-c("year","temp","precipitation")
attach(annualweather)


# The interaction between these variables and habitat variable(s) should show
# whether habitat affects species' response to weather.

###########################################################################
### Missing step NEW:

### 'temp' and 'precipitation' give the annual effect for species and site 
### selected in lines 44:47 from 1966 - 2011. 
### A loop will need to be added to produce these values for each site

### select years required for glmers:
### will need to be change to do it for each site
TEMP1<-temp[match(sp8$YEAR[sp8$SITE==site],year)]
PREC1<-precipitation[match(sp8$YEAR[sp8$SITE==site],year)]

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
