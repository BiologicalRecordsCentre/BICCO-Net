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


###########################################################################

###### converting local climate into difference from national averages ######

climate<-fakeweather
climate1<-fakeweather[,c(13:2,25:14)]


# setting up site-specific climate data in required format:
# each row contains temperature and precipitation from year of survey to five years before the survey
climate2<-cbind(Year[7:52],
                climate1[6:51,1:12],climate1[5:50,1:12],climate1[4:49,1:12],
                climate1[3:48,1:12],climate1[2:47,1:12],climate1[1:46,1:12],
                climate1[6:51,13:24],climate1[5:50,13:24],climate1[4:49,13:24],
                climate1[3:48,13:24],climate1[2:47,13:24],climate1[1:46,13:24])                
climate3<-matrix(0,46,144)


# transforming site-specific climate to difference from UK mean (scaled  by variance)
for(i in 1:46){ for(j in 1:144){
  climate3[i,j]<-(climate2[i,j+1]-transdata[1,j+1])/transdata[2,j+1]}}


### getting species response to each month's weather then averaging for annual effect:
weatherCovs<-matrix(0,46,144)
for(i in 1:46){ for(j in 1:144){
  weatherCovs[i,j]<-as.matrix(climate3[i,j]*spCovariates[sp,j+1])}}
year<-c(1966:2011)

temp<-rowMeans(weatherCovs[,(13-spCovariates[sp,146]):(72-spCovariates[sp,146])])
precipitation<-rowMeans(weatherCovs[,((13-spCovariates[sp,146]):(72-spCovariates[sp,146])+72)])


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
