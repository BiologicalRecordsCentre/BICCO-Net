###### BICCO-Net fine scale analysis ######
###### Gary Powney & Nick Isaac, July 2013 #######

##############################################################################
#
# Adding the climate results from the fourier analysis into the species data for use in the fine-scale analysis.
#
##############################################################################

rm(list=ls())    # clear R

############################################################################## FIND DATA

library("reshape2")

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- '../Analysis/Data'

# Option B: point to the Dropbox folder. The user should uncomment the appropriate line
#datadir <- 'C:/Users/gawn/Dropbox/BICCO-Net' #Gary
#datadir <- 'C:/Documents and Settings/NJBI/My Documents/Dropbox/work/BICCO-Net data' # Nick work
#datadir <- 'C:/Users/blaise/Dropbox/BICCO-Net' # Blaise to add

############################################################################## LOAD DATA
spp_data <- read.table(paste0(datadir,"/cleaned_UKBMS_site_indices.txt"),header=T,sep=",") # add species site data
spp_data <- spp_data[,c(1,4:5,6:8)]

# add climate data
# 'transdata' is the final version
# it is the national monthly mean and stdev of temp and prec from 1966 - 2011, going back to 1961 - 2006
# spCovariates: once we have the final version fron BioSS this will be final
transdata <- read.table(paste0(datadir,"/transdata.txt"),header=T)
spCovariates <- read.table(paste0(datadir,"/spCovariates"),header=T)

# NEED TO ADD THE BTO SITE CLIMATE DATA HERE IUNSTEAD OF UKBMS
site_climate <- read.table(paste0(datadir,"/ukbms_site_CIP_climate.txt"),header=T,sep="\t") ### Add the UKBMS site/ year climate data

names(site_climate)[3] <- "site"
names(site_climate)[6] <- "m_temp" # change problematic use of mean
names(site_climate)[4] <- "year" 
names(site_climate)[5] <- "month" 
names(site_climate)[7] <- "rainfall"

site_climate <- site_climate[site_climate$year<2012&site_climate$year>1965,] # ensure only including site climate years of interest

fakeweather<-read.table(paste0(datadir,"/fake_weather.txt"),header=T)  # just used for naming columns.

# drop sites from spp_data that do not have climate data #
nrow(spp_data) # 319571
spp_data <- spp_data[spp_data$SITE%in%site_climate$site,]
nrow(spp_data) # 312885

############################################################################## GENERAL STATS
length(unique(spp_data$COMMON_NAME)) # 60 species (only 50 in blaise's butterfly species list) # check which spp lost
length(unique(spp_data$SITE)) # 1724 sites

#################   calculate LOCAL Temp and Precipitation effects:   #######################
spp_list <- read.csv(paste0(datadir,"/Fourier species codes for butterflies.csv"),header=T) # spp code and name information

climate_impact <- data.frame(year=as.numeric(),temp=as.numeric(),precipitation=as.numeric(),SITE=as.numeric(),SPECIES=as.character())  # Create a dataframe to be filled

# test #
#sp <- as.character(spp_list[1,1])

## sort out the name differences ##
levels(spp_list$NAME)[28]<-"Orange-tip"
levels(spp_list$NAME)[5]<-"Chalkhill Blue"
levels(spp_list$NAME)[11]<-"Duke of Burgundy "
levels(spp_list$NAME)[13]<-"Gatekeeper"

#### MISSING SPECIES ####
# "Small/Essex Skipper"  "Black Hairstreak"    "Lulworth Skipper"     "Large Heath"  "Mountain Ringlet"     "Cryptic Wood White"   "Chequered Skipper"    "Glanville Fritillary"

######  For each species in spp_data we want to loop through each site where it has been recorded and use the climate data from Tom O to create the site/species/year climate effect. The table we create will include years where the species was not recorded but the excess years will be dropped when we merge the new cliamte impact data frame with the spp_data. ######
for (sp in as.character(spp_list[,"SPECIES"])){   # loop through each species
  # loop through sites
  site_list <- unique(spp_data[spp_data$COMMON_NAME==as.character(spp_list[spp_list$SPECIES==sp,"NAME"]),"SITE"])  # create a vector of the sites where species "sp" was found.
  cat(sp,"\n")
  
  for (si in site_list) {
    cat(si,"\n")
    test <- site_climate[site_climate$site==si,c("year","month","m_temp","rainfall")] # ID the climate of the site in question
    temper <- dcast(test,year~month,value.var="m_temp")  # convert the site climate data to match the format needed by the code below
    preci <- dcast(test,year~month,value.var="rainfall")
    climate<-cbind(temper,preci[,-1])
    names(climate)<-names(fakeweather) # sort out the names
    
    ### converting local climate into difference from national averages:
    climate1<-climate[,c(13:2,25:14)] # ensure this matches the format of "fakeweather"
    climate2<-cbind(climate$Year[7:41],
                climate1[6:40,1:12],climate1[5:39,1:12],climate1[4:38,1:12],
                climate1[3:37,1:12],climate1[2:36,1:12],climate1[1:35,1:12],
                climate1[6:40,13:24],climate1[5:39,13:24],climate1[4:38,13:24],
                climate1[3:37,13:24],climate1[2:36,13:24],climate1[1:35,13:24])         # create a dataframe that has rows as years and columns of climate climate -1 year, climate -2 years... etc to -5       
    names(climate2)[1] <- "Year"
    climate3<-matrix(0,35,144)
    for(i in 1:35){ for(j in 1:144){
      climate3[i,j]<-(climate2[i,j+1]-transdata[1,j+1])/transdata[2,j+1]}}                  # transdata are the monthly national means and SDs

    ### getting species response to each month's weather then averaging for annual effect:
    weatherCovs<-matrix(0,35,144)
    for(i in 1:35){ for(j in 1:144){
      weatherCovs[i,j]<-as.matrix(climate3[i,j]*spCovariates[spCovariates$species==sp,j+1])}}
    year<-c(1977:2011)

    temp<-rowMeans(weatherCovs[,(13-spCovariates[spCovariates$species==sp,146]):(72-spCovariates[spCovariates$species==sp,146])])
    precipitation<-rowMeans(weatherCovs[,((13-spCovariates[spCovariates$species==sp,146]):(72-spCovariates[spCovariates$species==sp,146])+72)])

    annualweather<-as.data.frame(cbind(c(1977:2011),temp,precipitation))
    names(annualweather) <- c("year","temp","precipitation")
    annualweather$SITE <- rep(si,nrow(annualweather))
    annualweather$SPECIES <- rep(sp,nrow(annualweather))

    climate_impact <- rbind(climate_impact,annualweather)    # build up a dataframe of annual impact of climate on species across all sites and species
  }   
}   

head(climate_impact)    # this file can then be matched to the spp_data for the final fine-scale models to be run

# Add the species name #
final_data <- merge(climate_impact,spp_list)

# merge this with the spp_data file #
names(final_data)[6]<-"COMMON_NAME"
names(final_data)[2]<-"YEAR"
names(final_data)[1]<-"bicconet_code"

final_table <- merge(spp_data,final_data)

### save the file and use this in the final fine-scale analysis ###
write.csv(final_table,file=paste0(datadir,"/UKBMS_spp_data_for_bicco_net_analysis.csv"),row.names=FALSE)

