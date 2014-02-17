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
# datadir <- '../Analysis/Data'

 datadir <- "C:\\Users\\gawn\\Documents\\Gary Work\\BRC Quantitative Ecologist\\BICCO\\Analysis\\data"  # GP working from home
# Option B: point to the Dropbox folder. The user should uncomment the appropriate line
# datadir <- 'C:/Users/gawn/Dropbox/BICCO-Net' #Gary
# datadir <- 'C:/Documents and Settings/NJBI/My Documents/Dropbox/work/BICCO-Net data' # Nick work
# datadir <- 'C:/Users/blaise/Dropbox/BICCO-Net' # Blaise to add

############################################################################## LOAD DATA
#spp_data <- read.table(paste0(datadir,"/cleaned_UKBMS_site_indices.txt"),header=T,sep=",") # add species site data
spp_data <- read.csv(paste0(datadir,"/CBC_for_plots_with 10_years_species_removed.csv"),header=T) # add species site data
#spp_data <- spp_data[,c(1,4:5,6:8)]
spp_data <- spp_data[,c("Gridref","Year","Spname","Count")]

# add climate data
# 'transdata' is the final version
# it is the national monthly mean and stdev of temp and prec from 1966 - 2011, going back to 1961 - 2006
# spCovariates: once we have the final version fron BioSS this will be final
transdata <- read.table(paste0(datadir,"/transdata.txt"),header=T)
spCovariates <- read.table(paste0(datadir,"/spCovariates"),header=T)

# site_climate <- read.table(paste0(datadir,"/ukbms_site_CIP_climate.txt"),header=T,sep="\t") ### Add the UKBMS site/ year climate data
site_climate <- read.table(paste0(datadir,"/cbc_site_CIP_climate.txt"),header=T,sep="\t") ### Add the CBC site/ year climate data

names(site_climate)[3] <- "site"
names(site_climate)[6] <- "m_temp" # change problematic use of mean
names(site_climate)[4] <- "year" 
names(site_climate)[5] <- "month" 
names(site_climate)[7] <- "rainfall"

# site_climate <- site_climate[site_climate$year<2012&site_climate$year>1965,] # ensure only including site climate years of interest #UKBMS
site_climate <- site_climate[site_climate$year<2003&site_climate$year>1951,] # ensure only including site climate years of interest #CBC

fakeweather<-read.table(paste0(datadir,"/fake_weather.txt"),header=T)  # just used for naming columns.

### sorting CBC column names ###
names(spp_data)[1] <- "SITE" 
names(spp_data)[3] <- "COMMON_NAME" 

# drop sites from spp_data that do not have climate data #
nrow(spp_data) # 319571
spp_data <- spp_data[spp_data$SITE%in%site_climate$site,]
nrow(spp_data) # 312885

############################################################################## GENERAL STATS
length(unique(spp_data$COMMON_NAME)) # UKBMS = 60 species (only 50 in blaise's butterfly species list) # check which spp lost
length(unique(spp_data$COMMON_NAME)) # CBC = 84 species () # check which spp lost

length(unique(spp_data$SITE)) # 1724 UKBMS sites, 323 CBC sites

#################   calculate LOCAL Temp and Precipitation effects:   #######################
spp_list <- read.csv(paste0(datadir,"\\Species_names.csv"),header=T)
#spp_list <- read.csv(paste0(datadir,"/Fourier species codes for butterflies.csv"),header=T) # spp code and name information
spp_list <- spp_list[spp_list$Taxa=="Bird",c("name","COMMON")]# extract the birds
names(spp_list)[1]<-"SPECIES"
names(spp_list)[2]<-"NAME"

climate_impact <- data.frame(year=as.numeric(),temp=as.numeric(),precipitation=as.numeric(),SITE=as.numeric(),SPECIES=as.character())  # Create a dataframe to be filled

# test #
#sp <- as.character(spp_list[1,1])

## sort out the name differences - UKBMS specific ##
#levels(spp_list$NAME)[28]<-"Orange-tip"
#levels(spp_list$NAME)[5]<-"Chalkhill Blue"
#levels(spp_list$NAME)[11]<-"Duke of Burgundy "
#levels(spp_list$NAME)[13]<-"Gatekeeper"

#### MISSING SPECIES - UKBMS SPECIFIC ####
# "Small/Essex Skipper"  "Black Hairstreak"    "Lulworth Skipper"     "Large Heath"  "Mountain Ringlet"     "Cryptic Wood White"   "Chequered Skipper"    "Glanville Fritillary"

######  For each species in spp_data we want to loop through each site where it has been recorded and use the climate data from Tom O to create the site/species/year climate effect. The table we create will include years where the species was not recorded but the excess years will be dropped when we merge the new cliamte impact data frame with the spp_data. ######

### TEMP REMOVE WHEN RUNNING FINAL MODELS ###

for (sp in as.character(spp_list[,"SPECIES"])){   # loop through each species
  # loop through sites
  site_list <- as.character(unique(spp_data[spp_data$COMMON_NAME==as.character(spp_list[spp_list$SPECIES==sp,"NAME"]),"SITE"]))  # create a vector of the sites where species "sp" was found.
  cat(sp,"\n")
  
  # si<-site_list[1] # test
  
  for (si in site_list) {
    cat(si,"\n")
    test <- site_climate[site_climate$site==si,c("year","month","m_temp","rainfall")] # ID the climate of the site in question
    temper <- dcast(test,year~month,value.var="m_temp")  # convert the site climate data to match the format needed by the code below
    preci <- dcast(test,year~month,value.var="rainfall")
    climate<-cbind(temper,preci[,-1])
    names(climate)<-names(fakeweather) # sort out the names
    
    ### converting local climate into difference from national averages:
    climate1<-climate[,c(13:2,25:14)] # ensure this matches the format of "fakeweather"
    climate2<-cbind(climate$Year[12:46],
                    climate1[11:45,1:12],climate1[10:44,1:12],climate1[9:43,1:12],
                    climate1[8:42,1:12],climate1[7:41,1:12],climate1[6:40,1:12],
                    climate1[5:39,1:12],climate1[4:38,1:12],climate1[3:37,1:12],
                    climate1[2:36,1:12],climate1[1:35,1:12],
                    
                    climate1[11:45,13:24],climate1[10:44,13:24],climate1[9:43,13:24],
                    climate1[8:42,13:24],climate1[7:41,13:24],climate1[6:40,13:24],
                    climate1[5:39,13:24],climate1[4:38,13:24],climate1[3:37,13:24],
                    climate1[2:36,13:24],climate1[1:35,13:24])        # create a dataframe that has rows as years and columns of climate climate -1 year, climate -2 years... etc to -5 # UPDATE # to -10       
    
    names(climate2)[1] <- "Year"
    
    climate3<-matrix(0,46,264)
    
    for(i in 1:46){ for(j in 1:264){
      climate3[i,j] <- (climate2[i,j+1]-transdata[1,j+1])/transdata[2,j+1]}}                  # transdata are the monthly national means and SDs
    # climate3[i,j]<-(climate2[i,j+1]-transdata[1,j+1])/transdata[2,j+1]}}
    ### getting species response to each month's weather then averaging for annual effect:
    weatherCovs<-matrix(0,46,264)
    
    for(i in 1:46){ for(j in 1:264){
      weatherCovs[i,j]<-as.matrix(climate3[i,j]*spCovariates[spCovariates$species==sp,j+1])}}
    
    year<-c(1966:2011)

    ### sum across all years
    temp<-rowSums(weatherCovs[,1:132])/120
    precipitation<-rowSums(weatherCovs[,133:264])/120

    annualweather<-as.data.frame(cbind(c(1966:2011),temp,precipitation))
    names(annualweather)<-c("year","temp","precipitation")
    annualweather$SITE <- rep(si,nrow(annualweather))
    annualweather$SPECIES <- rep(sp,nrow(annualweather))

    climate_impact <- rbind(climate_impact,annualweather)    # build up a dataframe of annual impact of climate on species across all sites and species
  }   
}   

head(climate_impact)    # this file can then be matched to the spp_data for the final fine-scale models to be run


### LOAD IN THE CLUSTER OUTPUT ###
climate_impact <- read.table("C:\\Users\\gawn\\Documents\\Gary Work\\BRC Quantitative Ecologist\\BICCO\\Analysis\\data\\CBC_BICCO_cluster_output.txt",header=T)

# Add the species name #
final_data <- merge(climate_impact,spp_list)

# merge this with the spp_data file #
# rename matching columns for the merge #
names(final_data)[6]<-"COMMON_NAME"
names(final_data)[2]<-"YEAR"
names(final_data)[1]<-"bicconet_code"
names(spp_data)[2] <- "YEAR"

final_table <- merge(spp_data,final_data)

### save the file and use this in the final fine-scale analysis ###
write.csv(final_table,file=paste0(datadir,"/CBC_spp_data_for_bicco_net_analysis.csv"),row.names=FALSE)

