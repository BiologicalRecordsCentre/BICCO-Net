###### BICCO-Net ID site climate anomaly  ######
###### Gary Powney & Nick Isaac, Dec 2013 #######

##############################################################################
#
# This is for use in the fine-scale analysis 
#
##############################################################################

rm(list=ls())    # clear R

############################################################################## FIND DATA

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- 'W:/PYWELL_SHARED/Pywell Projects/BRC/BICCO-Net/Phase 2/Analysis/Data'

############################################################################## LOAD DATA
site_climate <- read.table(paste0(datadir,"/all.sites.CIP.mean.temp.AND.rainfall.1971_2012_running.3month.means.txt"),header=T,sep="\t") ### Add the site/ year climate data
site_climate <- site_climate[site_climate$year<2012,] 
names(site_climate)[4] <- "m_temp" # change problematic use of mean

# drop sites below <10000 as these are not proper UKBMS sites #
site_climate <- site_climate[site_climate$site<10000,]

###### Run through each site and ID the yearly anomaly (climate anomaly = climate - mean climate across all sites) for each year ######

### Next, loop through sites to id the climate anomalies. ###
site_list <- unique(site_climate$site)
year_list <- unique(site_climate$year)

TEMP_ANO <- NULL
PREC_ANO <- NULL
SITE <- NULL
YEAR <- NULL

# loop through each site #
for (i in site_list){
  temp_table <- site_climate[site_climate$site==i,]
  
  for (j in year_list){ # loop through each year
    TEMP_ANO <- c(TEMP_ANO, mean(temp_table[temp_table$year==j,"m_temp"]) - mean(site_climate[site_climate$year==j,"m_temp"])) # ID the diff between local and overall annual temp
    PREC_ANO <- c(PREC_ANO, mean(temp_table[temp_table$year==j,"rainfall"]) - mean(site_climate[site_climate$year==j,"rainfall"])) # ID the diff between local and overall annual rainfall)
    SITE <- c(SITE,i)
    YEAR <- c(YEAR,j)
  }  
}

clim_ano <- data.frame(site=SITE,year=YEAR,temp_ano=TEMP_ANO,prec_ano=PREC_ANO)

write.csv(clim_ano,file=paste0(datadir,"/climate_anomaly.csv"),row.names=FALSE)

