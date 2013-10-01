rm(list=ls())    # clear R

############################################################################## FIND DATA

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- '../Data'

# Option B: point to the Dropbox folder. The user should uncomment the appropriate line
#datadir <- 'C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net' #Gary
#datadir <- 'C:/Documents and Settings/NJBI/My Documents/Dropbox/work/BICCO-Net data' # Nick work
#datadir <- NULL # Blaise to add

############################################################################## LOAD DATA
recs <- read.table(paste0(datadir,"/UKBMS site indices.txt"),sep=",",header=T) ## add in UKBMS site indices
head(recs)

### Data exclusion ###
# remove sites > 10000 as these are not normal transects (i.e. larval web counts, egg counts, etc...) 
recs <- recs[recs$SITE<10000,]

# drop pilot years
recs <- recs[recs$YEAR>=1976,]

# remove records with -1 
recs <- recs[recs$SINDEX>-1,]

# ensure there are no duplicate rows
recs <- unique(recs)


### Loop through each species and only include sites with at least 7 years of data for ### ~~~~ DISCUSS WITH NICK ET AL. what are the best exclusion criteria ~~~~
spp <- as.character(unique(recs$COMMON_NAME))  # ID species list

new_recs <- NULL # create the empty new_records object to be filled
######TEST####### spp<-spp[1:4]

for (i in spp){
  spp_table <- recs[recs$COMMON_NAME==i,] # create mini table for species of interest

  ### drop sites with less than seven years of data
  sites <- unique(spp_table$SITE)		# ID all sites 
  
  # loop through each site to identify the number of years of data
  counts <- NULL
  for (i in sites) {
    counts <- c(counts,nrow(spp_table[spp_table$SITE==i,])) # ID the number of counts for the species
  }
  
  site_counts <- as.data.frame(cbind(sites,counts))
  site_counts <- site_counts[site_counts$counts>6,]  			# only include sites with 7 or more years of counts
  
  spp_table <- spp_table[spp_table$SITE%in%site_counts$sites,]  # drop sites with less than 7 years of count data
  
  ### identify and add in the site index of the previous year in the main table
  sites <- unique(spp_table$SITE)								# ID all sites
  
  pre_sindex<-NULL											# create pre_sindex variable to be filled
  
  for (i in sites) {			# loop through each site
    temp <- spp_table[spp_table$SITE==i,]	# create a mini table for the site of interest
    temp_pre_sindex <- NULL					# create a temp_pre_sindex variable to be filled
    for (j in temp$YEAR) {				# loop through each year that has a count for the site in question 
      if(nrow(temp[temp$YEAR==j-1,]) < 1) {				# if no count then temp_pre_sindex = NA
        temp_pre_sindex <- c(temp_pre_sindex,NA)				
      } else {											# else add the sindex of the previous year
        temp_pre_sindex <- c(temp_pre_sindex,round(mean(temp[temp$YEAR==j-1,"SINDEX"])))  # we say mean here as some sites have more than 1 SINDEX for each year.
      }	
    }
    pre_sindex <- c(pre_sindex,temp_pre_sindex)	  # build up a list of site index values for the previous year
  }
  
  spp_table$pre_sindex <- pre_sindex		# add the previous years site index scores on the main dataframe
  new_recs<-rbind(new_recs,spp_table)		# re-build the records file but now excluding sites with less than 7 years of data and also adding in the sindex of the year before.
  
} # close the species loop

head(new_recs)


write.csv(new_recs,file=paste0(datadir,"/cleaned_UKBMS_site_indices.txt"),row.names=F)
