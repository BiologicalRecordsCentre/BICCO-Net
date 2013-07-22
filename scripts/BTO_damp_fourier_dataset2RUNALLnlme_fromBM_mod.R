# Script written by Blaise Martay
# Modified to read data on other machines by nick Isaac, 22 July


#LOOK HERE TO FIND OUT ABOUT AR1 CORRELATIONS,
#APPARENTLY ALLOWING NEGATIVE VALUES
#http://stat.ethz.ch/R-manual/R-devel/library/nlme/html/corAR1.html

#R PROGRAM FOR CONSTRAINING MONTHLY REGRESSION COEFFICIENTS TO LIE ON DAMPED FOURIER CURVES,
#TREATING MODEL HAS HAVING 1 NON-LINEAR PARAMETER, THE REMAINDER BEING LINEAR
#CONDITIONAL ON THE VALUE OF THE EXPONENTIAL DECAY PARAMETER
#AS SENT TO BTO
#

# THIS VERSION RUNS USING NLME

#Delete all existing structures, to ensure nothing gets carrier forwards from a
#previous analysis into the results
rm(list=ls(all=TRUE))

############## CROSS-USER FILE LOCATIONS

# Option A: data are stored locally in a folder on the same level as the folder containing the RStudio project
datadir <- '../Data'

# Option B: point to the Dropbox folder. The user should uncomment the appropriate line
#datadir <- 'C:\\Documents and Settings\\GAWN\\My Documents\\My Dropbox\\BICCO-Net' #Gary
#datadir <- 'C:/Documents and Settings/NJBI/My Documents/Dropbox/work/BICCO-Net data' # Nick work
#datadir <- "C:/Users/blaise/Desktop/BICCO-Net/Fourier bioss 04.04" #Blaise

############## END NI CODE

#
#CALL NECESSARY LIBRARIES
#
library(nlme)
library(MASS)
library(Hmisc)
library(gmodels)

#####
# Declare whether (s_n_c=1) or not (s_n_c=0) to use scale and centre the covariates
#
s_n_c <- 1

#####
# control parameters for...
drawgraphs <- 0   # drawing graphs (0 don't, 1 do)
pause <- "FALSE"  # pausing between graphs "TRUE" do, "FALSE" don't
fitE <- 1         # fitting the E model, with periods 12 months and 24 months
calcQuant <- 1    # calculating quantiles for the regression coefficients
gridmin <- -3     # lower limit for grid search with non-linear parameter
gridmax <- 10     # upper limit for grid search
no_decay_threshold <- 9     # upper theshold for grid search
years_back <- 40            #number of years weather data to include  
try_no_decay_model <- TRUE  # whether to force no decay if likelihood minimum is at gridmax


#
#DATA VARIABLES
###Declare directories to read from and years of first records
birdyr1 <- 1966   # first year with bird data
metyr1 <- 1910    # first year with met data
yrlag <- birdyr1 - metyr1 
nyr <- 46 #number of years 46 for BTO data with the 0s of 1966 in, 100 for michelle's paper
ydatafile <- "BTO_all_indices_all_years.csv"  #read species index data

#READ FILE OF SPECIES NAMES AND ABBREVIATIONS
setwd(datadir)
birds <- read.csv("BTOallspecies.csv",header=TRUE)
birdnames <- birds$speciesname

#NAMES OF MET VARIABLES
#READING OF MET VARIABLES TAKES PLACE WITHIN THE PROGRAM
metnames <- c("minT","maxT","meanT","prec")


#
#LOOP ROUND MET VARIABLES AND SPECIES
#AS SET UP NOW, NO NEED TO CAPTION OUT CLOSING }} IF PICKING OUT A SINGLE SP. FOR INVESTIGATION
# species 56 is an example of two different decay parameters giving similar models
#
#imet_single <- 1
#isp_single <- 56
#for (imet in imet_single:imet_single){
#for (isp in isp_single:isp_single){

for (imet in 3:3){
  
  ## isp in 1:494 won't run all the way through - I had to do it in chunks, 
  ### restarting when a species didn't work.
  ### (sp148 made something get stuck in a loop rather than stopping the program)
  for (isp in c(1:494)){
    
    #FLAGS AS TO HOW THE RESULTS WERE DERIVED 
    # -1 IF DECAY PARAMETER SET TO 1 (NO DECAY)
    # 0 IF FINAL ESTIMATION HAD DECAY PARAMETER CONSTRAINED (AS ESTIMATED FROM lme)
    # 1 IF FINAL ESTIMATION FROM NMLE
    final0 <- 1
    finalC <- 1
    finalD <- 1
    finalE <- 1
    
    metname <- metnames[imet]
    birdname <- birdnames[isp]
    
    savenameE <- paste(birdname,"_",metname,isp,",E.Rsave",sep="")
    savenamenoE <- paste(birdname,"_",metname,isp,",noE.Rsave",sep="")
    #save(imet,isp,file=savenameE)
    #save(imet,isp,file=savenamenoE)
    
    x_files <-c("min_temp","max_temp","mean_temp","precipitation")
    xdatafile <- paste(x_files[imet],"_from_1910.csv",sep="")
    xdatafile
    
    
    #setwd(datadir)
    ydata <- read.csv(ydatafile,header=TRUE)
    xdata <- read.csv(xdatafile,header=TRUE)
    
    head(ydata)
    head(xdata)
    tail(xdata)
    tail(ydata)
    
    
    
    #
    # FROM HERE ON, CODE IS ENTIRELY GENERAL
    # OUTPUTS TO SCREEN ONLY AT PRESENT
    #
    
    twopi <- 2*pi
    x60 <- 1:60
    nq <- 5
    prob_q <- c(0.025,0.25,0.5,0.75,0.975)
    
    # Other controlling variables, set automatically from user specifications above
    
    out.df <- as.data.frame(array(NA,dim=c(ncol(ydata[-1]),44)))
    
    ybird <- ydata[,isp+1]     # pick a response variable
    colmonmat <- 12*years_back  #NO. COLUMNS IN AUGMENTED WEATHER VARIABLE ARRAY
    
    #
    # DECLARE AND CREATE STRUCTURES
    #
    
    #year covariate
    year <- c(1:nyr)
    
    #
    #CREATE AUGMENTED DATA MATRIX
    #1ST COLUMN MONTH OF RECORDING
    #THEREAFTER SUCCESSIVE MONTHS BEFOREHAND
    #
    
    monmat_obs <- array(0,dim=c(nyr,colmonmat))
    monmat <- array(0,dim=c(nyr,colmonmat))
    
    for (i in 1:46){
      for (j in 1:years_back){
        for (k in 1:6 ){
          monmat_obs[i,k+12*(j-1)] <- xdata[yrlag+i-(j-1),8-k]
        }
        for (k in 7:12){
          monmat_obs[i,k+12*(j-1)] <- xdata[yrlag+i-(j-1),20-k]
        }
      }
    }
    
    for (j in 1:colmonmat) {
      if (s_n_c==0) {
        monmat[ ,j] <- monmat_obs[ ,j]
      }
      if (s_n_c==1){
        monmat[ ,j] <- (monmat_obs[ ,j] -mean(monmat_obs[ ,j]))/sqrt(var(monmat_obs[ ,j]))
      }
    }
    
    
    
    head(cbind(ybird,monmat_obs))
    head(xdata)
    
    tail(cbind(ybird,monmat_obs))
    tail(xdata)
    
    
    #
    #MATRIX OF POWERS FOR CORRELATION MATRIX
    #
    lcm <- array (dim=c(nyr,nyr))
    for (i in 1:nyr){
      for (j in 1:nyr){
        lcm[i,j] <- abs(i-j)
      }
    }
    head(lcm)
    tail(lcm)
    
    y <- ybird
    
    #PLOT DATA USED
    if (drawgraphs==1){
      par(ask=pause)
      ts.plot(y)}
    
    
    #REMOVE ROWS AND COLUMNS FOR YEAR WITH MISSING RESPONSE
    
    #IDENTIFY AND COUT MISSING VALUES
    npresy <- sum(complete.cases(y)) # Count of complete cases 
    nmissy <- sum(!complete.cases(y)) # Count of incomplete cases
    
    missy <- 0
    if (nmissy>= 1){
      missy <- which(!complete.cases(y))}

    
    
    #REMOVE ROWS AND COLUMNS FOR YEAR WITH MISSING RESPONSE
    nyr_red <- nyr-nmissy
    y_red <- y[-missy]
    monmat_red <- monmat[-missy,]
    
    year_red <- year[-missy]
    lcm_red <- lcm[-missy,-missy]
    
    
    
    #nyr_red <- nyr-1
    #y_red <- y[-36]
    #monmat_red <- monmat[-36,]
  #  
  #  year_red <- year[-36]
  #  lcm_red <- lcm[-36,-36]
    
    #
    #BASE FUNCTION, LINEAR REGRESSION ON YEAR
    #
    #
    
    #CREATE DATA STRUCTURE NEED FOR RANDOM EFFECT SMOOTHING
    g_red <- rep(1,nyr_red)   # grouping variable
    datagp2a <- groupedData(y_red~year_red+monmat_red|g_red) #data structure for use in lme
    
    #FIT BASE MODEL WITHOUT WEATHER COVARIATES
    #lmixres0 <- lme(y_red~year_red,random=pdIdent(~-1),data=datagp2a,na.action=na.exclude,
    #                correlation=corCAR1(0.83, form=~year_red),method="ML",lmeControl(msVerbose=TRUE))
    
    fs0.nlme.fn <- function(input,beta0,beta1){
      value <- beta0 + beta1*input 
      value
    }
    
    tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red)
    tempfit <- gls(y_red~year_red,correlation=corAR1(form=~year_red | g_red),data=tempnogrpdata,method="ML")
    nlmecar1.start <- c(as.numeric(coef(tempfit)))
    names(nlmecar1.start) <- c("beta0", "beta1")
    
    final0model <- try(nlme(y_red~fs0.nlme.fn(year_red,beta0,beta1),
                            fixed=beta0+beta1~1,
                            random=beta0 ~ 1|g_red,
                            start=nlmecar1.start,correlation=corAR1(form=~year_red | g_red),
                            control=nlmeControl(opt="nlminb",returnObject=TRUE,maxIter=100,pnlsMaxIter=100,msMaxIter=100),method="ML"))
    if(class(final0model)[1]=="try-error"){ final0 <- 0}
    
    
    summary(final0model)
    logLik0 <- as.numeric(summary(final0model)["logLik"][1])
    neg_logLik0 <- - logLik0
    
    corr0 <- as.numeric(final0model$modelStruct$corStruct)
    corr0 <- 2*((1/(1+exp(-corr0)))-0.5)
    sigma_0 <- final0model$sigma
    sigma_sq_0 <- sigma_0^2
    VarCorr(final0model)
    
    #
    #LINEAR REGRESSION ON YEAR WITH AR1 ERRORS, 1ST ORDER FOURIER TERM AND EXPONENTIAL DECAY
    #
    
    fs2.nlme.fn <- function(input,scalec,beta0,beta1,beta2,beta3,beta4,monmat_red){
      decay <- 1/(1+exp(-scalec))
      Xmat <- array(0,dim=c(nyr_red,3))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
      }
      value <- beta0 + beta1*input + beta2*Xmat[,1] + beta3*Xmat[,2] + beta4*Xmat[,3]
      value
    }
    
    # This part is tricky, so firstly explore range of values for the decay parameter (or scaled, to be precise)
    neglogliks <- NULL
    tempx <- NULL
    scalec <- gridmin
    while(scalec < gridmax){
      decay <- 1/(1+exp(-scalec))
      Xmat <- array(0,dim=c(nyr_red,3))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
      }
      tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
      tempgrpdata <- groupedData(y_red~year_red+X.1+X.2+X.3|g_red,data=tempnogrpdata)
      tempfit <- try(lme(y_red~year_red+X.1+X.2+X.3,random=~1|g_red,correlation=corAR1(form=~year_red | g_red),data=tempgrpdata,method="ML"))
      scalec <- scalec+0.25
      if(class(tempfit)=="try-error"){
        next
      }
      neglogliks <- c(neglogliks,-logLik(tempfit))
      tempx <- c(tempx,scalec-0.25)   
    }
    if (drawgraphs==1){plot(tempx,neglogliks)}
    whichscalec <- tempx[which.min(neglogliks)]
    if (which.min(neglogliks)==length(neglogliks)) {
      if(try_no_decay_model){whichscalec <- no_decay_threshold+1}}    
    
    
    if(whichscalec > no_decay_threshold){ # Run model without decay
      Xmat <- array(0,dim=c(nyr_red,3))
      for(j in 1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] + monmat_red[ ,j]
        Xmat[ ,2] <- Xmat[ ,2] + monmat_red[ ,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] + monmat_red[ ,j]*cos(twopi*1*j/12)
      }
      g_red <- rep(1,nyr_red)   # grouping variable
      datanongp2a <- data.frame(X=Xmat,y_red=y_red,year_red=year_red,g_red=g_red)
      datagp2a <- groupedData(formula=y_red~year_red+X.1+X.2+X.3 |g_red,data=datanongp2a)
      finalCmodel <- lme(y_red~year_red+X.1+X.2+X.3,random=pdIdent(~-1),data=datagp2a,
                         na.action=na.exclude,correlation=corAR1(form=~year_red | g_red),method="ML")
      finalC <- -1
    }else{ # Run decay model
      # Get some starting values
      decay <- 1/(1+exp(-whichscalec))
      Xmat <- array(0,dim=c(nyr_red,3))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
      }
      tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
      tempfit <- gls(y_red~year_red+X.1+X.2+X.3,correlation=corAR1(form=~year_red | g_red),data=tempnogrpdata,method="ML")
      gnls.start <- c(whichscalec,as.numeric(coef(tempfit)))
      names(gnls.start) <- c("scalec", "beta0", "beta1", "beta2", "beta3", "beta4")
      gnlsresC <- try(gnls(y_red~fs2.nlme.fn(year_red,scalec,beta0,beta1,beta2,beta3,beta4,monmat_red),
                           start=gnls.start,data=datagp2a,correlation=corAR1(form=~year_red | g_red),
                           control=nls.control(maxiter=5000)))
      if(class(gnlsresC)[1]!="try-error"){
        nlmecar1.start <- as.numeric(coef(gnlsresC))
        corrCtemp <- as.numeric(gnlsresC$modelStruct$corStruct)
        corrCstart <- 2*((1/ (1+ exp(-corrCtemp)))-0.5)
      }else{
        nlmecar1.start <- gnls.start
        corrCtemp <- as.numeric(tempfit$modelStruct$corStruct)
        corrCstart <- 2*((1/ (1+ exp(-corrCtemp)))-0.5)
      }
      finalCmodel <- try(nlme(y_red~fs2.nlme.fn(year_red,scalec,beta0,beta1,beta2,beta3,beta4,monmat_red),
                              fixed=scalec+beta0+beta1+beta2+beta3+beta4~1,
                              random=beta0 ~ 1|g_red,
                              start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrCstart, form=~year_red | g_red),
                              control=nlmeControl(opt="nlminb"),method="ML"))
      if(class(finalCmodel)[1]=="try-error"){
        finalCmodel <- nlme(y_red~fs2.nlme.fn(year_red,scalec,beta0,beta1,beta2,beta3,beta4,monmat_red),
                            fixed=scalec+beta0+beta1+beta2+beta3+beta4~1,
                            random=beta0 ~ 1|g_red,
                            start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrCstart, form=~year_red | g_red),
                            control=nlmeControl(opt="nlminb",returnObject=TRUE,msVerbose=TRUE,maxIter=100,pnlsMaxIter=0,msMaxIter=100),method="ML")
        finalC <- 0
      }
    }
    
    # Extract parameter estimates and standard errors
    parC1.out <- coef(finalCmodel)
    if(whichscalec > no_decay_threshold){
      parC1 <- c(as.numeric(parC1.out),999)
      names(parC1) <- c(names(coef(finalCmodel)),"No decay")
      print(parC1)
      seC1 <- c(sqrt(diag(vcov(finalCmodel))),0)
      names(seC1) <- c(names(coef(finalCmodel)),"No decay")
      print(seC1)
    }else{
      print(parC1 <- parC1.out[c(2:6,1)])
      print(seC1 <- sqrt(diag(vcov(finalCmodel)))[c(2:6,1)])
    }
    corrC1 <- as.numeric(finalCmodel$modelStruct$corStruct)
    corrC1 <- 2*((1/(1+exp(-corrC1)))-0.5)
    (neg_log_lik_C1 <- -logLik(finalCmodel))
    (sigma_C1 <- finalCmodel$sigma)
    (sigma_sq_C1 <- sigma_C1^2)
    (ann_decay_C1 <- as.numeric((1/(1+exp(-parC1[6])))^12))
    #
    #LINEAR REGRESSION ON YEAR WITH AR1 ERRORS, 2ND ORDER FOURIER TERM AND EXPONENTIAL DECAY
    #
    
    #FUNCTION FOR 4-TERM FOURIER SERIES
    fs4.nlme.fn <- function(scaled,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,input){
      decay <- 1/(1+exp(-scaled))
      Xmat <- array(0,dim=c(nyr_red,5))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
        Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*2*j/12)
        Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*2*j/12)
      }
      value <- beta0 + beta1*input + beta2*Xmat[,1] + beta3*Xmat[,2] + beta4*Xmat[,3] + beta5*Xmat[,4] + beta6*Xmat[,5]
      value
    }
    
    # This part is tricky, so firstly explore range of values for the decay parameter (or scaled, to be precise)
    neglogliks <- NULL
    tempx <- NULL
    scaled <- gridmin
    while(scaled < gridmax){
      decay <- 1/(1+exp(-scaled))
      Xmat <- array(0,dim=c(nyr_red,5))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
        Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*2*j/12)
        Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*2*j/12)
      }
      tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
      tempgrpdata <- groupedData(y_red~year_red+X.1+X.2+X.3+X.4+X.5|g_red,data=tempnogrpdata)
      tempfit <- try(lme(y_red~year_red+X.1+X.2+X.3+X.4+X.5,random=~1|g_red,correlation=corAR1(form=~year_red | g_red),data=tempgrpdata,method="ML"))
      scaled <- scaled+0.25
      if(class(tempfit)=="try-error"){
        next
      }
      neglogliks <- c(neglogliks,-logLik(tempfit))
      tempx <- c(tempx,scaled-0.25)
    }
    if (drawgraphs==1){plot(tempx,neglogliks)}
    
    whichscaled <- tempx[which.min(neglogliks)]
    if (which.min(neglogliks)==length(neglogliks)) {
      if(try_no_decay_model){whichscaled <- no_decay_threshold+1}}    
    
    if(whichscaled > no_decay_threshold){ # Run model without decay
      Xmat <- array(0,dim=c(nyr_red,5))
      for(j in 1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] + monmat_red[ ,j]
        Xmat[ ,2] <- Xmat[ ,2] + monmat_red[ ,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] + monmat_red[ ,j]*cos(twopi*1*j/12)
        Xmat[ ,4] <- Xmat[ ,4] + monmat_red[ ,j]*sin(twopi*2*j/12)
        Xmat[ ,5] <- Xmat[ ,5] + monmat_red[ ,j]*cos(twopi*2*j/12)
      }
      g_red <- rep(1,nyr_red)   # grouping variable
      datanongp2a <- data.frame(X=Xmat,y_red=y_red,year_red=year_red,g_red=g_red)
      datagp2a <- groupedData(formula=y_red~year_red+X.1+X.2+X.3+X.4+X.5 |g_red,data=datanongp2a)
      finalDmodel <- lme(y_red~year_red+X.1+X.2+X.3+X.4+X.5,random=pdIdent(~-1),data=datagp2a,
                         na.action=na.exclude,correlation=corAR1(form=~year_red | g_red),method="ML")
      finalD <- -1
    }else{ # Run decay model
      # Get some starting values
      decay <- 1/(1+exp(-whichscaled))
      Xmat <- array(0,dim=c(nyr_red,5))
      for (j in  1:colmonmat) {
        Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
        Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
        Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
        Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*2*j/12)
        Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*2*j/12)
      }
      tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
      tempfit <- gls(y_red~year_red+X.1+X.2+X.3+X.4+X.5,correlation=corAR1(form=~year_red | g_red),data=tempnogrpdata,method="ML")
      gnls.start <- c(whichscaled,as.numeric(coef(tempfit))) #DAE corrected
      names(gnls.start) <- c("scaled", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6")
      gnlsresD <- try(gnls(y_red~fs4.nlme.fn(scaled,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                           start=gnls.start,data=datagp2a,correlation=corAR1(form=~year_red | g_red),
                           control=nls.control(maxiter=5000)))
      if(class(gnlsresD)[1]!="try-error"){
        nlmecar1.start <- as.numeric(coef(gnlsresD))
        corrDtemp <- as.numeric(gnlsresD$modelStruct$corStruct)        
        corrDstart <- 2*((1/ (1+ exp(-corrDtemp)))-0.5)
      }else{
        nlmecar1.start <- gnls.start
        corrDtemp <- as.numeric(tempfit$modelStruct$corStruct)
        corrDstart <- 2*((1/ (1+ exp(-corrDtemp)))-0.5)
      }
      finalDmodel <- try(nlme(y_red~fs4.nlme.fn(scaled,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                              fixed=scaled+beta0+beta1+beta2+beta3+beta4+beta5+beta6~1,
                              random=beta0 ~ 1|g_red,
                              start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrDstart, form=~year_red | g_red),
                              control=nlmeControl(opt="nlminb"),method="ML"))
      if(class(finalDmodel)[1]=="try-error"){
        finalDmodel <- nlme(y_red~fs4.nlme.fn(scaled,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                            fixed=scaled+beta0+beta1+beta2+beta3+beta4+beta5+beta6~1,
                            random=beta0 ~ 1|g_red,
                            start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrDstart, form=~year_red | g_red),
                            control=nlmeControl(opt="nlminb",returnObject=TRUE,msVerbose=TRUE,maxIter=100,pnlsMaxIter=0,msMaxIter=100),method="ML")
        finalD <- 0
      }
    }
    
    
    # Extract parameter estimates and standard errors
    parD1.out <- coef(finalDmodel)
    if(whichscaled > no_decay_threshold){
      parD1 <- c(as.numeric(parD1.out),999)
      names(parD1) <- c(names(coef(finalDmodel)),"No decay")
      print(parD1)
      seD1 <- c(sqrt(diag(vcov(finalDmodel))),0)
      names(seD1) <- c(names(coef(finalDmodel)),"No decay")
      print(seD1)
    }else{
      print(parD1 <- parD1.out[c(2:8,1)])
      print(seD1 <- sqrt(diag(vcov(finalDmodel)))[c(2:8,1)])
    }
    (neg_log_lik_D1 <- -logLik(finalDmodel))
    (sigma_D1 <- finalDmodel$sigma)
    (sigma_sq_D1 <- sigma_D1^2)
    corrD1 <- as.numeric(finalDmodel$modelStruct$corStruct)
    corrD1 <- 2*((1/(1+exp(-corrD1)))-0.5)
    (ann_decay_D1 <- as.numeric((1/(1+exp(-parD1[8])))^12))
    
    #
    #LINEAR REGRESSION ON YEAR WITH AR1 ERRORS, 2ND ORDER FOURIER TERM AND EXPONENTIAL DECAY
    # 2nd order term with period 0.5 years
    # E block of code
    #
    
    if (fitE==1){
      #FUNCTION FOR 4-TERM FOURIER SERIES
      fs5.nlme.fn <- function(scalee,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,input){
        decay <- 1/(1+exp(-scalee))
        Xmat <- array(0,dim=c(nyr_red,5))
        for (j in  1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
          Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*0.5*j/12)
          Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*0.5*j/12)
        }
        value <- beta0 + beta1*input + beta2*Xmat[,1] + beta3*Xmat[,2] + beta4*Xmat[,3] + beta5*Xmat[,4] + beta6*Xmat[,5]
        value
      }
      
      # This part is tricky, so firstly explore range of values for the decay parameter (or scaled, to be precise)
      neglogliks <- NULL
      tempx <- NULL
      scalee <- gridmin
      while(scalee < gridmax){
        decay <- 1/(1+exp(-scalee))
        Xmat <- array(0,dim=c(nyr_red,5))
        for (j in  1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
          Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*0.5*j/12)
          Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*0.5*j/12)
        }
        tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
        tempgrpdata <- groupedData(y_red~year_red+X.1+X.2+X.3+X.4+X.5|g_red,data=tempnogrpdata)
        tempfit <- try(lme(y_red~year_red+X.1+X.2+X.3+X.4+X.5,random=~1|g_red,correlation=corAR1(form=~year_red | g_red),data=tempgrpdata,method="ML"))
        scalee <- scalee+0.25
        if(class(tempfit)=="try-error"){
          next
        }
        neglogliks <- c(neglogliks,-logLik(tempfit))
        tempx <- c(tempx,scalee-0.25)
      }
      if (drawgraphs==1){plot(tempx,neglogliks)}
      whichscalee <- tempx[which.min(neglogliks)]
      if (which.min(neglogliks)==length(neglogliks)) {
        if(try_no_decay_model){whichscalee <- no_decay_threshold+1}}   
      
      
      if(whichscalee > no_decay_threshold){ # Run model without decay
        Xmat <- array(0,dim=c(nyr_red,5))
        for(j in 1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] + monmat_red[ ,j]
          Xmat[ ,2] <- Xmat[ ,2] + monmat_red[ ,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] + monmat_red[ ,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] + monmat_red[ ,j]*sin(twopi*0.5*j/12)
          Xmat[ ,5] <- Xmat[ ,5] + monmat_red[ ,j]*cos(twopi*0.5*j/12)
        }
        g_red <- rep(1,nyr_red)   # grouping variable
        datanongp2a <- data.frame(X=Xmat,y_red=y_red,year_red=year_red,g_red=g_red)
        datagp2a <- groupedData(formula=y_red~year_red+X.1+X.2+X.3+X.4+X.5 |g_red,data=datanongp2a)
        finalEmodel <- lme(y_red~year_red+X.1+X.2+X.3+X.4+X.5,random=pdIdent(~-1),data=datagp2a,
                           na.action=na.exclude,correlation=corAR1(form=~year_red | g_red),method="ML")
        finalE <- -1
      }else{ # Run decay model
        # Get some starting values
        decay <- 1/(1+exp(-whichscalee))
        Xmat <- array(0,dim=c(nyr_red,5))
        for (j in  1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] +  (decay^(j-1))*monmat_red[,j]
          Xmat[ ,2] <- Xmat[ ,2] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] +  (decay^(j-1))*monmat_red[,j]*sin(twopi*0.5*j/12)
          Xmat[ ,5] <- Xmat[ ,5] +  (decay^(j-1))*monmat_red[,j]*cos(twopi*0.5*j/12)
        }
        tempnogrpdata <- data.frame(y_red=y_red,g_red=g_red,year_red=year_red,X=Xmat)
        tempfit <- gls(y_red~year_red+X.1+X.2+X.3+X.4+X.5,correlation=corAR1(form=~year_red | g_red),data=tempnogrpdata,method="ML")
        gnls.start <- c(whichscalee,as.numeric(coef(tempfit)))
        names(gnls.start) <- c("scalee", "beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6")
        gnlsresE <- try(gnls(y_red~fs5.nlme.fn(scalee,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                             start=gnls.start,data=datagp2a,correlation=corAR1(form=~year_red | g_red),
                             control=nls.control(maxiter=5000)))
        if(class(gnlsresE)[1]!="try-error"){
          nlmecar1.start <- as.numeric(coef(gnlsresE))
          corrEtemp <- as.numeric(gnlsresE$modelStruct$corStruct)
          corrEstart <- 2*((1/ (1+ exp(-corrEtemp)))-0.5)
        }else{
          nlmecar1.start <- gnls.start
          corrEtemp <- as.numeric(tempfit$modelStruct$corStruct)
          corrEstart <- 2*((1/ (1+ exp(-corrEtemp)))-0.5)
        }
        finalEmodel <- try(nlme(y_red~fs5.nlme.fn(scalee,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                                fixed=scalee+beta0+beta1+beta2+beta3+beta4+beta5+beta6~1,
                                random=beta0 ~ 1|g_red,
                                start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrEstart, form=~year_red | g_red),
                                control=nlmeControl(opt="nlminb"),method="ML"))
        if(class(finalEmodel)[1]=="try-error"){
          finalEmodel <- nlme(y_red~fs5.nlme.fn(scalee,beta0,beta1,beta2,beta3,beta4,beta5,beta6,monmat_red,year_red),
                              fixed=scalee+beta0+beta1+beta2+beta3+beta4+beta5+beta6~1,
                              random=beta0 ~ 1|g_red,
                              start=nlmecar1.start,data=datagp2a,correlation=corAR1(corrEstart, form=~year_red | g_red),
                              control=nlmeControl(opt="nlminb",returnObject=TRUE,msVerbose=TRUE,maxIter=100,pnlsMaxIter=0,msMaxIter=100),method="ML")
          finalE <- 0
        }
      }
      
      
      # Extract parameter estimates and standard errors
      parE1.out <- coef(finalEmodel)
      if(whichscalee > no_decay_threshold){
        parE1 <- c(as.numeric(parD1.out),999)
        names(parE1) <- c(names(coef(finalEmodel)),"No decay")
        print(parE1)
        seE1 <- c(sqrt(diag(vcov(finalEmodel))),0)
        names(seE1) <- c(names(coef(finalEmodel)),"No decay")
        print(seE1)
      }else{
        print(parE1 <- parE1.out[c(2:8,1)])
        print(seE1 <- sqrt(diag(vcov(finalEmodel)))[c(2:8,1)])
      }
      (neg_log_lik_E1 <- -logLik(finalEmodel))
      (sigma_E1 <- finalEmodel$sigma)
      (sigma_sq_E1 <- sigma_E1^2)
      corrE1 <- as.numeric(finalEmodel$modelStruct$corStruct)
      corrE1 <- 2*((1/(1+exp(-corrE1)))-0.5)
      (ann_decay_E1 <- as.numeric((1/(1+exp(-parE1[8])))^12))
      
      
      #
      # All Quantile calculations and graphs in here
      #
      if (calcQuant==1){
        
        
        # The Fourier parts (non-decaying and then decaying)
        # Note: deriving the CIs manually because the sin/cos breaks them in estimable
        if(whichscalec > no_decay_threshold){
          estf <- estimable(finalCmodel,cm=cbind(0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12)))
        }else{
          estf <- estimable(finalCmodel,cm=cbind(0,0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12)))
        }
        esft.pred <- estf$Estimate
        estf.df <- max(estf$DF)
        estf.se <- estf$"Std. Error"
        estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
        estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
        estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
        estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
        (quant_Fnd_C <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        
        # Now decaying
        quant_F_C <- quant_Fnd_C
        if(!(whichscalec > no_decay_threshold)){
          decayc <- 1/(1+exp(-coef(finalCmodel)$scalec))
          estf <- estimable(finalCmodel,cm=cbind(0,0,0,(decayc^(x60-1)),(decayc^(x60-1))*sin(twopi*(x60)/12),(decayc^(x60-1))*cos(twopi*(x60)/12)))
          esft.pred <- estf$Estimate
          estf.df <- max(estf$DF)
          estf.se <- estf$"Std. Error"
          estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
          estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
          estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
          estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
          (quant_F_C <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        }
        
        # The fitted values
        if(whichscalec > no_decay_threshold){
          decayc <- 1
        }
        Xmat <- array(0,dim=c(nyr,3))
        for(j in 1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] + (decayc^(j-1))*monmat[ ,j]
          Xmat[ ,2] <- Xmat[ ,2] + (decayc^(j-1))*monmat[ ,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] + (decayc^(j-1))*monmat[ ,j]*cos(twopi*1*j/12)
        }
        if(whichscalec > no_decay_threshold){
          fv.est <- estimable(finalCmodel,cm=cbind("(Intercept)"=1,year_red=year,X.1=Xmat[,1],X.2=Xmat[,2],X.3=Xmat[,3]))
        }else{
          fv.est <- estimable(finalCmodel,cm=cbind(scalec=0,beta0=1,beta1=year,beta2=Xmat[,1],beta3=Xmat[,2],beta4=Xmat[,3]))
        }
        fv.pred <- fv.est$Estimate
        fv.df <- max(fv.est$DF)
        fv.se <- fv.est$"Std. Error"
        fv50.lci <- fv.est$Estimate+qt(0.25,fv.df)*fv.se
        fv50.uci <- fv.est$Estimate+qt(0.75,fv.df)*fv.se
        fv95.lci <- fv.est$Estimate+qt(0.025,fv.df)*fv.se
        fv95.uci <- fv.est$Estimate+qt(0.975,fv.df)*fv.se
        (quant_fv_C <- cbind(fv95.lci,fv50.lci,fv.pred,fv50.uci,fv95.uci))
        
        
        
        # The Fourier parts (non-decaying and then decaying)
        # Note: deriving the CIs manually because the sin/cos breaks them in estimable
        if(whichscaled > no_decay_threshold){
          cm.est <- cbind(0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12),sin(2*twopi*(x60)/12),cos(2*twopi*(x60)/12))
        }else{
          cm.est <- cbind(0,0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12),sin(2*twopi*(x60)/12),cos(2*twopi*(x60)/12))
        }
        estf <- estimable(finalDmodel,cm=cm.est)
        esft.pred <- estf$Estimate
        estf.df <- max(estf$DF)
        estf.se <- estf$"Std. Error"
        estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
        estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
        estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
        estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
        (quant_Fnd_D <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        
        quant_F_D <- quant_Fnd_D
        if(!(whichscaled > no_decay_threshold)){
          # Now decaying
          decayd <- 1/(1+exp(-coef(finalDmodel)$scaled))
          cm.est <- cbind(0,0,0,(decayd^(x60-1)),(decayd^(x60-1))*sin(twopi*(x60)/12),(decayd^(x60-1))*cos(twopi*(x60)/12),(decayd^(x60-1))*sin(2*twopi*(x60)/12),(decayd^(x60-1))*cos(2*twopi*(x60)/12))
          estf <- estimable(finalDmodel,cm=cm.est)
          esft.pred <- estf$Estimate
          estf.df <- max(estf$DF)
          estf.se <- estf$"Std. Error"
          estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
          estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
          estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
          estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
          (quant_F_D <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        }
        
        # The fitted values
        if(whichscaled > no_decay_threshold){
          decayd <- 1
        }
        Xmat <- array(0,dim=c(nyr,5))
        for(j in 1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] + (decayd^(j-1))*monmat[ ,j]
          Xmat[ ,2] <- Xmat[ ,2] + (decayd^(j-1))*monmat[ ,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] + (decayd^(j-1))*monmat[ ,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] + (decayd^(j-1))*monmat[ ,j]*sin(twopi*2*j/12)
          Xmat[ ,5] <- Xmat[ ,5] + (decayd^(j-1))*monmat[ ,j]*cos(twopi*2*j/12)
        }
        if(whichscaled > no_decay_threshold){
          cm.est <- cbind("(Intercept)"=1,year_red=year,X.1=Xmat[,1],X.2=Xmat[,2],X.3=Xmat[,3],X.4=Xmat[,4],X.5=Xmat[,5])
        }else{
          cm.est <- cbind(scaled=0,beta0=1,beta1=year,beta2=Xmat[,1],beta3=Xmat[,2],beta4=Xmat[,3],beta5=Xmat[,4],beta6=Xmat[,5])
        }
        fv.est <- estimable(finalDmodel,cm=cm.est)
        fv.pred <- fv.est$Estimate
        fv.df <- max(fv.est$DF)
        fv.se <- fv.est$"Std. Error"
        fv50.lci <- fv.est$Estimate+qt(0.25,fv.df)*fv.se
        fv50.uci <- fv.est$Estimate+qt(0.75,fv.df)*fv.se
        fv95.lci <- fv.est$Estimate+qt(0.025,fv.df)*fv.se
        fv95.uci <- fv.est$Estimate+qt(0.975,fv.df)*fv.se
        (quant_fv_D <- cbind(fv95.lci,fv50.lci,fv.pred,fv50.uci,fv95.uci))
        
        
        
        
        # The Fourier parts (non-decaying and then decaying)
        # Note: deriving the CIs manually because the sin/cos breaks them in estimable
        if(whichscalee > no_decay_threshold){
          cm.est <- cbind(0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12),sin(0.5*twopi*(x60)/12),cos(0.5*twopi*(x60)/12))
        }else{
          cm.est <- cbind(0,0,0,1,sin(twopi*(x60)/12),cos(twopi*(x60)/12),sin(0.5*twopi*(x60)/12),cos(0.5*twopi*(x60)/12))
        }
        estf <- estimable(finalEmodel,cm=cm.est)
        esft.pred <- estf$Estimate
        estf.df <- max(estf$DF)
        estf.se <- estf$"Std. Error"
        estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
        estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
        estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
        estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
        (quant_Fnd_E <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        
        quant_F_E <- quant_Fnd_E
        if(!(whichscalee > no_decay_threshold)){
          # Now decaying
          decaye <- 1/(1+exp(-coef(finalEmodel)$scalee))
          cm.est <- cbind(0,0,0,(decaye^(x60-1)),(decaye^(x60-1))*sin(twopi*(x60)/12),(decaye^(x60-1))*cos(twopi*(x60)/12),(decaye^(x60-1))*sin(0.5*twopi*(x60)/12),(decaye^(x60-1))*cos(0.5*twopi*(x60)/12))
          estf <- estimable(finalEmodel,cm=cm.est)
          esft.pred <- estf$Estimate
          estf.df <- max(estf$DF)
          estf.se <- estf$"Std. Error"
          estf50.lci <- estf$Estimate+qt(0.25,estf.df)*estf.se
          estf50.uci <- estf$Estimate+qt(0.75,estf.df)*estf.se
          estf95.lci <- estf$Estimate+qt(0.025,estf.df)*estf.se
          estf95.uci <- estf$Estimate+qt(0.975,estf.df)*estf.se
          (quant_F_E <- cbind(estf95.lci,estf50.lci,esft.pred,estf50.uci,estf95.uci))
        }
        
        # The fitted values
        if(whichscalee > no_decay_threshold){
          decaye <- 1
        }
        Xmat <- array(0,dim=c(nyr,5))
        for(j in 1:colmonmat) {
          Xmat[ ,1] <- Xmat[ ,1] + (decaye^(j-1))*monmat[ ,j]
          Xmat[ ,2] <- Xmat[ ,2] + (decaye^(j-1))*monmat[ ,j]*sin(twopi*1*j/12)
          Xmat[ ,3] <- Xmat[ ,3] + (decaye^(j-1))*monmat[ ,j]*cos(twopi*1*j/12)
          Xmat[ ,4] <- Xmat[ ,4] + (decaye^(j-1))*monmat[ ,j]*sin(twopi*0.5*j/12)
          Xmat[ ,5] <- Xmat[ ,5] + (decaye^(j-1))*monmat[ ,j]*cos(twopi*0.5*j/12)
        }
        if(whichscalee > no_decay_threshold){
          cm.est <- cbind("(Intercept)"=1,year_red=year,X.1=Xmat[,1],X.2=Xmat[,2],X.3=Xmat[,3],X.4=Xmat[,4],X.5=Xmat[,5])
        }else{
          cm.est <- cbind(scalee=0,beta0=1,beta1=year,beta2=Xmat[,1],beta3=Xmat[,2],beta4=Xmat[,3],beta5=Xmat[,4],beta6=Xmat[,5])
        }
        fv.est <- estimable(finalEmodel,cm=cm.est)
        fv.pred <- fv.est$Estimate
        fv.df <- max(fv.est$DF)
        fv.se <- fv.est$"Std. Error"
        fv50.lci <- fv.est$Estimate+qt(0.25,fv.df)*fv.se
        fv50.uci <- fv.est$Estimate+qt(0.75,fv.df)*fv.se
        fv95.lci <- fv.est$Estimate+qt(0.025,fv.df)*fv.se
        fv95.uci <- fv.est$Estimate+qt(0.975,fv.df)*fv.se
        (quant_fv_E <- cbind(fv95.lci,fv50.lci,fv.pred,fv50.uci,fv95.uci))
        
        #close fitE
      }
      
      
      
      ######################################################################################################
      # OUTPUT PLOTS - start with single sin/cos
      
      # cut them all out
      if (drawgraphs==1){
        
        if(!(whichscalec > no_decay_threshold)){
          # Decaying Fourier series
          ymax <- max(quant_F_C[,5])
          ymin <- min(quant_F_C[,1])
          plot(y=quant_F_C[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
          lines(y=quant_F_C[,1],x=x60,col="grey")
          lines(y=quant_F_C[,5],x=x60,col="grey")
        }
        
        # Fourier series without decay
        ymax <- max(quant_Fnd_C[,5])
        ymin <- min(quant_Fnd_C[,1])
        plot(y=quant_Fnd_C[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
        lines(y=quant_Fnd_C[,1],x=x60,col="grey")
        lines(y=quant_Fnd_C[,5],x=x60,col="grey")
        
        # Fourier series without decay for one year only
        ymax <- max(quant_Fnd_C[,5])
        ymin <- min(quant_Fnd_C[,1])
        plot(y=quant_Fnd_C[1:12,"esft.pred"],x=x60[1:12],ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
        lines(y=quant_Fnd_C[1:12,1],x=x60[1:12],col="grey")
        lines(y=quant_Fnd_C[1:12,5],x=x60[1:12],col="grey")
        
        #PLOT DATA AND FITTED VALUES
        ymax <- max(quant_fv_C[,5],y_red)
        ymin <- min(quant_fv_C[,1],y_red)
        plot(y=y,x=year,xlab="Year",ylab="Index",col="red",main="C",ylim=c(ymin,ymax))
        lines(y=quant_fv_C[,"fv.pred"],x=year,col="black")
        lines(y=quant_fv_C[,1],x=year,col="grey")
        lines(y=quant_fv_C[,5],x=year,col="grey")
        
        
        # OUTPUT PLOTS - now the second order Fourier series
        
        if(!(whichscaled > no_decay_threshold)){
          # Decaying Fourier series
          ymax <- max(quant_F_D[,5])
          ymin <- min(quant_F_D[,1])
          plot(y=quant_F_D[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
          lines(y=quant_F_D[,1],x=x60,col="grey")
          lines(y=quant_F_D[,5],x=x60,col="grey")
        }
        
        # Fourier series without decay
        ymax <- max(quant_Fnd_D[,5])
        ymin <- min(quant_Fnd_D[,1])
        plot(y=quant_Fnd_D[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
        lines(y=quant_Fnd_D[,1],x=x60,col="grey")
        lines(y=quant_Fnd_D[,5],x=x60,col="grey")
        
        #12 Months plot without decay (pure Fourier form)
        ymax <- max(quant_Fnd_D[,5])
        ymin <- min(quant_Fnd_D[,1])
        plot(y=quant_Fnd_D[1:12,"esft.pred"],x=x60[1:12],ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
        lines(y=quant_Fnd_D[1:12,1],x=x60[1:12],col="grey")
        lines(y=quant_Fnd_D[1:12,5],x=x60[1:12],col="grey")
        
        #PLOT DATA AND FITTED VALUES
        ymax <- max(quant_fv_D[,5],y_red)
        ymin <- min(quant_fv_D[,1],y_red)
        plot(y=y,x=year,xlab="Year",ylab="Index",col="red",main="D",ylim=c(ymin,ymax))
        lines(y=quant_fv_D[,"fv.pred"],x=year,col="black")
        lines(y=quant_fv_D[,1],x=year,col="grey")
        lines(y=quant_fv_D[,5],x=year,col="grey")
        
        # OUTPUT PLOTS - now the second order Fourier series
        # 2nd term of period 0.5 years
        
        if (fitE==1) {
          if(!(whichscalee > no_decay_threshold)){
            # Decaying Fourier series
            ymax <- max(quant_F_E[,5])
            ymin <- min(quant_F_E[,1])
            plot(y=quant_F_E[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
            lines(y=quant_F_E[,1],x=x60,col="grey")
            lines(y=quant_F_E[,5],x=x60,col="grey")
          }
          
          # Fourier series with decay
          ymax <- max(quant_Fnd_E[,5])
          ymin <- min(quant_Fnd_E[,1])
          plot(y=quant_Fnd_E[,"esft.pred"],x=x60,ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
          lines(y=quant_Fnd_E[,1],x=x60,col="grey")
          lines(y=quant_Fnd_E[,5],x=x60,col="grey")
          
          #12 Months plot without decay (pure Fourier form)
          ymax <- max(quant_Fnd_E[,5])
          ymin <- min(quant_Fnd_E[,1])
          plot(y=quant_Fnd_E[1:12,"esft.pred"],x=x60[1:12],ylab="Parameter",xlab="Month",type="l",ylim=c(ymin,ymax))
          lines(y=quant_Fnd_E[1:12,1],x=x60[1:12],col="grey")
          lines(y=quant_Fnd_E[1:12,5],x=x60[1:12],col="grey")
          
          #PLOT DATA AND FITTED VALUES
          ymax <- max(quant_fv_E[,5],y_red)
          ymin <- min(quant_fv_E[,1],y_red)
          plot(y=y,x=year,xlab="Year",ylab="Index",col="red",main="E",ylim=c(ymin,ymax))
          lines(y=quant_fv_E[,"fv.pred"],x=year,col="black")
          lines(y=quant_fv_E[,1],x=year,col="grey")
          lines(y=quant_fv_E[,5],x=year,col="grey")
          
          #close skipping E block
        }
        
        
        #close of graph-skipping
      }
      
      
      #close quantile block
    }
    
    #LOG LIKELIHOODS FOR ALL MODELS
    if (fitE==1){
      cbind(neg_logLik0,neg_log_lik_C1,neg_log_lik_D1,neg_log_lik_E1)}
    if (fitE==0){
      cbind(neg_logLik0,neg_log_lik_C1,neg_log_lik_D1)}
    
    #PARAMETER ESTIMATES, CORRELATIONS AND ANNUAL DECAY FOR MODELS C AND D
    #correlation not printed if optimisation not secure
    # Check not needed if no nonlinear part
    if(whichscalec > no_decay_threshold){
      print(cbind(parC1,seC1))
      parC1O <- parC1
    }else{
      print(cbind(t(parC1),seC1))
      parC1O <- t(parC1)
    }
    print(corrC1)
    #print(corrC1_check) # if missing, indication that fit not very reliable
    print(cbind(sigma_sq_C1,sigma_C1,ann_decay_C1,neg_log_lik_C1))
    
    if(whichscaled > no_decay_threshold){
      print(cbind(parD1,seD1))
      parD1O <- parD1
    }else{
      print(cbind(t(parD1),seD1))
      parD1O <- t(parD1)
    }
    print(corrD1)
    #print(corrD1_check) # if missing, indication that fit not very reliable
    print(cbind(sigma_sq_D1,sigma_D1,ann_decay_D1,neg_log_lik_D1))
    
    if (fitE==1){
      if(whichscalee > no_decay_threshold){
        print(cbind(parE1,seE1))
        parE1O <- parE1
      }else{
        print(cbind(t(parE1),seE1))
        parE1O <- t(parE1)
      }
      print(corrE1)
      #print(corrE1_check) # if missing, indication that fit not very reliable
      print(cbind(sigma_sq_E1,sigma_E1,ann_decay_E1,neg_log_lik_E1))
      
      #close fitE
    }
    
    # save results for future use
    if (fitE==1){
      save( 
        neg_logLik0,neg_log_lik_C1,neg_log_lik_D1,neg_log_lik_E1,
        final0,finalC,finalD,finalE,
        whichscalec,whichscaled,whichscalee,
        sigma_0,sigma_sq_0,corr0,
        sigma_C1,sigma_sq_C1,corrC1,ann_decay_C1,
        sigma_D1,sigma_sq_D1,corrD1,ann_decay_D1,
        sigma_E1,sigma_sq_E1,corrE1,ann_decay_E1,
        parC1O,seC1,parD1O,seD1,parE1O,seE1,
        quant_fv_C,quant_Fnd_C,quant_F_C,
        quant_fv_D,quant_Fnd_D,quant_F_D,
        quant_fv_E,quant_Fnd_E,quant_F_E,
        final0model,finalCmodel,finalDmodel,finalEmodel,
        file=savenameE)
    }
    
    # don't save these, source of error# quant_F_C,quant_F_D,
    if (fitE==0){
      save(  
        neg_logLik0,neg_log_lik_C1,neg_log_lik_D1,
        final0,finalC,finalD,
        whichscalec,whichscaled,
        sigma_0,sigma_sq_0,corr0,
        sigma_C1,sigma_sq_C1,corrC1,ann_decay_C1,
        sigma_D1,sigma_sq_D1,corrD1,ann_decay_D1,
        parC1O,seC1,parD1O,seD1,
        quant_fv_C,quant_Fnd_C,quant_F_C,
        quant_fv_D,quant_Fnd_D,quant_F_D,
        final0model,finalCmodel,finalDmodel,
        file=savenamenoE)
    }
    
    
    
    #
    # close species and covariate loop
  }}


#
#check can recover objects
#load("bullf_meanT4,E.Rsave")

#
#}
#--
#Biomathematics and Statistics Scotland (BioSS) is formally part of The
#James Hutton Institute (JHI), a registered Scottish charity No. SC041796
#and a company limited by guarantee No. SC374831
