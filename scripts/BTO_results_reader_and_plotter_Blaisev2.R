#R PROGRAM FOR CONSTRAINING MONTHLY REGRESSION COEFFICIENTS TO LIE ON DAMPED FOURIER CURVES,
#TREATING MODEL HAS HAVING 1 NON-LINEAR PARAMETER, THE REMAINDER BEING LINEAR
#CONDITIONAL ON THE VALUE OF THE EXPONENTIAL DECAY PARAMETER
#AS SENT TO BTO
#

# THIS VERSION RUNS USING NLME

#
#CALL NECESSARY LIBRARIES
#

rm(list=ls(all=TRUE))

#library(Matrix)
library(nlme)
library(MASS)
#library(numDeriv)
library(Hmisc)
library(gmodels)




#READ DATA
###Declare directories to read from
#


#1 USE BTO DATA,
#other options to be declared for other data sets
datadir <- "C:/Users/blaise/Desktop/BICCO-Net/Fourier bioss 04.04/MeanTemp Results"
resdir <- "C:/Users/blaise/Desktop/BICCO-Net/Fourier bioss 04.04/MeanTemp Results"


#resdir <- c("C:/dae/temp_13_jan/bicconet/national_index_work/bto/min_runs2",
#            "C:/dae/temp_13_jan/bicconet/national_index_work/bto/max_runs2",
#            "C:/dae/temp_13_jan/bicconet/national_index_work/bto/mean_runs2",
#            "C:/dae/temp_13_jan/bicconet/national_index_work/bto/rain_runs2")

setwd(datadir)
birds <- read.csv("BTOallspecies.csv",header=TRUE)
birdnames <- as.character(birds$speciesname)
birdcodes <- as.character(birds$Decode)
cbind(birdnames,birdcodes)

#
###Choose a covariate file
# x_file_pick <- 3 #choose which weather variable to use
#

try_no_decay_model <- TRUE
no_decay_threshold <- 7


#
#loop round met variable and species, remember to remove / close { } if active or not
#
# imet <- 1 #choose which weather variable to use
# isp <- 1 #choose which bird species to use
# 69 to choose from
#

# READ RESULTS

spno<-c(1:7,9:147,149:162,164:276,278:351,353:375,377:385,387:409,411:418,420,421,423:427,429:476,478:494)
lspno<-length(spno)
vno<-1
taxa<-factor(c(rep("birds",84),rep("moths",length(85:346)),rep("aphids",length(347:419)),rep("butterflies",length(420:470)),rep("mammals",length(471:lspno))))

metnames <- "meanT"

all_neg_log_lik_0 <- array(0,dim=c(lspno,vno))
all_neg_log_lik_C <- array(0,dim=c(lspno,vno))
all_neg_log_lik_D <- array(0,dim=c(lspno,vno))
all_neg_log_lik_E <- array(0,dim=c(lspno,vno))
all_sigma_0 <- array(0,dim=c(lspno,vno))
all_sigma_sq_0 <- array(0,dim=c(lspno,vno))
all_corr_0 <- array(0,dim=c(lspno,vno))
all_sigma_C1 <- array(0,dim=c(lspno,vno))
all_sigma_sq_C1 <- array(0,dim=c(lspno,vno))
all_corr_C1 <- array(0,dim=c(lspno,vno))
all_ann_decay_C1 <- array(0,dim=c(lspno,vno))
all_sigma_D1 <- array(0,dim=c(lspno,vno))
all_sigma_sq_D1 <- array(0,dim=c(lspno,vno))
all_corr_D1 <- array(0,dim=c(lspno,vno))
all_ann_decay_D1 <- array(0,dim=c(lspno,vno))
all_sigma_E1 <- array(0,dim=c(lspno,vno))
all_sigma_sq_E1 <- array(0,dim=c(lspno,vno))
all_corr_E1 <- array(0,dim=c(lspno,vno))
all_ann_decay_E1 <- array(0,dim=c(lspno,vno))
all_par_C1 <- array(0,dim=c(lspno,vno,6))
all_par_C1_PCA <- array(0,dim=c(lspno,4))
all_par_D1 <- array(0,dim=c(lspno,vno,8))
all_par_D1_PCA <- array(0,dim=c(lspno,6))
all_par_E1 <- array(0,dim=c(lspno,vno,8))
all_par_E1_PCA <- array(0,dim=c(lspno,6))
all_final_C <- array(0,dim=c(lspno,vno))
all_final_D <- array(0,dim=c(lspno,vno))
all_final_E <- array(0,dim=c(lspno,vno))
all_final_Fnd_C <- array(0,dim=c(lspno,12))

imet<-1
  setwd(resdir)
  for (isp in spno){
    metname <- metnames[imet]
    birdname <- birdnames[isp]
    
    savenameE <- paste(birdname,"_",metname,isp,",E.Rsave",sep="")
    savenamenoE <- paste(birdname,"_",metname,isp,",noE.Rsave",sep="")
    
    load(file=savenameE)
    
    all_neg_log_lik_C[c(1:lspno)[spno==isp],imet] <- neg_log_lik_C1
    all_neg_log_lik_D[c(1:lspno)[spno==isp],imet] <- neg_log_lik_D1
    all_neg_log_lik_E[c(1:lspno)[spno==isp],imet] <- neg_log_lik_E1
    all_neg_log_lik_0[c(1:lspno)[spno==isp],imet] <- neg_logLik0
    all_sigma_0[c(1:lspno)[spno==isp],imet] <-    sigma_0
    all_sigma_sq_0[c(1:lspno)[spno==isp],imet] <-   sigma_sq_0
    all_corr_0[c(1:lspno)[spno==isp],imet] <-  corr0
    all_sigma_C1[c(1:lspno)[spno==isp],imet] <-   sigma_C1
    all_sigma_sq_C1[c(1:lspno)[spno==isp],imet] <-   sigma_sq_C1
    all_corr_C1[c(1:lspno)[spno==isp],imet] <-  corrC1
    all_ann_decay_C1[c(1:lspno)[spno==isp],imet] <-  ann_decay_C1
    all_sigma_D1[c(1:lspno)[spno==isp],imet] <-  sigma_D1
    all_sigma_sq_D1[c(1:lspno)[spno==isp],imet] <-   sigma_sq_D1 
    all_corr_D1[c(1:lspno)[spno==isp],imet] <- corrD1
    all_ann_decay_D1[c(1:lspno)[spno==isp],imet] <-  ann_decay_D1
    all_sigma_E1[c(1:lspno)[spno==isp],imet] <-  sigma_E1
    all_sigma_sq_E1[c(1:lspno)[spno==isp],imet] <-  sigma_sq_E1
    all_corr_E1[c(1:lspno)[spno==isp],imet] <-   corrE1     
    all_ann_decay_E1[c(1:lspno)[spno==isp],imet] <-  ann_decay_E1
    all_par_C1[c(1:lspno)[spno==isp],imet,1:6] <- parC1O[1:6]
    all_par_D1[c(1:lspno)[spno==isp],imet,1:8] <- parD1O[1:8]
    all_par_E1[c(1:lspno)[spno==isp],imet,1:8] <- parE1O[1:8]
   all_par_C1_PCA[c(1:lspno)[spno==isp],((4*(imet-1)+1):(4*(imet-1)+4))] <- parC1O[3:6]
   all_par_D1_PCA[c(1:lspno)[spno==isp],((6*(imet-1)+1):(6*(imet-1)+6))] <- parD1O[3:8]
   all_par_E1_PCA[c(1:lspno)[spno==isp],((6*(imet-1)+1):(6*(imet-1)+6))] <- parE1O[3:8]
    all_final_C[c(1:lspno)[spno==isp],imet] <- finalC
    all_final_D[c(1:lspno)[spno==isp],imet] <- finalD
    all_final_E[c(1:lspno)[spno==isp],imet] <- finalE
    if (finalC + finalD + finalE == -3) {print(cbind(c(1:lspno)[spno==isp],imet,finalC,finalD,finalE))}
    all_final_Fnd_C[c(1:lspno)[spno==isp],1:12] <- quant_Fnd_C[1:12,3]
    #whichscalec,whichscaled,whichscalee,
  }#}

# likelihood ratio tests
all_lrt_0C <- array(0,dim=c(lspno,1))
all_lrt_CD <- array(0,dim=c(lspno,1))
all_lrt_CE <- array(0,dim=c(lspno,1))

all_lrt_0C <- 2*(all_neg_log_lik_0-all_neg_log_lik_C)
all_lrt_CD <- 2*(all_neg_log_lik_C-all_neg_log_lik_D)
all_lrt_CE <- 2*(all_neg_log_lik_C-all_neg_log_lik_E)

#
#################

#check form of some outputs
print(all_final_C)
print(all_final_D)
print(all_final_E)



#check form of some outputs
cbind(parC1O,seC1)
cbind(parD1O,seD1)
cbind(parE1O,seE1)




#check some values, compare with values from model runs
cbind(birdnames[spno],all_lrt_0C[,1:vno])
cbind(birdnames[spno],all_lrt_CD[,1:vno])
cbind(birdnames[spno],all_lrt_CE[,1:vno])

cbind(birdnames[spno],all_neg_log_lik_C[,1:vno])
picksp <- 17
pickmet <- 1
cbind(birdnames[picksp],all_neg_log_lik_0[picksp,pickmet],all_neg_log_lik_C[picksp,pickmet],
      all_neg_log_lik_D[picksp,pickmet],all_neg_log_lik_E[picksp,pickmet])

cbind(birdnames[picksp],all_lrt_0C[picksp,pickmet],all_lrt_CD[picksp,pickmet],
      all_lrt_CE[picksp,pickmet])

#PCA on basis of significance
PCA_lik_C <- princomp(all_lrt_0C[,1:vno] ,cor=TRUE)
summary(PCA_lik_C)
PCA_lik_C$scores
PCA_lik_C$loadings
par(mfrow=c(1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,mai=c(1,0.5,1.5,0.5)) # mai order b,l,t,r
par_labs4 <- c("i","a","e","p")
par(mfrow=c(1,1))
biplot(PCA_lik_C,xlabs=birdcodes,exp=0.9,ylabs=par_labs4,main="PCA_of_LRT_0_F1",xlab="Axis 1",ylab="Axis 2")

#Check values for contrasting species according to the biplot
all_lrt_0C[33,] #little owl
all_lrt_0C[18,] #greenfinch
all_lrt_0C[67,] #wren


#PCA using parameter estimates
PCA_data_C <- array(0,dim=c(lspno,(4*vno)))
PCA_data_C[ ,1:(4*vno)] <- all_par_C1_PCA[,1:(4*vno)]
PCA_data_C[ ,4] <- all_ann_decay_C1[ ,1] 
#PCA_data_C[ ,8] <- all_ann_decay_C1[ ,2] 
#PCA_data_C[ ,12] <- all_ann_decay_C1[ ,3] 
#PCA_data_C[ ,16] <- all_ann_decay_C1[ ,4] 
PCA_par_C <- princomp(PCA_data_C[,1:(4*vno)] ,cor=TRUE)
summary(PCA_par_C)
PCA_par_C$scores
PCA_par_C$loadings
par(mfrow=c(1,1))
#par_labs4 <- c("e2","e3","e4","e5")
#par_labs16 <- c("a2","a3","a4","a5","i2","i3","i4","i5",
#                "e2","e3","e4","e5","p2","p3","p4","p5" )


biplot(PCA_par_C)#,xlabs=birdcodes,exp=0.75) #,ylabs=par_labs4)
PCA_data_C[57,] #tree creeper
PCA_data_C[65,] #woodcock

PCA_lik_C [ ,1:4] <- princomp(PCA_data_C[,1:4] ,cor=TRUE)
all_neg_log_lik_C


# MAX AND MIN LRT VALUES

minOC <- min(all_lrt_0C)
minCD <- min(all_lrt_CD)
minCE <- min(all_lrt_CE)
cbind(minOC,minCD,minCE)

maxOC <- max(all_lrt_0C)
maxCD <- max(all_lrt_CD)
maxCE <- max(all_lrt_CE)
cbind(maxOC,maxCD,maxCE)


# ISOLATE NEGATIVE IMPROVEMENTS

for (imet in 1:vno){
  for (isp in 1:lspno){
    
    if (all_lrt_0C[isp,imet] < 0) {print(cbind("0C", isp, imet, all_lrt_0C[isp,imet] ))} 
    if (all_lrt_CD[isp,imet] < 0) {print(cbind("CD", isp, imet, all_lrt_CD[isp,imet] ))} 
    if (all_lrt_CE[isp,imet] < 0) {print(cbind("CE", isp, imet, all_lrt_CE[isp,imet] ))} 
  }}

cbind(birdnames[spno],birdcodes[spno],all_lrt_0C[,1:vno])


#Histograms of likelihood ratio tests
#Remove stray zeros

all_lrt_CDH<- array(0,dim=c(lspno,vno))
all_lrt_CDH <- all_lrt_CD
all_lrt_CDH[all_lrt_CDH<0] <- 0

all_lrt_CEH<- array(0,dim=c(lspno,vno))
all_lrt_CEH <- all_lrt_CE
all_lrt_CEH[all_lrt_CEH<0] <- 0

#give histograms in blocks of 4
par(mfrow=c(2,2))
lim0C <- c(0,5,10,15,20,25,30,35,40,45,50)
hist(all_lrt_0C[,1],xlim=c(0,50),ylim=c(0,40),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=lim0C,main="Min_Temp_0_F1",xlab="LRT_Chisq_4")
hist(all_lrt_0C[,2],xlim=c(0,50),ylim=c(0,40),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=lim0C,main="Max_Temp_0_F1",xlab="LRT_Chisq_4")
hist(all_lrt_0C[,3],xlim=c(0,50),ylim=c(0,40),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=lim0C,main="Mean_Temp_0_F1",xlab="LRT_Chisq_4")
hist(all_lrt_0C[,4],xlim=c(0,50),ylim=c(0,40),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=lim0C,main="Prec_0_F1",xlab="LRT_Chisq_4")

par(mfrow=c(2,2))
limCD <- c(-0.001,3,6,9,12,15,18,21,24,27,30)
hist(all_lrt_CDH[,1],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Min_Temp_F1_F2",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CDH[,2],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Max_Temp_F1_F2",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CDH[,3],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Mean_Temp_F1_F2",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CDH[,4],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Prec_F1_F2",xlab="LRT_Chisq_2",freq=TRUE)

par(mfrow=c(2,2))
limCD <- c(-0.001,3,6,9,12,15,18,21,24,27,30)
hist(all_lrt_CEH[,1],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Min_Temp_F1_F2h",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CEH[,2],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Max_Temp_F1_F2h",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CEH[,3],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Mean_Temp_F1_F2h",xlab="LRT_Chisq_2",freq=TRUE)
hist(all_lrt_CEH[,4],xlim=c(0,30),ylim=c(0,50),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,breaks=limCD,main="Prec_F1_F2h",xlab="LRT_Chisq_2",freq=TRUE)


#SCATTERPLOTS
pairs(~all_lrt_0C[,1]+all_lrt_0C[,2]+all_lrt_0C[,3]+all_lrt_0C[,4],main="LRT_0_F1_Scatterplot",labels=c("Min","Max","Mean","Prec"),
      xlim=c(0,40),ylim=c(0,40))
pairs(~all_lrt_CD[,1]+all_lrt_CD[,2]+all_lrt_CD[,3]+all_lrt_CD[,4],main="LRT_F1_F2_Scatterplot",labels=c("Min","Max","Mean","Prec"),  
      xlim=c(0,25),ylim=c(0,25) )
pairs(~all_corr_C1[,1]+all_corr_C1[,2]+all_corr_C1[,3]+all_corr_C1[,4],main="F1_AR1_Correlation_Scatterplot",labels=c("Min","Max","Mean","Prec"),xlim=c(0,1),ylim=c(0,1))
pairs(~all_ann_decay_C1[,1]+all_ann_decay_C1[,2]+all_ann_decay_C1[,3]+all_ann_decay_C1[,4],main="F1_Annual_Decay_Scatterplot",labels=c("Min","Max","Mean","Prec"))

maxlrt0C <- max(all_lrt_0C[,])
maxlrt0C
par(mfrow=c(2,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
plot(y=all_lrt_0C[,1],x=all_ann_decay_C1[,1],main="Min_Temp_0_F1",ylab="LRT_Chisq_4",xlab="Annual_scaling_F1",ylim=c(0,30))
plot(y=all_lrt_0C[,2],x=all_ann_decay_C1[,2],main="Max_Temp_0_F1",ylab="LRT_Chisq_4",xlab="Annual_scaling_F1",ylim=c(0,30))
plot(y=all_lrt_0C[,3],x=all_ann_decay_C1[,3],main="Mean_Temp_0_F1",ylab="LRT_Chisq_4",xlab="Annual_scaling_F1",ylim=c(0,30))
plot(y=all_lrt_0C[,4],x=all_ann_decay_C1[,4],main="Prec_0_F1",ylab="LRT_Chisq_4",xlab="Annual_scaling_F1",ylim=c(0,30))

maxlrtCD <- max(all_lrt_CD[,])
maxlrtCD
maxlrtCE <- max(all_lrt_CE[,])
maxlrtCE
#show the next 2 sets of graphs with the same y scale
par(mfrow=c(2,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
plot(y=all_lrt_CD[,1],x=all_ann_decay_D1[,1],main="Min_Temp_F1_F2",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2",ylim=c(0,25))
plot(y=all_lrt_CD[,2],x=all_ann_decay_D1[,2],main="Max_Temp_F1_F2",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2",ylim=c(0,25))
plot(y=all_lrt_CD[,3],x=all_ann_decay_D1[,3],main="Mean_Temp_F1_F2",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2",ylim=c(0,25))
plot(y=all_lrt_CD[,4],x=all_ann_decay_D1[,4],main="Prec_F1_F2",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2",ylim=c(0,25))

maxlrtCE <- max(all_lrt_CE[,])
maxlrtCE
par(mfrow=c(2,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
plot(y=all_lrt_CE[,1],x=all_ann_decay_E1[,1],main="Min_Temp_F1_F2h",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2h",ylim=c(0,25))
plot(y=all_lrt_CE[,2],x=all_ann_decay_E1[,2],main="Max_Temp_F1_F2h",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2h",ylim=c(0,25))
plot(y=all_lrt_CE[,3],x=all_ann_decay_E1[,3],main="Mean_Temp_F1_F2h",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2h",ylim=c(0,25))
plot(y=all_lrt_CE[,4],x=all_ann_decay_E1[,4],main="Prec_F1_F2h",ylab="LRT_Chisq_2",xlab="Annual_scaling_F2h",ylim=c(0,25))


#now do something on timing during year 1 of effects
#Check for wren (isp=67), min temp (imet=1)
all_par_C1[67,1, ] 
twopi <- 2*3.14159
month12 <- c(1:12)
beta12 <- ((1/(1+exp(-all_par_C1[67,1,6])))^month12)*
  (all_par_C1[67,1,3]+all_par_C1[67,1,4]*sin(twopi*month12/12)+all_par_C1[67,1,5]*cos(twopi*month12/12))

month60 <- c(1:60)
beta60 <- ((1/(1+exp(-all_par_C1[67,1,6])))^month60)*
  (all_par_C1[67,1,3]+all_par_C1[67,1,4]*sin(twopi*month60/12)+all_par_C1[67,1,5]*cos(twopi*month60/12))

par(mfrow=c(1,1))
plot(y=beta60,x=month60)

#roll out to all species, all weather variables 
###########################

all_beta_C <- array(0,dim=c(lspno,48))

#for (isp in (1:lspno)){
#  for (imet in (1:vno)) {
#    beta12 <- ((1/(1+exp(-all_par_C1[isp,imet,6])))^month12)*
#      (all_par_C1[isp,imet,3]+all_par_C1[isp,imet,4]*sin(twopi*month12/12)+all_par_C1[isp,imet,5]*cos(twopi*month12/12))  
#    all_beta_C[isp,(12*(imet-1)+(1:12))] <- beta12[1:12] 
#  }}

all_beta_C <- array(0,dim=c(lspno,12))

for (isp in (1:lspno)){
    beta12 <- ((1/(1+exp(-all_par_C1[isp,,6])))^month12)*
      (all_par_C1[isp,,3]+all_par_C1[isp,,4]*sin(twopi*month12/12)+all_par_C1[isp,,5]*cos(twopi*month12/12))  
    all_beta_C[isp,1:12] <- beta12[1:12] 
  }

all_beta_C[1,]
"check wren again"
plot(y=all_beta_C[67,1:12],x=month12)

"variable labels for PCA"
#par_labs24 <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12",
#                "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12")

par_labs12T <- c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
#par_labs12P <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12")
par_labs12T <- c("Jun","May","Apr","Mar","Feb","Jan","Dec","Nov","Oct","Sep","Aug","Jul")

#use mean temperature and precipitation
#PCA_beta_C_TP <- princomp(all_beta_C[,25:48] ,cor=TRUE)
#summary(PCA_beta_C_TP)
#PCA_beta_C_TP$scores
#PCA_beta_C_TP$loadings

#par(mfrow=c(1,1),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,mai=c(1,0.5,1.5,0.5)) # mai order b,l,t,r
#biplot(PCA_beta_C_TP,xlabs=birdcodes,exp=0.75,ylabs=par_labs24,main="PCA_prec_mean_temp_beta_F1",xlab="Axis 1",ylab="Axis 2")



#use mean temperature only
PCA_beta_C_T <- princomp(all_beta_C ,cor=TRUE)
summary(PCA_beta_C_T)
PCA_beta_C_T$scores
PCA_beta_C_T$loadings

library(car)
par(mfrow=c(1,1),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,mai=c(1.5,1.5,0.5,1)) # mai order b,l,t,r
#biplot(PCA_beta_C_T,xlabs=birdcodes[spno],exp=0.75,ylabs=par_labs12T,main="PCA_mean_temp_beta_F1",xlab="Axis 1",ylab="Axis 2")
plot(PCA_beta_C_T$scores[,1],PCA_beta_C_T$scores[,2],col=as.numeric(taxa)+1,xlab="Axis 1",ylab="Axis 2",pch=16,cex=0.5)
arrows(x0=0,y0=0,PCA_beta_C_T$loadings[,1]*20,PCA_beta_C_T$loadings[,2]*20)#,xlabs=birdcodes[spno],exp=0.75,ylabs=par_labs12T,main="PCA_mean_temp_beta_F1",xlab="Axis 1",ylab="Axis 2")
text(PCA_beta_C_T$loadings[,1]*20,PCA_beta_C_T$loadings[,2]*20,par_labs12T,cex=0.8)#,xlabs=birdcodes[spno],exp=0.75,ylabs=par_labs12T,main="PCA_mean_temp_beta_F1",xlab="Axis 1",ylab="Axis 2")
dataEllipse(PCA_beta_C_T$scores[,1],PCA_beta_C_T$scores[,2],groups=taxa,add=T,plot.points = F,levels = c(0.95),center.pch=F,col=c(2:6),cex=0.9)

"check four species"
"coal tit 9"
"lesser redpoll 28"
"whitethroat 62"
"wren, 67"

#graph with y axes equal
miny <- min(all_beta_C[9,1:12],all_beta_C[28,1:12],all_beta_C[62,1:12],all_beta_C[67,1:12])
maxy <- max(all_beta_C[9,1:12],all_beta_C[28,1:12],all_beta_C[62,1:12],all_beta_C[67,1:12])
cbind(miny,maxy)
par(mfrow=c(2,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
cylim <-c(-0.11,0.11)
plot(y=all_beta_C[9,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Coal_tit_F1_mean_temp" ,xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[28,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Lesser_redpoll_F1_mean_temp",xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[62,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Whitethroat_F1_mean_temp",xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[67,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Wren_F1_mean_temp",xlim=c(0,13),ylim=cylim)

#graph with y axes unconstrained
par(mfrow=c(2,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
plot(y=all_beta_C[9,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Coal_tit_F1_mean_temp" ,xlim=c(0,13))
plot(y=all_beta_C[28,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Lesser_redpoll_F1_mean_temp",xlim=c(0,13))
plot(y=all_beta_C[62,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Whitethroat_F1_mean_temp",xlim=c(0,13))
plot(y=all_beta_C[67,1:12],x=month12,ylab="Beta",xlab="Months before survey",main="Wren_F1_mean_temp",xlim=c(0,13))



#use precipitation only
PCA_beta_C_P <- princomp(all_beta_C[,37:48] ,cor=TRUE)
summary(PCA_beta_C_P)
PCA_beta_C_P$scores
PCA_beta_C_P$loadings


par(mfrow=c(1,1),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,mai=c(1,0.5,1.5,0.5)) # mai order b,l,t,r
biplot(PCA_beta_C_P,xlabs=birdcodes,exp=0.75,ylabs=par_labs12P,main="PCA_prec_beta_F1",xlab="Axis 1",ylab="Axis 2")

#"check four species"
#"bullfinch 4"
#"buzzard 5"
#"little owl, 33"
#"tufted duck 60" 

#graph with y axes equal
miny <- min(all_beta_C[4,37:48],all_beta_C[5,37:48],all_beta_C[33,37:48],all_beta_C[60,37:48])
maxy <- max(all_beta_C[4,37:48],all_beta_C[5,37:48],all_beta_C[33,37:48],all_beta_C[60,37:48])
cbind(miny,maxy)
cylim <-c(-0.12,0.12)
par(mfrow=c(2,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
plot(y=all_beta_C[4,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Bullfinch_F1_prec" ,xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[5,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Buzzard_F1_prec",xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[33,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Little_owl_F1_prec",xlim=c(0,13),ylim=cylim)
plot(y=all_beta_C[60,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Tufted_duck_F1_prec",xlim=c(0,13),ylim=cylim)


#graph with y axes unconstrained
par(mfrow=c(2,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
plot(y=all_beta_C[4,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Bullfinch_F1_prec" ,xlim=c(0,13))
plot(y=all_beta_C[5,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Buzzard_F1_prec",xlim=c(0,13))
plot(y=all_beta_C[33,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Little_owl_F1_prec",xlim=c(0,13))
plot(y=all_beta_C[60,37:48],x=month12,ylab="Beta",xlab="Months before survey",main="Tufted_duck_F1_prec",xlim=c(0,13))

######################################################

birdnames[spno]
par(mfrow=c(1,1))
boxplot(all_lrt_0C[taxa=="birds"],all_lrt_0C[taxa=="moths"],all_lrt_CD[taxa=="birds"],all_lrt_CD[taxa=="moths"],all_lrt_CE[taxa=="birds"],all_lrt_CE[taxa=="moths"])

boxplot(all_lrt_0C~taxa)
plot(all_final_Fnd_C[-1616])

maxmonth<-c();minmonth<-c()
for(i in 1:lspno){maxmonth[i]<-which.max(all_final_Fnd_C[i,])
                  minmonth[i]<-which.min(all_final_Fnd_C[i,])}

monthmeans<-rbind(colMeans(all_final_Fnd_C[taxa=="birds",]),colMeans(all_final_Fnd_C[taxa=="moths",]),
                  colMeans(all_final_Fnd_C[taxa=="aphids",]))

monthsd<-rbind(apply(all_final_Fnd_C[taxa=="birds",],2,sd),apply(all_final_Fnd_C[taxa=="moths",],2,sd),
               apply(all_final_Fnd_C[taxa=="aphids",],2,sd))

plot(monthmeans[1,],ylim=c(min(monthmeans)/5,max(monthmeans)/20),type="l")
lines(monthmeans[2,]/5,col=2)
lines(monthmeans[3,]/20,col=4)

lines(monthmeans[1,]+monthsd[1,])
lines(monthmeans[1,]-monthsd[1,])
lines(monthmeans[2,]/5+monthsd[2,]/5,col=2)
lines(monthmeans[2,]/5-monthsd[2,]/5,col=2)
lines(monthmeans[3,]+monthsd[2,],col=4)
lines(monthmeans[3,]-monthsd[2,],col=4)


plot(all_final_Fnd_C~taxa)
plot(monthmeans)  

library(circular)
circular()

prcomp(all_final_Fnd_C)
mine<-princomp(all_final_Fnd_C,cor=TRUE)
biplot(mine)

plot(mine$scores[,1],mine$scores[,2],type="n")
text(mine$scores[,1],mine$scores[,2],c(1:lspno),col=as.numeric(taxa))

names(mine$scores[,1])

#-- 
#  Biomathematics and Statistics Scotland (BioSS) is formally part of The
#James Hutton Institute (JHI), a registered Scottish charity No. SC041796
#and a company limited by guarantee No. SC374831