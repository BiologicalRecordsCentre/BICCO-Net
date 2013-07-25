#Fake weather data:

fakeweather<-read.table("c:\\Users\\blaise\\desktop\\BICCO-Net\\Fine Scale\\fake_weather.txt",header=T)
attach(fakeweather) 

climate<-fakeweather

##########

climate1<-cbind(climate[2:51,7],climate[2:51,6],climate[2:51,5],climate[2:51,4],climate[2:51,3],climate[2:51,2],climate[1:50,13],climate[1:50,12],climate[1:50,11],climate[1:50,10], climate[1:50,9],climate[1:50,8],climate[2:51,19],climate[2:51,18],climate[2:51,17],climate[2:51,16], climate[2:51,15],climate[2:51,14],climate[1:50,25], climate[1:50,24],climate[1:50,23],climate[1:50,22],climate[1:50,21],climate[1:50,20])
climate2<-cbind(year[6:51],climate1[5:50,1:12],climate1[4:49,1:12],climate1[3:48,1:12],climate1[2:47,1:12],climate1[1:46,1:12],climate1[5:50,13:24],climate1[4:49,13:24],climate1[3:48,13:24],climate1[2:47,13:24],climate1[1:46,13:24])                

transdata<-read.table("c:\\Users\\blaise\\desktop\\BICCO-Net\\Fine Scale\\transdata.txt",header=T)
spCovariates<-read.table("c:\\Users\\blaise\\desktop\\BICCO-Net\\Fine Scale\\spCovariates.txt",header=T)


climate3<-matrix(0,46,120)
for(i in 1:46){ for(j in 1:120){
  climate3[i,j]<-(climate2[i,j+1]-transdata[1,j+1])/transdata[2,j+1]}}

## Example with blackbirds:

sp<-1
spname<-spCovariates[sp,1];spname

weatherCovs<-matrix(0,46,120)
for(i in 1:46){ for(j in 1:120){
  weatherCovs[i,j]<-as.matrix(climate3[i,j]*spCovariates[sp,j+1])}}

annualweather<-cbind(c(1966:2011),rowMeans(weatherCovs[,1:60]),rowMeans(weatherCovs[,61:120]))
colnames(annualweather)<-c("year","temp","precipitation")
