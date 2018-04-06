
# first import the soil
n<-4
Data_soil <- read.xlsx("Data.xlsx", sheet = n, colNames = TRUE)
refCol<-2
nr<-dim(Data_soil)[1]
nc<-2151
Reflectance <-as.matrix(Data_soil[1:nr,refCol:nc])
wav <- (as.numeric(colnames(Reflectance)) >= 400 & as.numeric(colnames(Reflectance)) <= 2500)
Reflectance <- Reflectance[,wav]
wav <- (as.numeric(colnames(Reflectance)) >= 1300 & as.numeric(colnames(Reflectance)) <= 1450)
Reflectance<-Reflectance[,!wav]
wav <- (as.numeric(colnames(Reflectance)) >= 1750 & as.numeric(colnames(Reflectance)) <= 2000)
SoilReflectance<-Reflectance[,!wav]
Wavelength <- as.numeric(colnames(Reflectance))

plot(Wavelength,Reflectance[1,],ty="l")
for (i in 2:nr){
  lines(Wavelength,Reflectance[i,],ty="l")
}
  

Smoothed_Wavelength_plot <- read.xlsx("Data.xlsx", sheet = 11, colNames = FALSE)
Smoothed_Wavelength_plot <- as.numeric(Smoothed_Wavelength_plot[,])

wav<-as.numeric(colnames(Reflectance))
Spectra <- resample(Reflectance,wav,Smoothed_Wavelength_plot, "spline")
smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
trans_Spectra <- smoothing
trans_Wavelength <- as.numeric(colnames(trans_Spectra))
soil<-colMeans(trans_Spectra)

###################
Data_leaf<-Data_Prep(Type = "DryLeaf",atr = "PercentN",outlier = FALSE,transformation = "smooth",resampling =F) 
mean_leaf<-colMeans(Data_leaf$Spectra)
Data_shrub<-Data_Prep(Type = "shrub",atr = "LAI",outlier = FALSE,transformation = "smooth",resampling =F)
mean_shrub<-colMeans(Data_shrub$Spectra)
Data_plot<-Data_Prep(Type = "plot",atr = "PercentN",outlier = FALSE,transformation = "smooth",resampling =F)
mean_plot<-colMeans(Data_plot$Spectra)
par (mar=c (4.1, 4.6, 1.2, 1.2) + 0.1)
plot(Data_leaf$Wavebands,mean_leaf,type="l",lty=4,lwd=2.5,ylim=c(0.01,0.7),xlab="Wavelength (nm)",ylab="Reflectance",cex.lab=1.4,col="blue")
lines(Data_shrub$Wavebands,mean_shrub,type="l",lty=3,lwd=2.5,col="red")
lines(Data_plot$Wavebands,mean_plot,type="l",lty=2,lwd=2.5,col="brown")
lines(trans_Wavelength,soil,type="l",lty=1,lwd=2.5,"black")

legend(1800,0.65,lty=c(4,3,2,1),c("Dry leaf","Canopy","Plot","Soil"),lwd=c(2.5,2.5,2.5,2.5),col=c("blue","red","brown","black"),cex=1.1)
