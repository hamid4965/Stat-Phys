# StatPaper
# Paper Title: Estimation of vegetation nitrogen content: uncertainity in upscaling and feature selection. 
# Phases: 
  # phase 1; Data preparation (leaf, canopy and plot scales)
  # Phase 2; Run the PLS using the optimum number of components
  # Phase 3; Variable selection, predictivity analysis and variable overlap
  # Phase 4; Do upscaling using LAI
  # Phase 5 Variable selection again with new features


#########################################################
                                                        #
              #Phase 1 (data preparation)               #
# Import data from leaf, canopy and plots               #
                                                        #
#########################################################
# install.packages("openxlsx")
# install.packages("hsdar")
# install.packages("prospectr")
# install.packages("hyperSpec")
#install.packages("pls")
# install.packages("e1071")
# install.packages ("randomForest")




  rm(list=ls(all=TRUE))                     # clear workspace
  graphics.off()                            # close any open graphics
  setwd("N:/Data02/bcal/Personal/hamid/StatPaper/final_working/PNAS")        # set the directory 
  
  ###### Load library ####
  library("openxlsx")
  library("hsdar")
  library("pls")
  library("prospectr")
  library("hyperSpec")
  library("randomForest")
  library("reshape2")
  
  
  source("Data_Prep.R")
  #source("Data_Prep_v2.R")
  source("My_KS.R") 
source("Random_Selection.R")
  source("VIP.R")
  source("Pred_R2.R")
  source("specplot.R")
  #source("R2adj.R")
  source("rse09382-mmc2.R")
  source("mean_predict2.R")
  source("Plotting.R")
  

################ Reading  data and ploting ###################
#Type: shrub, leaf, plot, dryleaf_INL, dryleaf_ACCP, OPT, lopex, "DryLeaf"
# transformation : "log", "derv1", "raw","smooth", "cr", "derv1log"
# atr : PercentN; LAI for shrubs and LAI_plot(Inds) for plots
Data<-Data_Prep(Type = "shrub",atr = "LAI",outlier = FALSE,transformation = "derv1",resampling =F) 
# Preparing data for the Hsdar and ploting based on different subsets
NewSpecLib<-speclib(Data$Spectra,Data$Wavebands)
#SI(NewSpecLib) <- Data$other

plot(NewSpecLib, FUN =  (1:dim(Data$Spectra)[1]), main = "Spectrum") # plot all the bands based on the Hsdar package  (1:dim(Data$Spectra)[1])

################ Extracting data for canopy specral test#####################
#atr<-Data$other

#SP<-Data$Spectra       # Shrubs for N-sage17 plot
#plot(Data$Wavebands, SP)

#write.table(SP, file="plot_spectra.txt", row.names=FALSE, col.names=FALSE,sep = "\t")
#write.table(Data$Wavebands, file="plot_wave.txt", row.names=FALSE, col.names=FALSE)
#write.table(atr, file="plot_atr.txt", row.names=FALSE, col.names=colnames(atr))




#############################################################################
                                                                            #
      # Phase 2 Run the PLS using optimum number of PC and do the CV        #
                                                                            #
#############################################################################
# Phase 2: prediction capability of the PLSR at different scales
# perturb the original data
data <- as.matrix(cbind(Data$atr,Data$Spectra))
set.seed(09)    #Setting seeds starts with 03 for leaf smooth, then leaf derv leaf log.....
seed<-09
shuffled_data <- data[sample(nrow(data)),]
X<-shuffled_data[,2:ncol(shuffled_data)]
Y<-shuffled_data[,1]
dframe_all <- data.frame(N=I(Y), spectra=I(X))
wavesbands <- as.numeric(colnames(X))
# test_Ratio <- 0.1
# ratio <-ceiling(test_Ratio*nrow(X))

# Devide the target dataset(leaf,shrub, plot) scale dataset into Cal, test based on the KS algorithm.
# KS_Rank <- My_KS(X)
# N_test <- Y[KS_Rank[1:ratio]]
# Spectra_test <- X[KS_Rank[1:ratio],]
# dframe_test <- data.frame(N=I(N_test), spectra=I(Spectra_test))
# N_cal <- Y[KS_Rank[(ratio+1):nrow(X)]]
# Spectra_cal <- X[KS_Rank[(ratio+1):nrow(X)],]
# dframe_cal <- data.frame(N=I(N_cal),spectra=I(Spectra_cal))


# Repeated double CV
itr <- 1000
In_Ratio=0.2
ncomp <- 15
ncomp.onesigma<-vector(length = itr)
R2_plain <-vector(length = itr)
intercept <-vector(length = itr)
VIPs <- matrix(nrow = ncol(dframe_all$spectra), ncol = itr)
RMSEP_Cal <- vector(length = itr)
Reg_Coef <- matrix(nrow= (ncol(dframe_all$spectra)+1), ncol = itr)
pls_scale <- matrix(nrow = (ncol(dframe_all$spectra)),ncol=itr)
loading_cof <- matrix(nrow = ncol(dframe_all$spectra), ncol = itr)                    

for (i in 1:itr){
  
       newList <- Random_Selection(dframe_all,In_Ratio)
       In_df_cal<-data.frame(N=I(newList$Ycal),spectra=I(newList$Xcal))
       In_df_val<-data.frame(N=I(newList$Yval),spectra=I(newList$Xval))
      
      Train_pls<-plsr(N ~ spectra,scale=T, ncomp = ncomp, data = In_df_cal,method = "oscorespls", validation = "LOO" )
      ncomp.onesigma[i]<- selectNcomp(Train_pls, method = "onesigma", plot = F, ylim = c(.18, .6))
        if (ncomp.onesigma[i]==0){
          next
        }
      pred<- drop(predict(Train_pls,ncomp= ncomp.onesigma[i], newdata=In_df_val))
      R2_plain[i]<-caret::R2(pred,In_df_val$N)
      rms<-RMSEP(Train_pls,newdata = In_df_val,ncomp = ncomp.onesigma[i])
      RMSEP_Cal[i] <- (rms$val[2]*100)/mean(In_df_val$N)
      S1<-VIP(Train_pls)
      VIPs[,i] <- S1[ncomp.onesigma[i],]
      Reg_Coef[,i] <- coef(Train_pls, ncomp = ncomp.onesigma[i], intercept = TRUE)
      loading_cof[,i] <- Train_pls$loading.weights[,ncomp.onesigma[i]]
      pls_scale[,i] <-Train_pls$scale
      base::print(i)
}

# Following is the final model is based on the number of components selected fromt he repeated double Cv. 

#final_comp <-as.numeric(tail(names(sort(table(ncomp.onesigma))),1))
#if (final_comp==0){ final_comp<- as.numeric(tail(names(sort(table(ncomp.onesigma))),2))[1]}
#final_pls <- plsr(N ~ spectra, scale=T, ncomp = ncomp, data = dframe_cal,method = "oscorespls", validation = "LOO" )
#Rccv_pred <- drop(predict(final_pls,ncomp = final_comp, newdata=dframe_test ))
#final_plain_R2<-caret::R2(Rccv_pred,dframe_test$N)
#RMSEPs_Rdcv <- (sqrt((sum((Rccv_pred - dframe_test$N)^2))/length(dframe_test$N))/mean(dframe_test$N))*100

# Based on the mean mnodel 
# B<- base::rowMeans(Reg_Coef,na.rm = T)
# Mean_pred<-mean_predict(B,unclass(dframe_test$spectra),pls_scale)
# final_plain_R2_Mean<-drop(caret::R2(Mean_pred,dframe_test$N))
# RMSEPs_mean <- sqrt((sum((Mean_pred - dframe_test$N)^2))/length(dframe_test$N))/mean(dframe_test$N)*100
# #final_adj_R2_Mean<-R2adj(length(dframe_test$N),final_comp,final_plain_R2_Mean)
# final_Pred_R2_Mean <- Pred_R2(dframe_test$N,Mean_pred)
# loadings_final<-base::rowMeans(loading_cof,na.rm = T)


#-----------------------------------------------------------------------------------------------------------------------
# Phase 3: Variable selection plsRef

#Based on the bootstraping method suggested in (http://www.sciencedirect.com/science/article/pii/S0169743909001841)
Band_Mean<-rowMeans(VIPs,na.rm = T)
Band_sd <- apply(VIPs,1,sd,na.rm=T) 
VIP_Mean <-  Band_Mean+Band_sd >1 & Band_Mean-Band_sd >1 

plot(wavesbands,Band_Mean,type="l",xlab="Wavelength (nm)", ylab="VIP value", cex.lab=1.5)  # plot the mean and SD of the VIPs
polygon(c(wavesbands,rev(wavesbands)),c(Band_Mean-Band_sd,rev(Band_Mean+Band_sd)),col = "grey75", border = FALSE)
lines(wavesbands,Band_Mean, lwd = 2)
abline(1,0)


#------------------------------------------
# Fit based on Asner method " Multi-method ensemble selection of spectral bands related to leaf biochemistry" (Ensemble method) 
itr2 <- 10
r2_pls_plain <- vector(length = itr2)

RMSEPs_pls <- vector(length = itr2)
r2_rf_plain <- vector(length = itr2)

RMSEPs_rf <- vector(length = itr2)
r2_svm_plain <- vector(length = itr2)

RMSEPs_svm <- vector(length = itr2)

flds <- createFolds(Y, k = 10, list = TRUE, returnTrain = FALSE)
se<-seq(from=1, to=length(Y),by=1)

for (j in 1:itr2){
  
  indexx<-flds[[j]]
  y_cal<-Y[!se %in% indexx ]
  x_cal<-X[!se %in% indexx ,]
  y_val<-Y[se %in% indexx ]
  x_val<-X[se %in% indexx ,]
  en <- ensemble(x_cal,y_cal,wavesbands)  
  
  en_val_df <- data.frame(cbind(y_val,x_val))
  ncomp_pls<- selectNcomp(en$PLS, method = "onesigma", plot = FALSE, ylim = c(.18, .6))
  pls_pred<- drop(predict(en$PLS,ncomp = ncomp_pls, newdata=en_val_df))
  r2_pls_plain[j]<-caret::R2(pls_pred,y_val)
  RMSEPs_pls[j] <- (sqrt((sum((y_val - pls_pred)^2))/length(pls_pred))*100)/mean(y_val)
  
  rf_pred <- predict(en$RF,x_val)
  r2_rf_plain[j]<-caret::R2(rf_pred,y_val)
  RMSEPs_rf[j] <- (sqrt((sum((y_val - rf_pred)^2))/length(rf_pred))*100)/mean(y_val)
  
  svm_pred <- predict(en$SVM,x_val)
  r2_svm_plain[j]<-caret::R2(svm_pred,y_val)
  RMSEPs_svm[j] <- (sqrt((sum((y_val - svm_pred)^2))/length(svm_pred))*100)/mean(y_val)
  
  base::print(j)
}

# use the ensemble function with all the data and compare the graphs with the mean of 1000

en_all <- ensemble(X,Y,wavesbands)
en_all$vipRef<-VIP_Mean
# ------------------------------------------------------------------------
save.image('shrub_LAI_derv1.rds')

# Read matlab Bayesian selected bands for plotting
rm(list=ls(all=TRUE))                     # clear workspace
graphics.off()   
Bayes_sel <- unlist(readMat("SelectedBands_shrub_LAI_log_bayes.mat"))
load("Corrected_canopy_LAI_Log.rds")
mean(R2_plain)
mean(RMSEP_Cal)
final_plain_R2_Mean
RMSEPs_mean
mean(r2_pls_plain)
mean(r2_svm_plain)
mean(r2_rf_plain)
mean(RMSEPs_pls)
mean(RMSEPs_svm)
mean(RMSEPs_rf)

source("Plotting.R")
Plotting(en_all,Bayes_sel,Band_Mean,Band_sd,type = "canopy-LAI-log",legend = F)







#-----------------------------------------------------
#Graphing the results
# Significant variables 
#Plotting(en_all,type = "Plot",legend = F)
#My.plot.ensemble(en_all,legends = F,type="Canopy")    #All observations in ensemble function
#save.image('canopy_N_Dervlog.rds')
#-----------------------------------------------------------------
#Exporting the selected wavelenghts 
a<-wavesbands[VIP_Mean]
b<-wavesbands[as.logical(en_all$selection[6,])]
save(b, file = "plot_band_selected_DervLog.RData")


load("Leaf_band_selected_DervLog.RData")
a4<-b
Reduce(intersect, list(a1,a3,a4))

#c<-wavesbands[VIP_Loocv>1]
Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip") 
write.xlsx(b,file="test.xlsx")
#-----------------------------------------------------------------


mean(R2_plain)
mean(RMSEP_Cal)
final_plain_R2_Mean
RMSEPs_mean
mean(r2_pls_plain)
mean(r2_svm_plain)
mean(r2_rf_plain)
mean(RMSEPs_pls)
mean(RMSEPs_svm)
mean(RMSEPs_rf)


#-----------------------------------------------------------------
#Plot the R2 between each band and the N
r2_bands<- vector(length = length(Data$Wavebands))
for (k in 1:dim(Data$Spectra)[2]){
mod1 <- lm(Data$atr ~ Data$Spectra[,k])
if (mod1$coefficients[2]<0){
  r2_bands[k]<- -1*caret::R2(Data$atr,mod1$fitted.values)
  } else {
    r2_bands[k]<- caret::R2(Data$atr,mod1$fitted.values)
  }
}
plot(Data$Wavebands, r2_bands)
#------------------------------------------------------------------
## calculating the Red edge and ploting it
rd<-rededge(NewSpecLib)
REIP<-median(rd[,3],na.rm = T)
NIR <- Data$Spectra[,100]
temp_df<-data.frame(x=Data$atr,y=NIR)
par (mar=c (4.1, 4.6, 1.2, 1.2) + 0.1)
with(temp_df,plot(x,y,pch=1,col="blue",ylab=paste("Canopy BRF (",round(REIP),"nm)"),
                  xlab="LAI",cex.lab=1.6),xlim=c(2,3))
abline(fit<-lm(y~x,data=temp_df),col='red')
legend("top", bty="n",cex=1.5)

# Based on the Ollinger 2005 method

pls_Loocv<-plsr(N ~ spectra,scale=T, ncomp = ncomp, data = dframe_all,method = "oscorespls", validation = "LOO" )
ncomp_pls<- selectNcomp(pls_Loocv, method = "onesigma", plot = F, ylim = c(.18, .6))
VIP_Loocv<-VIP(pls_Loocv)
pred_Loocv<-pls_Loocv$validation$pred[,1,ncomp_pls]
r2_pls_Loocv<-caret::R2(pred_Loocv,Y)
RMSEP_Loocv <- (sqrt((sum((Y - pred_Loocv)^2))/length(pred_Loocv))*100)/mean(Y)
plot(wavesbands,pls_Loocv$loadings[,ncomp_pls])
plot(wavesbands,VIP_Loocv[ncomp_pls,],type="l")  # plot the mean and SD of the VIPs
lines(wavesbands,Band_Mean, lwd = 2)
abline(1,0)

VIP_Loocv_sel<-(VIP_Loocv[ncomp_pls,] > 1)
VIP_Loocv_final<-wavesbands[VIP_Loocv_sel]
en_all$VIP_Loocv_sel<-(VIP_Loocv[ncomp_pls,] > 1)


#-----------------------------------------------------------------
# ploting the LMA vs N (load the Dryleaf dataset)
a1<-Data$other
yy<-as.numeric(a1[,11])
xx<-Data$atr
rr2<-caret::R2(xx,yy,na.rm = T)
rr2<-round(rr2,digits = 2)
par (mar=c (4.1, 4.6, 1.2, 1.2) + 0.1)
plot(xx,yy, xlab="Nitrogen (%)", ylab="Leaf mass per area (g/cm2)",cex=1.2,cex.lab=1.4)
abline(fit<-lm(yy~xx),col="red")
legend("topright", bty="n",cex=1.5, legend=paste("R2 =", rr2))
#--------------------------------------------------------------------
# ploting the N_Height (Load the shrub dataset)
a1<-Data$other
yy<-as.numeric(a1[,7])
height<-Data$atr
rr2<-caret::R2(height,yy,na.rm = T)
rr2<-round(rr2,digits = 2)
par (mar=c (4.1, 4.6, 1.2, 1.2) + 0.1)
plot(height,yy, xlab="Nitrogen (%)", ylab="Height (m)",cex=1.2,cex.lab=1.4)
abline(fit<-lm(yy~height),col="red")
legend("topright", bty="n",cex=1.5, legend=paste("R2 =", rr2))
#--------------------------------------------------------------------
# ploting the LAI_N (Load the shrub dataset)
a1<-Data$other
N<-as.numeric(a1[,7])
LAI<-as.numeric(a1[,9])
rr2<-caret::R2(N,LAI,na.rm = T)
rr2<-round(rr2,digits = 2)
par (mar=c (4.1, 4.6, 1.2, 1.2) + 0.1)
plot(N,LAI, xlab="Nitrogen (%)", ylab="LAI",cex=1.2,cex.lab=1.4)
abline(fit<-lm(LAI~N),col="red")
legend("topright", bty="n",cex=1.5, legend=paste("R2 =", rr2))



















