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
  setwd("C:/PNAS")       # set the directory 
  
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
Data<-Data_Prep(Type = "plot",atr = "PercentN",outlier = FALSE,transformation = "log",resampling =F) 
# Preparing data for the Hsdar and ploting based on different subsets
NewSpecLib<-speclib(Data$Spectra,Data$Wavebands)
#SI(NewSpecLib) <- Data$other
#names(SI(NewSpecLib))
# 
# RC <- subset(NewSpecLib, SiteName == "RCEW")
# RC14 <- subset(RC, Year == "2014")   #20 samples
# RC15 <- subset(RC, Year == "2015")   #19 samples
# Hol14 <- subset(NewSpecLib, SiteName == "Hollister")   #15samples
# 
# BP <- subset(NewSpecLib, SiteName == "BP")
# BP14 <- subset(BP, Year == "2014")    #22 samples
# BP15 <- subset(BP, Year == "2015")    #21 samples
# 
# LP <- subset(NewSpecLib, SiteName == "LP")
# LP14 <- subset(LP, Year == "2014")    #10 samples
# LP15 <- subset(LP, Year == "2015")    #11 samples
plot(NewSpecLib, FUN =  (1:dim(Data$Spectra)[1]), main = "Spectrum") # plot all the bands based on the Hsdar package  (1:dim(Data$Spectra)[1])
################ Extracting data for cross site study #####################
data <- as.matrix(cbind(Data$atr,Data$Spectra))
atr<-Data$other
#SP<-Data$Spectra       # Shrubs for N-sage17 plot
#plot(Data$Wavebands, SP)

IND <- atr$SiteName=="RCEW" & atr$Year=="2014"
test_site <- data[IND,]
cal_data <- data[!IND,]

#############################################################################
                                                                            #
      # Phase 2 Run the PLS using optimum number of PC and do the CV        #
                                                                            #
#############################################################################
# Phase 2: prediction capability of the PLSR at different scales
# perturb the original data

set.seed(05)    #Setting seeds starts with 03 for leaf smooth, then leaf derv leaf log.....
shuffled_data <- cal_data[sample(nrow(cal_data)),]
X<-shuffled_data[,2:ncol(shuffled_data)]
Y<-shuffled_data[,1]
dframe_cal <- data.frame(N=I(Y), spectra=I(X))

dframe_test <- data.frame(N=I(test_site[,1]), spectra=I(test_site[,2:dim(test_site)[2]]))
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
VIPs <- matrix(nrow = ncol(dframe_cal$spectra), ncol = itr)
RMSEP_Cal <- vector(length = itr)
Reg_Coef <- matrix(nrow= (ncol(dframe_cal$spectra)+1), ncol = itr)
pls_scale <- matrix(nrow = (ncol(dframe_cal$spectra)),ncol=itr)
loading_cof <- matrix(nrow = ncol(dframe_cal$spectra), ncol = itr)                    

for (i in 1:itr){
  
       newList <- Random_Selection(dframe_cal,In_Ratio)
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
 B<- base::rowMeans(Reg_Coef,na.rm = T)
 Mean_pred<-mean_predict(B,unclass(dframe_test$spectra),pls_scale)
 final_plain_R2_Mean<-drop(caret::R2(Mean_pred,dframe_test$N))
 RMSEPs_mean <- sqrt((sum((Mean_pred - dframe_test$N)^2))/length(dframe_test$N))/mean(dframe_test$N)*100
# loadings_final<-base::rowMeans(loading_cof,na.rm = T)


#-----------------------------------------------------------------------------------------------------------------------
# Phase 3: Variable selection plsRef

#Based on the bootstraping method suggested in (http://www.sciencedirect.com/science/article/pii/S0169743909001841)
# Band_Mean<-rowMeans(VIPs,na.rm = T)
# Band_sd <- apply(VIPs,1,sd,na.rm=T) 
# VIP_Mean <-  Band_Mean+Band_sd >1 & Band_Mean-Band_sd >1 
# 
# plot(wavesbands,Band_Mean,type="l",xlab="Wavelength (nm)", ylab="VIP value", cex.lab=1.5)  # plot the mean and SD of the VIPs
# polygon(c(wavesbands,rev(wavesbands)),c(Band_Mean-Band_sd,rev(Band_Mean+Band_sd)),col = "grey75", border = FALSE)
# lines(wavesbands,Band_Mean, lwd = 2)
# abline(1,0)


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

en_cal_all <- ensemble(X,Y,wavesbands)
ncomp_pls<- selectNcomp(en_cal_all$PLS, method = "onesigma", plot = FALSE, ylim = c(.18, .6))
en_test_dframe <- data.frame(cbind(test_site[,1],test_site[,2:dim(test_site)[2]]))

pls_pred<- drop(predict(en_cal_all$PLS,ncomp = ncomp_pls, newdata=en_test_dframe))
r2_pls_final<-caret::R2(pls_pred,test_site[,1])
RMSEPs_pls_final <- (sqrt((sum((test_site[,1] - pls_pred)^2))/length(pls_pred))*100)/mean(test_site[,1])

rf_pred <- predict(en$RF,test_site[,2:dim(test_site)[2]])
r2_rf_final<-caret::R2(rf_pred,test_site[,1])
RMSEPs_rf_final <- (sqrt((sum((test_site[,1] - rf_pred)^2))/length(rf_pred))*100)/mean(test_site[,1])

svm_pred <- predict(en$SVM,test_site[,2:dim(test_site)[2]])
r2_svm_final<-caret::R2(svm_pred,test_site[,1])
RMSEPs_svm_final <- (sqrt((sum((test_site[,1] - svm_pred)^2))/length(svm_pred))*100)/mean(test_site[,1])


save.image('Cross_N_smooth_LP15.rds')





