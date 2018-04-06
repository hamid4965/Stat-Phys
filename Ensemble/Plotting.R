## This is for ploting the results


Plotting <- function (en,bayes,pls_mean_cof, main="",type,legend) {

graphics.off()
wl <- en[[1]][1,]
fit <- en[[2]]
cf <- en[[1]][2:4,]
cf <- cf * fit
fit <- round (fit, 2)
encf <- en[[1]][5,]
sel <- en[[1]][6,]
sel[sel==0] <- NA
vip_Ref<-en$vipRef
vip <- 1*vip_Ref
vip[vip==0] <- NA
p_bayes <- unlist(bayes[1])   # this is the probabilities of bands selected
beta_bma <- unlist(bayes[2])  # this is the coefficients estimated using the Bayesian regression
p_bayes_final <- 1*(p_bayes > mean(p_bayes)+sd(p_bayes)) # select the bands that their Bayesian probability is more than mean+std

p_bayes_final[p_bayes_final==0] <-NA
maxy <- max(cf)
miny <- min(cf)
max_cf<-max(cf)

#-----------plot Ensemble-----------------------------------
par (mfrow = c(3,1),mar=c (0.5, 5, 2, 1.2) + 0.1, xpd=F)
plot(wl, cf[1,], col="black", lwd=3,ty="l",ylim=c(miny*1.2,maxy*1.8),ylab="Coefficient values",xlab="",
     cex.lab=1.5,cex.axis=1.5,xaxt='n')
#grid (NULL,NULL) 
title(type, outer=FALSE,cex.main=2)
abline(0,0)
lines (wl, cf[2,], col="blue", lwd=3)
lines (wl, cf[3,], col="brown", lwd=3)
points(wl,sel*miny*1.2, yaxt = "n",ann=FALSE, pch=8,cex.axis=2,cex=1,col="red")
text(max (wl)*0.5, max (cf)*1.5,"Ensemble",cex=1.5)
if (legend==T){
legend( max (wl)*0.60, max (cf)*1.8, 
        legend=c("PLS","RF","SVM","Selected bands"),
        col=c("black","blue","brown","red"), lwd=1, lty=c(1,1,1,NA), 
        pch=c(NA,NA,NA,8),bty = "n",cex=1.2,text.font=2,y.intersp=0.75)
}

#--------- Plot the bayesian----------------------------------
par (mar=c (1, 5, 1, 1.2) + 0.1, xpd=F)
plot (wl,beta_bma, col="green",ty="l", lwd=3,ylab="Coefficient values",xlab="",
      cex.lab=1.5,cex.axis=1.5,xaxt='n')
abline(0,0)
points(wl,p_bayes_final*min(beta_bma)*1.2, yaxt = "n",ann=FALSE, pch=8,cex.axis=2,cex=1,col="red")
text(max (wl)*0.5, max(beta_bma)*0.8,"BR",cex=1.5)
if (legend==T){
  legend( max (wl)*0.60, max (beta_bma)*1, 
          legend=c("Bayesian Coefficients","Selected bands"),
          col=c("green","red"), lwd=1, lty=c(1,NA), 
          pch=c(NA,8),bty = "n",cex=1.2,text.font=2,y.intersp=0.75)
}


#----------- Plot PLS model-------------------------------------
par (mar=c (4.1, 5, 0, 1.2) + 0.1, xpd=F)
plot(wl,pls_mean_cof,ylim=c(min(pls_mean_cof)*1.2,max(pls_mean_cof)*1.1),type="l",xlab="Wavelength (nm)", ylab="VIP values", 
     cex.lab=1.5, cex.axis=1.5)  # plot the mean and SD of the VIPs

abline(0,0)
#polygon(c(wl,rev(wl)),c(pls_mean_cof-Band_sd,rev(pls_mean_cof+Band_sd)),col = "grey75", border = FALSE)
lines(wl,pls_mean_cof, lwd = 2)
text(max (wl)*0.5, max (pls_mean_cof)*1.2,"PLS_VIP",cex=1.5)

points(wl,vip*min(pls_mean_cof)*1.2,ann=FALSE, pch=8,cex.axis=1.5,cex=1,col="red")
if (legend==T){
  legend( max (wl)*0.65, max (pls_mean_cof)*1.2, 
          legend=c("PLS_ref Coefficients","Selected bands"),
          col=c("black","red"), lwd=1, lty=c(1,NA), 
          pch=c(NA,8),bty = "n",cex=1.2,text.font=2,y.intersp=0.9)

  }


}


