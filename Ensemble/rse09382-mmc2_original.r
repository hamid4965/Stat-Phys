################################################################################
##                                                                            ##
## Multi-method ensemble selection of spectral bands                          ##
##                                                                            ##
## This function performs a band selection based on a multi-method ensemble   ##
## assessment of the variable importance and regression coefficients of three ##
## different model types: Partial Least Squares regression, Random Forest     ##
## regression, and Support Vector Machine regression. Interactive user input  ##
## is required regarding the number of latent vectors included in the PLSR-   ##
## model.                                                                     ##
##                                                                            ##
## The script was tested with R version 3.1.2 on a x86_64-pc-linux-gnu        ##
## (64-bit) platform using the packages pls 2.4-3, randomForest 4.6-7, and    ##
## e1071 1.6-3. The script is not maintained and comes without any warranty.  ##
## In no event shall the authors be liable for any damages arising in any way ##
## out of the use of this code.                                               ##
##                                                                            ##
## Arguments:                                                                 ##
## x        Numeric matrix containing the spectra (samples as rows)           ##
## y        Numeric vector containing the response variable                   ##
## wl       Numeric vector containing the wavelength information of the bands ##
##                                                                            ##
## Contact: hannes.feilhauer@fau.de                                           ##
##                                                                            ##
################################################################################

ensemble <- function (x, y, wl=NA) {

  if (is.na (wl[1]) | length (wl)!=ncol (x))
     wl <- 1:ncol (x)
  
  ## load required libraries
  library (e1071)
  library (randomForest)
  library (pls)
  
  ## Model 1: Partial Least Squares regression
  dat <- data.frame (cbind (y, x)) 
  colnames (dat)[1] <- "y"
  pls <- plsr (y~., data=dat, jackknife=TRUE, scale=T, validation="CV")
  ## PRESS statistic for LV selection
  rmse <- RMSEP (pls, "all") 
  rms <- rmse$val[1:2, , -1][,1:25]
  barplot (rms[2,], space=0, col=4, ylim=c (0, max (rms[2,])), las=2, 
           cex.names=0.6, xlab="# latent vectors", ylab="RMSE")
  barplot (rms[1,], space=0, col=2, add=T, axes=F, axisnames=F)
  legend ("topright", legend=c ("cal", "val"), pch=15, col=c (2, 4), bty="n",  
           horiz=T)
  nlv <- as.numeric (readline ("Select # latent vectors for PLSR:"))
  plscf <- as.vector (coef (pls, ncomp=nlv, intercept=F)) ## extract coeff.
  plscf <- plscf / sd (plscf) ## scale regression coefficients
  plsrsq <- R2 (pls, "CV")$val[1, , -1][nlv] ## extract model fit 
  
  ## Model 2: Random Forest regression
  rf <- randomForest (x, y, importance=T)
    rfrsq <- rf$rsq[500] ## extract model fit for weighting
    rfcf <- rf$importance[,1] ## extract variable imprortances
    rfcf <- as.vector (rfcf / sd (rfcf)) ## scale variable importances
  
  ## Model 3: Support Vector Machine regression
  ## Grid search for C and gamma parameterization
  ccoef <- 2^seq (-5, 15, 2) ## vector of tested C coefficients
  gcoef <- 2^seq (-15, 3, 2) ## vector of tested gamma coefficients
  gsfit <- matrix (0, 11, 10) ## storage
  for (i in 1:11){
    for (j in 1:10){
      svmrun <- svm (x, y, type="eps-regression", cross=10, gamma=gcoef[j], 
                  cost=ccoef[i])
      gsfit[i, j]<-svmrun$scorrcoef
    }} ## run svr for all possible combinations
  svm.param <- which (gsfit==max (gsfit), arr.ind=TRUE)[1,] ## identify best 
  svr <- svm (x, y, type="eps-regression", cross=10, gamma=gcoef[svm.param[2]], 
              cost=ccoef[svm.param[1]]) ## final model
    svr.alpha <- t (svr$coefs) ## extract alpha vector
    svr.index <- svr$index ## extract alpha index
    ## calculate pseudo-regression coefficients from the alpha vector
        svrcf <- numeric (ncol (x))
        for(i in 1:ncol(x)) 
            svrcf[i] <- svr.alpha %*% x[svr.index, i]
    svrcf <- svrcf / sd (svrcf) ## scale pseudo-coefficients
    svrrsq <- svr$scorrcoef ## extract model fit for weighting
    
  ## get ensemble from all models and identify important variables
  ensemblecf <- abs (plscf) * plsrsq + abs (svrcf) * svrrsq + abs (rfcf) * rfrsq
  th <- mean (ensemblecf) + sd (ensemblecf) ## calculate threshold
  selbands <- ensemblecf > th ## apply threshold
  
  ## prepare output
  cf <- rbind (wl, plscf, rfcf, svrcf, ensemblecf, selbands)
    colnames (cf) <- colnames (x)
  
  fit <- c (plsrsq, rfrsq, svrrsq)
    names (fit) <- c ("PLS R2", "RF R2", "SVR R2")
  output <- list (cf, fit, th, pls, rf, svr)
    names (output) <- c ("selection", "fits", "threshold", "PLS", "RF", "SVM")
    class (output) <- "ensemble"
    output
}  

################################################################################
############################# END function ensemble ############################
################################################################################

################################################################################
##                                                                            ##
## plot.ensemble: visualization of ensemble objects                           ##
##                                                                            ##
## Arguments:                                                                 ##
## en       Ensemble object; output of the function ensemble()                ##
## main     Title for the graph                                               ##
##                                                                            ##
## Contact: hannes.feilhauer@fau.de                                           ##
##                                                                            ##
################################################################################

plot.ensemble <- function (en, main="Ensemble selection") {
  wl <- en[[1]][1,]
  fit <- en[[2]]
  cf <- en[[1]][2:4,]
  cf <- cf * fit
  fit <- round (fit, 2)
  encf <- en[[1]][5,]
  z1 <- matrix (rep (encf, 100), ncol=100)
  sel <- as.logical (en[[1]][6,])
  z2 <- matrix (rep (sel, 100), ncol=100)
      z2[,11:100] <- NA
      z2[z2==0] <- NA
  par (mar=c (5, 4, 4, 7) + 0.1, xpd=NA)
  image (wl, seq (min (cf) * 1.1, max (cf) * 1.1, length.out=100), z1, las=1, 
         xlab="Wavelength / nm", ylab="Weighted coefficients",  main=main,
         col=gray.colors (100, 1, 0, 1), cex.axis=0.7, cex.lab=0.7)
  image (wl, seq (min (cf) * 1.1, max (cf) * 1.1, length.out=100), z2, col=7, 
         xlab="", ylab="", add=T)
  lines (wl, cf[1,], col=2, lwd=2)
  lines (wl, cf[2,], col=4, lwd=2)
  lines (wl, cf[3,], col=3, lwd=2)
  
  labels <- c (paste (c ("PLS ", "RF ", "SVM "), "R2", "=", fit, sep=""), NA, 
               "ENSEMBLE", "Important", "No importance", "Selected")
  legend (max (wl)*1.05, max (cf), bty="n", col=c (2, 4, 3, NA, NA, rep (1, 3)), 
          pt.bg=c (rep (NA, 5), 1, "white", 7), lwd=c (rep (2, 3), rep (NA, 5)),
          pch=c (rep (NA, 5), rep (22, 3)), cex=0.7, pt.cex=1, legend=labels)
  box ()
}

################################################################################
########################## END function plot.ensemble ##########################
################################################################################
