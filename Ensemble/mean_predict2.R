mean_predict <- function(B, newdata,pls_scale){
  
  
  if (is.matrix(newdata)) {
    ## For matrices, simply check dimension:
    if (ncol(newdata) != dim(newdata)[2])
      stop("'newdata' does not have the correct number of columns")
    newX <- newdata
  }
  nobs <- dim(newX)[1]
  
  scale<- base::rowMeans(pls_scale,na.rm=T)
  newX <- newX / rep(scale, each = nobs)

  B <- B
  nr<-length(B)
  
  pred <- newX %*% B[-1] + rep(B[1], each = nobs)
  return(pred)
  }