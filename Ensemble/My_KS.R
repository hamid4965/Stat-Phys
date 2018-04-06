My_KS <- function (X){

  nr <- dim(X)[1]
  nc <- dim(X)[2]
  Rankk <- numeric(nr)
  out <- 1:nr
  D <- distli(X)
  Rankk[1:2] <- which(D==max(max(D)),arr.ind = TRUE)[1:2]
  out <- out[! out %in% Rankk[1:2]]
  
      iter <- 3
      while (iter<=nr){
        inn <- Rankk[which(Rankk>0)]
        Dsub <- D[inn,out]
        if (is.matrix(Dsub)) {
          minD <- apply(Dsub,2,min)
        } else if (is.vector(Dsub)) {
          minD <- min(Dsub)
        } 
        
          
        if (is.matrix(Dsub)) {
          indexmin <- apply(Dsub,2,min)
        } else if (is.vector(Dsub)) {
          indexmin <- which.min(Dsub)
        } 
        
       # indexmin <- apply(Dsub,2,which.min)
        maxD <- max(minD)
        indexmax <- which.max((minD))
        Vadd <- out[indexmax]
        Rankk[iter] <- Vadd
        out <- out[! out %in% Rankk[1:iter]]
        iter <- iter+1
      }
  
      return(Rankk)
}

distli <- function(X){
  X <- t(X)
  D <- dim(X)[1]
  N <- dim(X)[2]
  X2 <- apply(X^2,2,sum)
  D <- repmat(rbind(X2),N,1) + repmat(cbind(X2),1,N)-2*t(X)%*%X
}

repmat <- function(X,m,n)
  {
    # % REPMAT R equivalent of repmat (matlab)
    # % FORMAT
    # % DESC 
    # % description not available.
    
    if (is.matrix(X))
    {
      mx <- dim(X)[1]
      nx <- dim(X)[2]
      out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
    }
    else if (is.vector(X)) {
      mx <- 1
      nx <- length(X)
      out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
    }
    else if (length(X) == 1)
    {
      out <- matrix(X,m, n)
    }
    return (out)
  }

# test the code (I took the X from the NIR in spectra dataset matlab)
#X <- read.xlsx("remove.xlsx", sheet = 2, colNames = FALSE)
