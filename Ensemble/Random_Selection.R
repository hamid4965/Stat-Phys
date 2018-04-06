Random_Selection <- function(objc,val_ratio){
  #This script is just for randomly select a ratio of the data.
  # Y is the vector of response variable (e.g. nitrogen content)
  
  
  nr=length(objc$N)
  Ind <- 1:nr
  Ind_Selected <- sample(Ind, ceiling(val_ratio*nr),replace = FALSE) 
  Yval <- objc$N[Ind_Selected]                   # Select ratio% of data for validation
  Xval <- as.matrix(objc$spectra[Ind_Selected,])
  Ycal <- objc$N[setdiff(Ind,Ind_Selected)]      # Select % for calibration
  Xcal <- as.matrix(objc$spectra[setdiff(Ind,Ind_Selected),])
  newList <- list("Yval"=Yval,"Xval"=Xval,"Ycal"= Ycal,"Xcal"=Xcal)
  
  return(newList)
  
}

