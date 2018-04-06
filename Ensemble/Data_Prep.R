Data_Prep <- function(Type,atr,outlier,transformation,resampling){
  # This function is to import the Data file 
  #Arguments: n is the number of sheet in excel file 
  # transformation : "log", "derv1", "snvt", "raw".
  
  if (Type=="leaf"){              # These are the position of starting of the rflectance data 
    n <- 1
    refCol <- 16
  } else if (Type=="shrub"){
    n <- 2
    refCol <- 16
  } else if (Type=="plot"){
    n <- 3
    refCol <- 14
  } else if (Type=="soil"){
    n<-4
  refCol<-2
  } else if (Type=="dryleaf_INL"){
    n <- 6
    refCol <- 8
  }else if (Type=="dryleaf_ACCP"){
    n <- 7
    refCol <- 16
  }else if (Type=="OPT"){
    n <- 8
    refCol <- 35
  } else if (Type=="lopex"){
    n <- 9
    refCol <- 29
  } else if (Type=="DryLeaf"){
    n <- 10
    refCol <- 15
  }
    
  
  Data <- read.xlsx("Data.xlsx", sheet = n, colNames = TRUE)
  nr <- nrow(Data)
  nc<- ncol(Data)
  y <- Data[[atr]]
  index<-complete.cases(y)  #To remove the NA valuse (plots with no N measurments)
  Reflectance <-as.matrix(Data[1:nr,refCol:nc])
  aa<-which(Reflectance < 0, arr.ind=TRUE)[,2]
  other <- Data[,1:refCol-1]   # other attributes
  
  if (length(aa>0)){
    Reflectance<- Reflectance[,-aa] #Remove negative columns  
  }
  
  

  if (n==3){
    Y <- as.numeric(y[index]) # other attributes
    other <- data.frame(other[index,])
    wav <- (as.numeric(colnames(Reflectance)) >= 400 & as.numeric(colnames(Reflectance)) <= 2500)
    Reflectance <- Reflectance[,wav]
    
    Spectra <- as.matrix(Reflectance[index,])    
    
    if (transformation=="derv1"){
      smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)      
      trans_Spectra <- savitzkyGolay(smoothing, p = 2, w = 11, m = 1)
      wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
      trans_Spectra<-trans_Spectra[,!wav]
      wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
      trans_Spectra<-trans_Spectra[,!wav]
    
          } else if (transformation == "smooth") {
            smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
            trans_Spectra <- smoothing
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
            trans_Spectra<-trans_Spectra[,!wav]
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
            trans_Spectra<-trans_Spectra[,!wav]
          } else if (transformation=="log"){
            smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
            trans_Spectra <- log(1/smoothing)
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
            trans_Spectra<-trans_Spectra[,!wav]
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
            trans_Spectra<-trans_Spectra[,!wav]
          } else if (transformation=="raw"){
            smoothing <- Spectra
            trans_Spectra <- smoothing
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
            trans_Spectra<-trans_Spectra[,!wav]
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
            trans_Spectra<-trans_Spectra[,!wav]
          } else if (transformation=="cr"){
            smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
            NewSpecLib<-speclib(smoothing,(as.numeric(colnames(smoothing))))
            ch_ratio <- transformSpeclib(NewSpecLib, method = "ch", out = "ratio")                    
            trans_Spectra <- spectra(ch_ratio)
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
            trans_Spectra<-trans_Spectra[,!wav]
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
            trans_Spectra<-trans_Spectra[,!wav]
            z <- which(colSums(trans_Spectra)==0, arr.ind=TRUE)    # Since the variance of zero vector is zero, I am trying to get rid of zero vectors. This will actually impacts the scaling
                if (length(z>0)){
                  trans_Spectra<- trans_Spectra[,-z] #Remove negative columns  
                  colnames(trans_Spectra) <- paste(as.numeric(colnames(smoothing[,-z])))
                } else {
                  colnames(trans_Spectra) <- paste(as.numeric(colnames(smoothing)))                    
                }
            
          } else if (transformation=="derv1log"){
           
            smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
            smoothing <-  log10(1/smoothing)
            trans_Spectra <- savitzkyGolay(smoothing, p = 2, w = 11, m = 1)
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
            trans_Spectra<-trans_Spectra[,!wav]
            wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
            trans_Spectra<-trans_Spectra[,!wav]
             
          }
          
    trans_Wavelength <- as.numeric(colnames(trans_Spectra))
    
          #write.xlsx(trans_Wavelength,file="Smoothed_Plot_Wave.xlsx")
          #Sys.setenv(R_ZIPCMD= "C:/Rtools/bin/zip")  # in case getting Rtools error when using write
          if (outlier == TRUE){
            x.out <- pcout(trans_Spectra, makeplot = TRUE)  # http://www.statistik.tuwien.ac.at/forschung/CS/CS-2006-3complete.pdf
            I <- x.out$wfinal01==1
            
            newList <- list("Spectra" = trans_Spectra[I,], "Wavebands" =trans_Wavelength, "atr"=Y[I],"other"=other )
          }
          else if(outlier == FALSE) {
            newList <- list("Spectra" = trans_Spectra, "Wavebands" =trans_Wavelength, "atr"=Y,"other"=other )
          }
  }
  else if (n==1 | n==2 | n==6 | n==7 | n==8 | n==9 | n==10 ){
        Y <- as.numeric(y[index])
        other <- data.frame(other[index,])  # other attributes
        wav <- (as.numeric(colnames(Reflectance)) >= 400 & as.numeric(colnames(Reflectance)) <= 2400)
        
            if (resampling==TRUE){
            Smoothed_Wavelength_plot <- read.xlsx("Smoothed_Plot_Wave.xlsx", sheet = 1, colNames = FALSE)
            Smoothed_Wavelength_plot <- as.numeric(Smoothed_Wavelength_plot[,])
            Reflectance<-as.matrix(Reflectance[index,])
            wav<-as.numeric(colnames(Reflectance))
            Spectra <- resample(Reflectance,wav,Smoothed_Wavelength_plot, "spline")
            Spectra<- Spectra[,-which(Spectra < 0, arr.ind=TRUE)[,2]] #Remove negative columns
            }
        else { 
          Spectra <- as.matrix(Reflectance[index,wav])
        }
        if (transformation=="derv1"){
          smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
          trans_Spectra <- savitzkyGolay(smoothing, p = 2, w = 11, m = 1)
                 if (n==2 ){
                   wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                   trans_Spectra<-trans_Spectra[,!wav]
                   wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                   trans_Spectra<-trans_Spectra[,!wav]
                }
        } else if (transformation == "smooth") {
          smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
          trans_Spectra <- smoothing
                 if (n==2){
                   wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                   trans_Spectra<-trans_Spectra[,!wav]
                   wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                   trans_Spectra<-trans_Spectra[,!wav]
                 }
          
        } else if (transformation=="log"){
          smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
          trans_Spectra <- log(1/smoothing)
                if (n==2){
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                  trans_Spectra<-trans_Spectra[,!wav]
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                  trans_Spectra<-trans_Spectra[,!wav]
                }
        } else if (transformation=="raw"){
          smoothing <- Spectra
          trans_Spectra <- smoothing
                if (n==2){
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                  trans_Spectra<-trans_Spectra[,!wav]
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                  trans_Spectra<-trans_Spectra[,!wav]
                }
        } else if (transformation=="cr"){    # continuum removal
          smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
          NewSpecLib<-speclib(smoothing,(as.numeric(colnames(smoothing))))
          ch_ratio <- transformSpeclib(NewSpecLib, method = "sh", out = "ratio")                    
          trans_Spectra <- spectra(ch_ratio)
                  if (n==2){
                    wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                    trans_Spectra<-trans_Spectra[,!wav]
                    wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                    trans_Spectra<-trans_Spectra[,!wav]
                  }
          z <- which(colSums(trans_Spectra)==0, arr.ind=TRUE)    # Since the variance of zero vector is zero, I am trying to get rid of zero vectors. This will actually impacts the scaling
            if (length(z>0)){
              trans_Spectra<- trans_Spectra[,-z] #Remove negative columns  
              colnames(trans_Spectra) <- paste(as.numeric(colnames(smoothing[,-z])))
            } else {
              colnames(trans_Spectra) <- paste(as.numeric(colnames(smoothing)))                    
            }
        } else if (transformation=="derv1log"){
          
          smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
          smoothing <-  log10(1/smoothing)
          trans_Spectra <- savitzkyGolay(smoothing, p = 2, w = 11, m = 1)
                if (n==2){
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
                  trans_Spectra<-trans_Spectra[,!wav]
                  wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
                  trans_Spectra<-trans_Spectra[,!wav]
                }
          
        }
          
        trans_Wavelength <- as.numeric(colnames(trans_Spectra))
        
          if (outlier == TRUE){
            x.out <- pcout(trans_Spectra, makeplot = TRUE)  # http://www.statistik.tuwien.ac.at/forschung/CS/CS-2006-3complete.pdf
            I <- x.out$wfinal01==1
            
            newList <- list("Spectra" = trans_Spectra[I,], "Wavebands" =trans_Wavelength, "atr"=Y[I],"other"=other )
          }
          else if(outlier == FALSE) {
            newList <- list("Spectra" = trans_Spectra, "Wavebands" =trans_Wavelength, "atr"=Y,"other"=other)
          }
        
  } else if (n==4){
    Y <- (y[index])
    wav <- (as.numeric(colnames(Reflectance)) >= 400 & as.numeric(colnames(Reflectance)) <= 2400)
    Spectra <- as.matrix(Reflectance[index,wav])
    smoothing <- savitzkyGolay(Spectra, p = 2, w = 11, m = 0)
    trans_Spectra <- smoothing
    wav <- (as.numeric(colnames(trans_Spectra)) >= 1300 & as.numeric(colnames(trans_Spectra)) <= 1450)
    trans_Spectra<-trans_Spectra[,!wav]
    wav <- (as.numeric(colnames(trans_Spectra)) >= 1750 & as.numeric(colnames(trans_Spectra)) <= 2000)
    trans_Spectra<-trans_Spectra[,!wav]
    trans_Wavelength <- as.numeric(colnames(trans_Spectra))
    newList <- list("Spectra" = trans_Spectra, "Wavebands" =trans_Wavelength, "atr"=Y,"other"=other)
  }
 
       
   
       return(newList)
  
}


