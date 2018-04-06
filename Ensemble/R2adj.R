R2adj<-function(n,m,R2){
  
  R2adj <- 1-((n-1)*(1-R2)/(n-m-1))
  return(R2adj)
  
}