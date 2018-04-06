Pred_R2 <- function(y,yhat){
  # y is the vector of observation (e.g. %N)
  # yhat is the vector of predicted values
  a<- y-yhat
  a<- a^2
  press <- sum(a)
  b<- y-mean(y)
  b <- b^2
  TSS<- sum(b)
  Pred_R2 <- 1-press/(TSS)
  return(Pred_R2)
}