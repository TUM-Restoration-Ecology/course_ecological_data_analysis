#Correlation ratio: a measure of the relationship between the statistical dispersion 
#within individual categories and the dispersion across the whole population or sample.

#The correlation ratio takes values between 0 and 1. 
#The limit = 0 represents the special case of no dispersion among the means of the 
#different categories, while  = 1 refers to no dispersion within the respective categories

#Wilfried  THUILLER
#LECA




cor.ratio <- function(X,Xfac, weights){
  nr <- nrow(X)
  if(nrow(Xfac) != nr)
    stop("non convenient dimension")
  if(length(weights) != nr)
    stop("non convenient dimension")
  
  weights <- weights/sum(weights)
  rcor <- matrix(0, ncol(Xfac), ncol(X))
  floc <- function(x, fac, weights) {
    xwt <- x * weights
    poicla <- unlist(tapply(weights, fac, sum))
    moytot <- sum(xwt)
    z <- unlist(tapply(xwt, fac, sum))/poicla
    return(sum(poicla * (z-moytot)^2)/sum(weights * (x-moytot)^2))
  }
  for(i in 1:ncol(Xfac)){
    for(j in 1:ncol(X)){
      rcor[i,j] <- floc(X[,j],Xfac[,i], weights)
    }
  }
  
  rcor <- data.frame(rcor)
  row.names(rcor) <- names(Xfac)
  names(rcor) <- names(X)
  return(sqrt(rcor))
}
