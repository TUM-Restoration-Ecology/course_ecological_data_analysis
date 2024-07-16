Inference_modelset<-function(Explanatory=NULL){

################################

   #Keep the number of explanatory variables
  
  NbVar<-ncol(Explanatory)
  combin<-c(colnames(Explanatory))
    
    i<-1
    while(i<=NbVar){
        if(!is.factor(Explanatory[,i]))
          combin[i]<-c(paste("s(", combin[i], ",2)", sep = ""))
        else combin[i]<-c(paste(combin[i]))
        i<-i+1
    }
  
   
  Mat<-as.data.frame(matrix(0, nrow=NbVar, ncol=NbVar))
  Mat[1,1:NbVar]<- t(combin)
  Temp<-combin

  Dm=seq(from=2, to=NbVar-1)
  a=1
  while(a<=length(Dm)){
    Maxi=dim(combn(Temp,Dm[a]))[2]
    i<-1
    while(i<=Maxi){
      combin = c(combin, paste(combn(Temp,Dm[a])[,i], collapse="+"))
      i<-i+1
    }
    Mat<-cbind(Mat, rbind(combn(Temp,Dm[a]), matrix(0, nrow=NbVar-nrow(combn(Temp,Dm[a])),
       ncol=ncol(combn(Temp,Dm[a])))))      
    a<-a+1
  }

  combin = c(combin, paste(Temp,collapse="+"))  
  Mat = cbind(Mat, Temp)

  
  return(list(combin, Mat))
  
}

