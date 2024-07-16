
##  Function
partBio<-function(monoYield, obsYield, propSown){
  if(any(is.na(monoYield))){
    stop("NA not allowed in monoculture yield")
  }
  
  if(!all.equal(rownames(obsYield), rownames(propSown))){
    stop("All rownames of obsYield and propSown must be equal ")
  }
  
  propSown<-propSown/rowSums(propSown)
  presAbs<-ifelse(propSown>0,1,0)
  
  hector<-function(M,Yoi,P){
    M<-M*ifelse(P>0,1,NA)
    N<-sum(ifelse(P>0,1,0))
    Yo<-sum(Yoi)
    RYei<-P
    RYoi<-Yoi/M
    Yei<-RYei*M
    Ye<-sum(Yei,na.rm=T)
    deltaY<-Yo-Ye
    deltaRYi<-RYoi-RYei
    Mhat<-Ye
    Comp<-N*mean(as.numeric(deltaRYi),na.rm=T)*Mhat
    Sel<-N*cov(as.numeric(deltaRYi),as.numeric(M),use="pairwise.complete.obs")*(N-1)/N
    Trans<-Yo/max(M,na.rm=T)
    maxM<-max(M,na.rm=T)
    #spDom<-colnames(max(M,na.rm=T))
    result<-c(Comp,Sel,deltaY,Trans,maxM)#,spDom
    return(result)
  }
  obsYieldSpl<-split(obsYield,1:nrow(obsYield)) %>%
    lapply(., as.numeric)
  propSownSpl<-split(propSown,1:nrow(propSown)) %>%
    lapply(., as.numeric)
  myResult<-mapply(FUN = hector, list(monoYield), obsYieldSpl, propSownSpl)
  myResult<-t(myResult)
  colnames(myResult)<-c("Comp","Sel","NetEffect","Trans", "maxM") # ,"spDom"
  rownames(myResult) <- rownames(obsYield)
  return(myResult)
}