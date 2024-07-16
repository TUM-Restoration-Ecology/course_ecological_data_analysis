##MEan DIssimilarity Components##
#samp:  community matrix; sites in lines, species in columns
#dis:   dissimilarity matrix
#type:  "both" for results with abundances weighted and non-weighted
#       "abundance" for results with abundances weighted
#       "presence" for results with abundances non-weighted
  
melodic <- function(samp,dis,type="both"){
	if(class(samp)!="matrix"){samp <- as.matrix(samp)}
  if(class(dis)!="matrix"){dis <- as.matrix(dis)}
  N<-dim(samp)[1]
	melodic<-list()
	if (type=="both"){
    melodic$abundance<-list()
    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
    melodic$presence<-list()
    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  }
	if (type=="abundance"){ 
    melodic$abundance<-list()
    melodic$abundance$mpd<-melodic$abundance$rao<-melodic$abundance$simpson<-numeric(N)
  }
	if (type=="presence"){ 
    melodic$presence<-list()
    melodic$presence$mpd<-melodic$presence$rao<-melodic$presence$simpson<-numeric(N)
  }
	for (i in 1:N){
    sppInSample<-names(samp[i,samp[i, ]>0])
	  melodic$richness[i]<-length(sppInSample)
		if (length(sppInSample)>1){
      sample.dis<-dis[sppInSample,sppInSample]
			abund.w<-numeric(length(sppInSample))
			if (type=="both" | type=="abundance"){
  			abund.w <- samp[i , sppInSample] / sum(samp[i , sppInSample])
        sample.weights <- outer(abund.w , abund.w)
        melodic$abundance$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
		    melodic$abundance$rao[i] <- sum(sample.weights * sample.dis)
		    melodic$abundance$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
			} 	
      if (type=="both" | type=="presence"){
        abund.nw <- rep(1 , length(sppInSample)) / length(sppInSample)
        sample.weights <- outer(abund.nw , abund.nw)
        melodic$presence$mpd[i] <- weighted.mean(sample.dis[lower.tri(sample.dis)] , sample.weights[lower.tri(sample.weights)])
        melodic$presence$rao[i] <- sum(sample.weights * sample.dis)
        melodic$presence$simpson[i] <- sum(2*sample.weights[lower.tri(sample.weights)])
			}	
		}	else {
		  if (type=="both" | type=="abundance"){
        melodic$abundance$mpd[i] <- melodic$abundance$rao[i] <- melodic$abundance$simpson[i] <-NA
      }
      if (type=="both" | type=="presence"){
        melodic$presence$mpd[i] <- melodic$presence$rao[i] <- melodic$presence$simpson[i] <-NA
		  }
    }
  }  	
	out<-melodic
	return(out)	
}

############
############ EXAMPLES
  library(picante)
  data(phylocom)
  distances<-cophenetic(phylocom$phylo)/max(cophenetic(phylocom$phylo))
  test<-mpd(phylocom$sample,distances, abundance.weighted =T) ### Rao values
  test.ab<-melodic(phylocom$sample,distances, type="abundance") ### mpd, rao and simpson values, abundance weighted
  test.pr<-melodic(phylocom$sample,distances, type="presence") ### mpd, rao and simpson values, not abundance weighted
  test.both<-melodic(phylocom$sample,distances, type="both") ### mpd, rao and simpson values, not abundance weighted

