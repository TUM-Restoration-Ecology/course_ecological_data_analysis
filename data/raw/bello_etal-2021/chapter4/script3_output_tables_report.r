#Appendix of Kleyer et al. XXXXXXXXXXXXX
########## The CLUS-MOD Procedure ##############
#Copyright of Barbara Strauss & Michael Kleyer 
#Landscape Ecology Group, University of Oldenburg, D-26111 Oldenburg, Germany
#e-mail: barbara.strauss@uni-oldenburg.de, michael.kleyer@uni-oldenburg.de
########################################################
#SCRIPT 3 - CREATE OUTPUT TABLES WITH MODELLING RESULTS
########################################################

#read in required objects
source("cluster.output.txt")	
source("modelling.output.txt")

#extract names of traits used for clustering
names.traits<-names(traits)[1:length(names(traits))]   

#extract names of independent variables present in the models
names.indep.var<-vector(mode="character")
for(i in 1:length(mod.clust$beta.mean)){
	for(j in 1:length(mod.clust$beta.mean[[i]])){
		names.indep.var<-c(names.indep.var,rownames(mod.clust$beta.mean[[i]][[j]]))	
	}
}
names.indep.var<-unique(names.indep.var)

#create empty matrix for output
output.table<-matrix(ncol=6+length(names.traits)+length(names.indep.var)+length(names.indep.var))
colnames(output.table)<-c("no.combi","trait.combi","group.no","R2","no.spec.modelled","names.spec.modelled",
	names.traits,paste("beta",names.indep.var, sep="_"),paste("weight",names.indep.var, sep="_"))


#write results in output table for all combinations i
for(i in 1:length(final.combinations)){
	no.groups<-length(unlist(mod.clust$r2.groups[i]))
	no.stable.groups<-as.numeric(unlist(strsplit(final.result.cluster.boot$stable.clust[i],split=",")))
	#temporary output matrix for each group
	out.temp<-matrix(ncol=(ncol(output.table)-length(names.indep.var)),nrow=no.groups)
	colnames(out.temp)<-c("no.combi","trait.combi","group.no","R2","no.spec.modelled","names.spec.modelled",
	names.traits,names.indep.var)

	#temporary output matrix for variable weights of each group
	out.temp.2<-matrix(ncol=length(names.indep.var),nrow=no.groups)
	colnames(out.temp.2)<-names.indep.var

	#write combination number
	out.temp[,"no.combi"]<- rownames(final.result.cluster.boot)[i]

	#write involved traits
	out.temp[,"trait.combi"]<-  paste(unlist(final.combinations[i]),collapse=", ")

	#write R2 of averaged model for each group
	out.temp[,"R2"]<- unlist(mod.clust$r2.groups[i])

	#write number of species for each group
	out.temp[,"no.spec.modelled"]<- unlist(mod.clust$no.modelled.group[i])

	#write names of modelled species and group number for each group
	for(j in 1:no.groups){
		out.temp[j,"names.spec.modelled"]<-	unlist(mod.clust$spec.modelled[[i]][j])	
		out.temp[j,"group.no"]<-j
	}

	#write median and interquartile range for values of each trait t for each group
	for(t in (1:(ncol(final.trait.group.values[[i]])-2))){
		temp<-final.trait.group.values[[i]][no.stable.groups,]
		temp<-temp[,-c(1,2),drop=F]
		out.temp[,colnames(temp)[t]]<-temp[,t]
	}

	#write coefficients of averaged models
	for(j in 1:length(mod.clust$beta.mean[[i]])){
		if(is.null(mod.clust$beta.mean[[i]][[j]])==T) next
		temp<-t(mod.clust$beta.mean[[i]][[j]])
		out.temp[j,colnames(temp)]<-temp[1,]
	  }
	
	#write variable weights of averaged models
	  for(j in 1:length(mod.clust$beta.mean[[i]])){
		if(is.null(mod.clust$beta.mean[[i]][[j]])==T) next
		if(is.na(mod.clust$beta.mean[[i]][[j]])==T) next
    temp<-t(mod.clust$beta.mean[[i]][[j]])
		out.temp.2[j,colnames(temp)]<-temp[2,]
	}
	
	#combine variable weights with rest of the output
	out.temp.all<-cbind(out.temp,out.temp.2)

	#combine output for variable combination i with those of the previous combinations
	output.table<-rbind(output.table,out.temp.all)
} #end of combinations loop

#delete first line of the output (line is empty)
output.table<-output.table[-1,]

#put species names in last column (reads more easily)
output.table<-cbind(output.table[,-6],output.table[,6])
colnames(output.table)[ncol(output.table)]<-"names.spec.modelled"


#write output.table.1 to file
output.table[,"R2"]<-round(as.numeric(output.table[,"R2"]),digits=2)
write.table(output.table,file=name.outputfile.1,sep="\t",row.names=F)


#write output.table.2 to file
r2.all<-vector(mode="character")
for(i in 1:length(final.combinations)){
	r2.all<-c(r2.all,paste(round(unlist(mod.clust$r2.groups[[i]]),digits=2),collapse=", "))
}

r2.output<-cbind(round(mod.clust$r2.av,digits=2),round(mod.clust$min.r2,digits=2),
	r2.all)
colnames(r2.output)<-c("r2.av","r2.min","r2.all")
output.table.2<-cbind(result.cluster.boot[,1:4],r2.output,
	result.cluster.boot[,5:ncol(result.cluster.boot)])
write.table(output.table.2,file=name.outputfile.2,sep="\t",row.names=F)

#dump important object into file
dump(list=c("names.indep.var","output.table"),file="output.table.txt")

#clear workspace
#rm(list = ls())



