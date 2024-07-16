#Appendix of Kleyer et al. XXXXXXXXXXXXX
########## The CLUS-MOD Procedure ##############
#Copyright of Barbara Strauss & Michael Kleyer 
#Landscape Ecology Group, University of Oldenburg, D-26111 Oldenburg, Germany
#e-mail: barbara.strauss@uni-oldenburg.de, michael.kleyer@uni-oldenburg.de
################################################
#SCRIPT 5 - SPECIES ASSIGNMENTS TO GROUPS
################

#This script produces a table where, for selected combinations,
#group assignements for each species and combination can
#be compared.

#A .txt file is produced that can be read in with Excel.
#Species marked "NA" in a certain combination were either 
#not in a cluster at all or not in a stable cluster.

#Note that for combinations that lead to almost identical
#groups, comparable group do not necessarily have the same number
#(even though they often do)!


#read in required objects from previous steps
#source("cluster.output.txt")	
#source("modelling.output.txt")
#source("output.table.txt")

spec.group.check<-matrix(nrow=nrow(traits),ncol=length(combis.to.compare))
rownames(spec.group.check)<-rownames(traits)            #Barbaras Original Script:    rownames(spec.group.check)<-traits[,1]
colnames(spec.group.check)<-combis.to.compare
no.of.groups<-vector(mode="numeric")

#loop over all combinations that are to be compared
for(i in 1:length(combis.to.compare)){

	no<-combis.to.compare[i]
	combi.info<-output.table[output.table[,"no.combi"]==no,]

	#how many groups does the combination have
	no.of.groups[i]<-nrow(combi.info)

	#assign group numbers to species names
	spec.group<-matrix(nrow=nrow(spec.plots.t),ncol=1)      
	row.names(spec.group)<-rownames(spec.plots.t)   

for(j in 1:nrow(combi.info)){
		temp<-strsplit(combi.info[j,"names.spec.modelled"],", ")[[1]]
		spec.group[temp,]<-j
	}

	#write species numbers for combi i
	spec.group.check[,i]<-spec.group[,1]
}

out<-rbind(combis.to.compare,no.of.groups,spec.group.check)
write.table(out,file=name.table,sep="\t",row.names=T,col.names=F)

#clear workspace
#rm(list = ls())

