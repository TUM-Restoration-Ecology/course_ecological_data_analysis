
#Appendix of Kleyer et al. XXXXXXXXXXXXX
########## The CLUS-MOD Procedure ##############
#Copyright of Barbara Strauss & Michael Kleyer 
#Landscape Ecology Group, University of Oldenburg, D-26111 Oldenburg, Germany
#e-mail: barbara.strauss@uni-oldenburg.de, michael.kleyer@uni-oldenburg.de
################################################
#SCRIPT 1 - CLUSTERING
#############################################

#Version: June 6, 2009

#- This script reads in a matrix containing species and their trait values. 

#- Cluster analyses are performed to group the 
#  species according to their trait values. Clustering is done
#  based on different numbers and combinations of traits.

#- The optimal number of groups per cluster is calculated.

#- The resulting clusterings and their clusters
#  are bootstrapped to access cluster stability.


#Required INPUT TABLE (tab-delimeted .txt-file, stored in the search path): 
#trait table: traits (columns) x species (rows) matrix
#1st column: Species abbreviation ; no empty spaces are allowed!  Use "." or "_" instead
#			No double names!
#2nd column onward: traits 1 to n

#######################################################################
#Required LIBRARIES that need to be INSTALLED:
##########################################################################
#library(MASS),library(cluster),library(vegan),library(fpc),library(clusterSim)

#Readable OUTPUT is written to .txt file (set name with
#name.output.table). More output is produced to file
#"cluster.output.txt". Do not delete this file, it is
#required for the following steps.

################################
####_____Start___Script__________
################################

#if max.traits > no. of traits available in data, set to no. of traits
if(ncol(traits)-1 > max.traits) max.traits<-ncol(traits)-1

#select traits as set in traits.to.consider
if(mode(traits.to.consider)=="character"){
	if(traits.to.consider[1]=="all") {
		traits.to.consider<-c(2:ncol(traits))
		traits<-cbind(traits[,1,drop=F],traits[,traits.to.consider,drop=F])	
	} else traits<-cbind(traits[,1,drop=F],traits[,traits.to.consider,drop=F])
}
if(mode(traits.to.consider)=="numeric") {
	traits<-cbind(traits[,1,drop=F],traits[,traits.to.consider,drop=F])	
}

#Check, if all trait information is available for all species
#If not, script stops and an error message is printed.
stop.script<-0
if(sum(is.na(traits))>0){
	missing.traits<-traits
	missing.traits[,2:ncol(missing.traits)]<-NA
	temp<-vector("numeric",length=nrow(traits))
	for(i in 1:nrow(traits)){
		if(sum(is.na(traits[i,]))>0) {
			missing.traits[i,][is.na(traits[i,])]<-1
			temp[i]<-1
		}
	}
	missing.traits<-missing.traits[temp==1,]
	missing.traits<-cbind(rownames(missing.traits),missing.traits)
	write.table(missing.traits,file="missing.traits.txt",sep="\t",
		row.names=F,col.names=T)
	stop.script<-1
}

#continue only if all trait information is available for all species
if(stop.script!=1) {

####################################################################
#All Clusters with different combinations of traits
####################################################################

subsets<-function(n, k, set){
	if(k <= 0) NULL
	else if(k >= n) set
	else rbind(cbind(set[1], Recall(n-1, k-1, set[-1])), Recall(n-1, k, set[-1]))
}

n.traits<-ncol(traits)      
names.traits<-colnames(traits)[1:ncol(traits)]  
trait.combinations<-vector("list")
for(i in 1:max.traits){	
	temp<-subsets(n.traits,max.traits-i+1,names.traits)
	if(is.vector(temp)==T) temp<-t(as.matrix(temp))
	temp<-apply(temp,1,as.list)
	trait.combinations<-c(trait.combinations,temp)		
}


####################################################################
#Calculate the cophenetic correlation for all clusters
#Cophenetic correlation is the correlation coefficient between the values
#of the cophenetic matrix (the matrix describing the distances between points in
#the dendrogramm) and the distance matrix. It measures the extent to which the 
#clustering result corresponds to the original resemblance matrix.
#Spearman rank correaltion is used (Legendre 2006,p.375f).
#If the correlation is low, clusters are usually not stable when
#they are bootstrapped.
####################################################################
cophenetic.corr<-vector(mode="numeric")
traits.standard<-decostand(traits,method="standardize") 
#print("Calculate cophenetic correlation for")
for(i in 1:length(trait.combinations)){
	#print(paste("...Clustering ",i,sep=""))
	matrix.to.cluster<-traits.standard[,unlist(trait.combinations[i])]
	matrix.to.cluster.dist<-dist(matrix.to.cluster,method="euclidean")
	cluster<-hclust(matrix.to.cluster.dist,method="ward")
	cophenetic.corr[i]<-cor(matrix.to.cluster.dist,cophenetic(cluster),method="spearman")	
}

#Delete all clusters that have a cophenetic correlation 
#below the defined threshold.
combi.tracker<-1:length(trait.combinations)	#track the no. of combis retained 
good.combinations<-trait.combinations[cophenetic.corr>=min.cophenetic.corr]
combi.tracker<-combi.tracker[cophenetic.corr>=min.cophenetic.corr]

####################################################################
#Calculate which number of clusters is appropriate.
#Maximum number of clusters which is tested is set with max.no.groups
####################################################################

#result matrix for G1 (Gordon 1999, p.62): Calinski and Harabasz's index
g1<-matrix(nrow=length(good.combinations),ncol=max.no.groups-1)
colnames(g1)<-c(2:max.no.groups)

#empty result matrix, if clusters with numbers of cases <= the noisecut.value
#were excluded as "noise"
noise<-matrix(nrow=length(good.combinations),ncol=max.no.groups-1,0)
colnames(noise)<-c(2:max.no.groups)

for(clust in 1:length(good.combinations)){
	#print(paste("Find optimal number of groups for clustering ",clust,sep=""))
	matrix.to.cluster<-traits.standard[,unlist(good.combinations[clust])]
	matrix.to.cluster<-as.matrix(matrix.to.cluster)
	rownames(matrix.to.cluster)<-row.names(traits.standard)
	for(i in 2:max.no.groups){
		matrix.to.cluster.dist<-dist(matrix.to.cluster,method="euclidean")
		cluster<-disthclustCBI(matrix.to.cluster.dist,method="ward",
			cut="number",k=i,noisecut=noisecut.value)
		if(cluster$noise==F){
			g1[clust,i-1]<-index.G1(as.matrix(matrix.to.cluster.dist),cluster$partition)
		}
		if(cluster$noise==T) {
			if(cluster$nc<i) break
			noise[clust,i-1]<-1
			if(i>2){
				spec.exclude<-unlist(cluster$clusterlist[cluster$nc])
				matrix.to.cluster2<-matrix.to.cluster[spec.exclude==F,]
				matrix.to.cluster.dist2<-dist(matrix.to.cluster2,method="euclidean")
				cluster<-disthclustCBI(matrix.to.cluster.dist2,method="ward",
					cut="number",k=i-1)
				#if(cluster$noise==T) noise[clust,i-1]<-noise[clust,i-1]+1
				#temp.g1<-index.G1(as.matrix(matrix.to.cluster.dist2),cluster$partition)	
				#if(g1[clust,i-1]<temp.g1) g1[clust,i-1]<-temp.g1
				
				g1[clust,i-1]<-index.G1(as.matrix(matrix.to.cluster.dist2),cluster$partition)	
			}
			
		}
	}
}

#round the g1-result-matrix
g1.round<-round(g1,digits=0)


#Extract which is the best number of groups
#If g1 is the same for higher numbers of groups, the smallest is chosen.
#(for the sake of parsimony...)
best.no<-vector(mode="numeric")
g1.round[is.na(g1.round)]<-0
for(i in 1:nrow(g1.round)){
	temp<-colnames(g1.round)[g1.round[i,]==max(g1.round[i,])]
	best.no[i]<-min(as.numeric(temp))
}
best.no<-as.numeric(best.no)

#result-vector, if the best number of clusters contains "noise" cases
noise.vector<-vector(mode="numeric")
for(i in 1:length(best.no)){
	no<-as.character(best.no[i])
	if(as.vector(noise[i,no])==1) noise.vector[i]<-1 else noise.vector[i]<-0	
}

##########################################################################
#Bootstrap each clustering with the "best" numbers of clusters to
#see which clusters are stable and how many species are in stable clusters
##########################################################################

#result matrix for the bootstrapping procedure
result.cluster.boot<-matrix(ncol=9,nrow=length(good.combinations))
colnames(result.cluster.boot)<-c("clust.stable","clust.instable","spec.total","spec.stable.clust",
	"spec.instable.clust","spec.not.clust","stable.clust","spec.per.clust","mean.jacc.clust.stable")
rownames(result.cluster.boot)<-combi.tracker
result.cluster.boot<-as.data.frame(result.cluster.boot)
result.cluster.boot$spec.total<-nrow(traits.standard)

trait.group.values.all<-vector("list")

spec.exclude.list<-vector(length=length(good.combinations))
for(i in 1:length(good.combinations)){
	#print(paste("Bootstrapping of Clustering ",i, " (out of ",length(good.combinations),")",sep=""))
	matrix.to.cluster<-as.matrix(traits.standard[,unlist(good.combinations[i])])
	rownames(matrix.to.cluster)<-rownames(traits.standard)
	matrix.to.cluster.dist<-dist(matrix.to.cluster,method="euclidean")
 
	if(noise.vector[i]==0){
		cluster<-disthclustCBI(matrix.to.cluster.dist,method="ward",cut="number",k=best.no[i]) 
		temp<-clusterboot(matrix.to.cluster.dist,B=no.boot,distances=T,bscompare=T,
			bootmethod="boot",clustermethod=disthclustCBI,method="ward",
			cut="number",k=best.no[i],count=F)  
		result.cluster.boot$spec.not.clust[i]<-0
	}

	if(noise.vector[i]==1){
		cluster<-disthclustCBI(matrix.to.cluster.dist,method="ward",cut="number",k=best.no[i],
			noisecut=noisecut.value)
		spec.exclude<-unlist(cluster$clusterlist[cluster$nc])
		spec.exclude.list[i]<-paste(rownames(matrix.to.cluster)[spec.exclude==T],collapse=",")
		matrix.to.cluster2<-matrix.to.cluster[spec.exclude==F,]
		matrix.to.cluster.dist2<-dist(matrix.to.cluster2,method="euclidean")
		cluster<-disthclustCBI(matrix.to.cluster.dist2,method="ward",cut="number",k=best.no[i]-1,
			noisecut=noisecut.value)
		temp<-clusterboot(matrix.to.cluster.dist2,B=no.boot,distances=T,bscompare=T,
			bootmethod="boot",clustermethod=disthclustCBI,method="ward",
			cut="number",k=cluster$nccl,count=F)#$bootmean
		result.cluster.boot$spec.not.clust[i]<-sum(spec.exclude)
	}

	temp2<-vector(mode="numeric")
	clustering.traits<-traits[,unlist(good.combinations[i])]
	clustering.traits<-as.matrix(clustering.traits)
	colnames(clustering.traits)<-unlist(good.combinations[i])
	rownames(clustering.traits)<-rownames(traits)
	trait.group.values<-matrix(ncol=length(unlist(good.combinations[i]))+2,nrow=cluster$nc)
	colnames(trait.group.values)<-c("bootmean","no.species",unlist(good.combinations[i]))
	trait.group.values[,1]<-round(temp$bootmean,digits=2)
	for(g in 1:cluster$nc){
		temp2[g]<-sum(cluster$partition==g)
		trait.group.values[g,2]<-temp2[g]
		group.traits<-clustering.traits[names(cluster$partition[cluster$partition==g]),]
		group.traits<-as.matrix(group.traits)
		for(t in 1:ncol(group.traits)){
			quantiles<-round(quantile(group.traits[,t],probs=c(0.5,0.25,0.75)),digits=2)
			trait.group.values[g,t+2]<-paste(quantiles[1],
				"[",quantiles[2],"|",quantiles[3],"]",collapse="")
			
		}
	}
	result.cluster.boot$spec.per.clust[i]<-paste(temp2,collapse=",")
	
	clust.stable<-cbind(seq(1:temp$nc),temp$bootmean>=mean.jac,temp$bootmean)
	jacc.clust.stable<-clust.stable[,3][clust.stable[,2]==1]
	clust.stable<-clust.stable[,1][clust.stable[,2]==1]

	result.cluster.boot$stable.clust[i]<-paste(clust.stable,collapse=",")
	result.cluster.boot$mean.jacc.clust.stable[i]<-paste(round(jacc.clust.stable,digits=2),collapse=",")

	result.cluster.boot$clust.stable[i]<-length(clust.stable)
	result.cluster.boot$clust.instable[i]<-temp$nc-length(clust.stable)

	result.cluster.boot$spec.stable.clust[i]<-sum(is.element(temp$partition,clust.stable))
	result.cluster.boot$spec.instable.clust[i]<-
			length(temp$partition)-result.cluster.boot$spec.stable.clust[i]

	trait.group.values.all[i]<-list(trait.group.values)

}


result.cluster.boot<-cbind(result.cluster.boot,best.no,noise.vector)

#Proportion of species in stable clusters
prop.spec<-result.cluster.boot$spec.stable.clust/result.cluster.boot$spec.total
stable<-result.cluster.boot$stable.clust>1

no.stable<-result.cluster.boot$clust.stable

final.combinations<-good.combinations[prop.spec>=min.prop&stable&no.stable>1]
final.result.cluster.boot<-result.cluster.boot[prop.spec>=min.prop&stable&no.stable>1,]
final.trait.group.values<-trait.group.values.all[prop.spec>=min.prop&stable&no.stable>1]
combi.tracker<-combi.tracker[prop.spec>=min.prop&stable&no.stable>1]

#######################
#some more output
#write which traits are involved in respective trait combination
involved.traits<-matrix(ncol=3,nrow=length(final.combinations))
colnames(involved.traits)<-c("no.combi","no.traits","involved.traits")
for(i in 1:length(final.combinations)){
	involved.traits[i,"involved.traits"]<-paste(unlist(final.combinations[[i]]),collapse=", ")
	involved.traits[i,"no.traits"]<-length(final.combinations[[i]])
	involved.traits[i,"no.combi"]<-rownames(final.result.cluster.boot)[i]
}

#write checklist for traits in combination
trait.checklist<-matrix(ncol=length(names.traits),nrow=length(final.combinations))
colnames(trait.checklist)<-names.traits

for(i in 1:length(final.combinations)){
	temp<-unlist(final.combinations[i])
	for(j in 1:length(temp)){
		trait.checklist[i,temp[j]]<-1
	}
}

#make output table
result.cluster.boot<-cbind(involved.traits,cophenetic.corr[combi.tracker],
	final.result.cluster.boot,trait.checklist)
colnames(result.cluster.boot)[colnames(result.cluster.boot)==
		"cophenetic.corr[combi.tracker]"]<-"cophenetic.corr"
result.cluster.boot[,"cophenetic.corr"]<-round(as.numeric(result.cluster.boot[,"cophenetic.corr"]),digits=2)

result.cluster.boot<-result.cluster.boot[,colnames(result.cluster.boot)!="noise.vector"]
result.cluster.boot[,"best.no"]<-result.cluster.boot[,"clust.stable"]+
			result.cluster.boot[,"clust.instable"]
#colnames(result.cluster.boot)[colnames(result.cluster.boot)=="best.no"]<-"no.clusters"
result.cluster.boot<-cbind(result.cluster.boot[,1:4],result.cluster.boot[,14], result.cluster.boot[,5:13],result.cluster.boot[,15:ncol(result.cluster.boot)])
colnames(result.cluster.boot)[5]<-"no.clusters"


write.table(result.cluster.boot,file=name.output.table,sep="\t",row.names=F)


#Dump important objects in file.
#This file is needed for partII of the skripts, where models are calculated.
dump(list=c("traits","traits.standard","final.result.cluster.boot",	
	"final.combinations","final.trait.group.values","result.cluster.boot",
	"noisecut.value"),file="cluster.output.txt")

}

if(stop.script==1) warning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		There are species with missing trait values. 
		First complete the traits file and then restart the script.
		A table with an overview of the missing values has been 
		written to the file 'missing.traits.txt' 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

