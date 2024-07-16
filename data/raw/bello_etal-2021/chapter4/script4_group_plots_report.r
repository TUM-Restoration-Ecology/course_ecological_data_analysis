#Appendix of Kleyer et al. XXXXXXXXXXXXX
########## The CLUS-MOD Procedure ##############
#Copyright of Barbara Strauss & Michael Kleyer 
#Landscape Ecology Group, University of Oldenburg, D-26111 Oldenburg, Germany
#e-mail: barbara.strauss@uni-oldenburg.de, michael.kleyer@uni-oldenburg.de

#read in required objects from previous steps
source("cluster.output.txt")	
source("modelling.output.txt")
source("output.table.txt")

#open empty pdf file to write output

par(mfrow=c(4,2)) # Subsequent figures will be drawn in an nr-by-nc array on the device mfrow. Delete if you want to have single graphs                  

#loop over all combinations that are to be plotted
for(i in 1:length(combis.to.plot)){

	no<-combis.to.plot[i]
	combi.info<-output.table[output.table[,"no.combi"]==no,]
	traits.involved<-strsplit(combi.info[1,"trait.combi"],", ")[[1]]

	#start new page for each combi
	for(a in 1:8){
		plot.new()
	}
	par(mfg=c(1,1))             #i and j indicate which figure in an array of figures is to be drawn next.  Delete if you want to have single graphs
	#plot combi info
	temp<-as.matrix(strsplit(combi.info[1,"trait.combi"],", ")[[1]])
	colnames(temp)<-"involved traits"
	textplot(temp,show.rownames=F,show.colnames=T,halign="left",
		hadj=0,cex=1,valign="top")
	title(paste("Combi No.",no),cex=2,adj=0)

	#assign group numbers to species names
	spec.group<-matrix(nrow=nrow(spec.plots.t),ncol=1)
  row.names(spec.group)<-row.names(spec.plots.t)      


	for(j in 1:nrow(combi.info)){
    temp<-strsplit(combi.info[j,"names.spec.modelled"],", ")[[1]]
		spec.group[temp,]<-j
	}

	#boxplot of trait values for each group and each trait
	for (j in 1:length(traits.involved)){
		boxplot.n(traits[,traits.involved[j]]~spec.group,top=T)
		title(traits.involved[j])
	}

	#plot group models for each environmental variable
	###################################################

	#plot legend with group colours and R2 of groups
	group.names<-vector(mode="character")
	for(k in 1:nrow(combi.info)){
		group.names[k]<-paste("Group",k)
	}

	temp=matrix(ncol=2,nrow=length(group.names))
	colnames(temp)<-c("Group","R^2")
	temp[,1]<-group.names
	temp[,2]<-round(as.numeric(combi.info[,"R2"]),digits=2)

	textplot(temp,cex=1,halign="left",show.rownames=F,hadj=0, cmar=5)
	legend(x=0,y=1,legend=group.names,
		col=mycol[1:length(group.names)],lty=1,lwd=5)

	
	#matrix for medians of environmental variables that 
	#are to be held constant (i.e. not the variable to be plotted)
	dummy.indep<-matrix(ncol=2*ncol(env),nrow=100)               
	colnames(dummy.indep)<-c(colnames(env),                
			paste("I(",colnames(env),"^2)",sep=""))      
	temp<-c(apply(env,2,median))	             
	temp<-as.vector(c(temp,temp^2))
	for(r in 1:nrow(dummy.indep)){
		dummy.indep[r,]<-temp
	}

	#loop over all environmental variables j
	for(j in 1:ncol(env)){                   

		#values for the env. var. to be plotted range
		#from their min to max values in the data
		env.var<-colnames(env)[j]               
		x.plot<-seq(from=min(env[,j]),to=max(env[,j]),         
				length.out=100) 

		temp.dummy.indep<-dummy.indep
		temp.dummy.indep[,j]<-x.plot
		temp.dummy.indep[,ncol(env)+j]<-x.plot^2      

		#matrix for y-values for all groups
		y.plot<-matrix(ncol=nrow(combi.info),nrow=100)	

		#loop over all groups k
		for(k in 1:nrow(combi.info)){	
			#extract model coefficients	
			betas<-combi.info[k,paste("beta",names.indep.var, 
				sep="_")[-1]]
			betas[is.na(betas)]<-0
			temp<-names(betas)
			betas<-as.numeric(betas)
			names(betas)<-temp

			#calculate products of coefficients and env. vars
			products<-temp.dummy.indep
			colnames(products)<-paste("beta",colnames(temp.dummy.indep),sep="_")
			for(b in 1:length(betas)){
				products[,names(betas)[b]]<-products[,names(betas)[b]]*betas[b]
			}
			product.sums<-rowSums(products)+
				as.numeric(combi.info[k,"beta_Intercept"])

			#calculate predicted values
			y.plot[,k]<- 1/(1+exp(-1*product.sums))
		}

		#plot predicions for each group
		plot(1, type="n", 
			ylim=c(0,1),xlim=c(x.plot[1],x.plot[100]),
			xlab=colnames(env[j]),                       
			ylab="Frequency / Coverage of group",
			main=colnames(env[j]))                 
		for(k in 1:ncol(y.plot)){
			if(combi.info[k,paste("weight",env.var, sep="_")]>min.weight&&
			 is.na(combi.info[k,paste("weight",env.var, sep="_")])==F){
				points(x.plot,y.plot[,k],type="l",col=mycol[k],lwd=3)
			}
		}

		#if plot.points==1, plot data points of original data
		if(plot.points==1){
			for(k in 1:ncol(y.plot)){
			   if(y.plot[1,k]!=y.plot[100,k]){
				temp<-rownames(spec.group)[spec.group==k]
				temp<-is.element(spe[,1],temp)             
				spec.group.plot<-colSums(spe[temp==T,][,-1])   
				points(env[,j], spec.group.plot/100,col=mycol[k])  
			   }
			}
		}	
	}	


}




