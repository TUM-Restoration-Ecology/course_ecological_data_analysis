#Appendix of Kleyer et al. XXXXXXXXXXXXX
########## The CLUS-MOD Procedure ##############
#Copyright of Barbara Strauss & Michael Kleyer 
#Landscape Ecology Group, University of Oldenburg, D-26111 Oldenburg, Germany
#e-mail: barbara.strauss@uni-oldenburg.de, michael.kleyer@uni-oldenburg.de
################################################
#SCRIPT 2 - MODELLING
################################################

#- Calculates models for functional groups (group frequencies 
#  or coverages depending on environmental parameters). 

#- First, univariate models are estimated to find significant variables
#  and determine form of relationships (monotonous or unimodal). Second,
#  multiple models are estimated. Third, model averaging is applied
#  to all significant univariate and multiple models.

#- Output from part1 of the Script is required (read in from file "cluster.output.txt")

#- It is possible to use the results of a clustering that considered
#  additional species that are not present in the plot data.

#- To create readable OUTPUT, run script 3! Unreadable ouptut is written
#  to file "modelling.output.txt". Do not delete this file, it is 
#  required for later steps.

#Required INPUT TABLES (tab-delimeted .txt-files, stored in the search path)
#1. A species x plots table 
#   Species (rows) x plots (columns) matrix with species frequencies
#   or coverage [%, 0-100].
#   1st column: species number/abbreviation, may be the same
#		as 2nd column, no double names!
#   2nc column: species full names, no spaces allowed!
#   1st row: plot names. No spaces!
#   The order of the species

#2. A plots x environment table 
#   Plots (rows) x environment (columns) matrix containing environmental parameters.
#   1st column: plot names. No spaces!
#   1st row: names of envrionmental variables. No spaces!



################################################################
#Requried LIBRARIES that need to be installed:
################################################################
#library(fpc), library(cluster), library(vegan), library(Hmisc), 
#library(lmtest)

###############################################################

#read in results from script 1
source("cluster.output.txt")			

#species occurences & species names  (comment this line out if species are already in rows and plot in colummns)
spec.plots.t<-t(spe)    # Transpose to species in rows and plots in columns      


spec.full.names<-matrix(rownames(spec.plots.t))
rownames(spec.full.names)<-rownames(spec.plots.t) #species numbers / abbreviations   
#exclude column with species names and scale covers to [0,1]
spec.plot.covers<-spec.plots.t/100	     
spec.plot<-rownames(spec.plot.covers) 	#species found on the plot

#################################################################
#Check, if trait information is available for all the species in 
#the plots x species matrix. If not, script stops and 
#an error message is printed.
##################################################################
stop<-0
match<-is.element(row.names(spec.plots.t),row.names(traits))
if(sum(as.numeric(match))!=nrow(spec.plots.t)) {
	missing<-spec.plots.t[match=="FALSE",]
	missing<-cbind(row.names(missing),missing[,1])
	write.table(missing,file="missing.species.txt",sep="\t",
		col.names=c("SBS.nr.layer","SBS.name"),row.names=F)	#hier muss ein neuer Eintrag geschehen
	stop<-1	
}

if(stop==0){										#0o (start modelling if all 
												#trait info is available)

#############################################################
#1. Recreate the clustering for final.combinations[i]
#    and store which species belong to which cluster-group
#############################################################
#list for output
mod.clust<-vector("list",8)		 
names(mod.clust)<-c("final.combi","r2.groups","r2.av","min.r2","beta.mean",
			"no.modelled.group","no.modelled.clustering","spec.modelled") 

#start loop over all clusterings i
#1o (i, clustering loop)
for(i in 1:length(final.combinations)){ 					
	#print(paste("---------- CLUSTERING ",i," ------------",sep=""))

	matrix.to.cluster<-as.matrix(traits.standard[,unlist(final.combinations[i])])
	rownames(matrix.to.cluster)<-row.names(traits.standard)
	matrix.to.cluster.dist<-dist(matrix.to.cluster,method="euclidean")
	no.stable.groups<-as.numeric(unlist(strsplit(final.result.cluster.boot$stable.clust[i],split=",")))

	cluster<-disthclustCBI(matrix.to.cluster.dist,method="ward",cut="number",
		k=final.result.cluster.boot$clust.stable[i]+final.result.cluster.boot$clust.instable[i]+
		final.result.cluster.boot$noise.vector[i],noisecut=noisecut.value)
	

	if(final.result.cluster.boot$noise.vector[i]==1){
		spec.exclude<-unlist(cluster$clusterlist[cluster$nc])
		matrix.to.cluster2<-matrix.to.cluster[spec.exclude==F,]
		matrix.to.cluster.dist2<-dist(matrix.to.cluster2,method="euclidean")
		cluster<-disthclustCBI(matrix.to.cluster.dist2,method="ward",cut="number",
		k=final.result.cluster.boot$clust.stable[i]+final.result.cluster.boot$clust.instable[i])
	}

	spec.groups<-cluster$partition	#which species belong to which group
	

############################################################
#2. For each stable cluster-group, calculate models
#############################################################
	#vectors and matrices for storing output
	out.group.r2<-matrix(nrow=ncol(env),ncol=length(no.stable.groups))        
	rownames(out.group.r2)<-colnames(env)                                     
	colnames(out.group.r2)<-no.stable.groups
	out.group.shape<-out.group.r2
	out.group.lr<-out.group.r2
	out.group.coeff<-rbind(out.group.r2,out.group.r2)
	rownames(out.group.coeff)<-c(rownames(out.group.r2),
		sapply(rownames(out.group.r2),function(x) paste("I(",x,"^2)",sep="")) )
	av.r2<-vector(length=length(no.stable.groups))
	beta.mean.groups<-vector("list",length(no.stable.groups))
	spec.modelled<-vector(mode="character",length=length(no.stable.groups))
	no.spec.modelled<-vector(mode="numeric",length=length(no.stable.groups))
	names(av.r2)<-names(beta.mean.groups)<-names(spec.modelled)<-names(no.spec.modelled)<-
		as.character(no.stable.groups)

	#start loop over all groups j
	for(j in no.stable.groups[1]:no.stable.groups[length(no.stable.groups)]){	#2o(j, group loop,
	if(is.element(j,no.stable.groups)) {							#only if group is a staple group)
		#print(paste("..GROUP - ",j," -",sep=""))		
		try(detach(data),silent=T) #avoid conflicting definitions of variables
		

		#################################
		#prepare DATA
		################################
		#print(paste(".....Group ",j," - prepare DATA",sep=""))

		#extract which species of the group are actually
		#present in the data
		##########################
		group<-names(spec.groups[spec.groups==j])	#names of all species in group
							#(from overall clustering - not necessarily present on
							#this particular plot)
		spec.group.plot<-spec.plot[is.element(spec.plot,group)]	#names of species of the group found 
												#on the plot
		#go to next group, if no species of the current group are present in the data
		if(length(spec.group.plot)==0){						#3o
			av.r2[as.character(j)]<-NA			# assign NA in output.file
			beta.mean.groups[as.character(j)]<-NA	## assign NA in output.file
			#print("...............no species of this group are present in the plot data")
			next
		}											#3c
		spec.group.covers<-spec.plot.covers[spec.group.plot,]		#covers of group species on plot
		
		spec.modelled[as.character(j)]<-paste(spec.full.names[spec.group.plot,],collapse=", ") #full names of group species
		no.spec.modelled[as.character(j)]<-length(spec.group.plot)

		#calculate cumulative covers for all species of the group,
		#scaled between 0 and 1 
		###########################		
		cum.covers<-rep(0,ncol(spec.group.covers))

		for(c in 1:ncol(spec.group.covers)){					#3o
			covers<-spec.group.covers[,c][spec.group.covers[,c]>0]
			if(is.na(sum(covers)))cum.covers[c]<-NA
			else{									#4o
				if(sum(covers)==0) cum.covers[c]<-0
				if(sum(covers)>0){					#5o
					if(length(covers)==1) temp<-covers
					else{							#6o
					   temp<-covers[1]
					   for(r in 2:length(covers)){ 		#7o
						temp<-temp+(1-temp)*covers[r]
					   }							#7c
					}							#6c
					cum.covers[c]<-temp
				}								#5c
			}									#4c
		}										#3c			

	
		
#		#stop loop if cumulative coverages of the group do not exceed min.cum.cov in any plot	
#		if(max(round(cum.covers,digits=5),na.rm=T)<min.cum.cov){
#			print(".......group does not reach cumulative coverages > min.cum.cov in any plot" )
#			av.r2[as.character(j)]<-NA
#			beta.mean.groups[as.character(j)]<-NA
#			next
#		}
#		#stop loop if prevalence of group (= no. of occurences) < min.prev
#		if(sum(round(cum.covers,digits=5)>0,na.rm=T)/nrow(plots.env)< min.prev){
#			print(".......prevalence of group is too low" )
#			av.r2[as.character(j)]<-NA
#			beta.mean.groups[as.character(j)]<-NA
#			next
#		}
                 
		data<-cbind(cum.covers,env)	    
		rm(cum.covers)
	
		data<-na.exclude(data)		#all data for group j of cluster i
		data<-data.frame(data)
		attach(data)			#necessary so that variable names will be known
		
		y<-round(cbind(data[,1],1-data[,1])*100,digits=0)	
                                
		#stop loop if cumulative coverages of the group do not exceed min.cum.cov in any plot	
		if(max(y[,1])<min.cum.cov){
			#print(".......group does not reach cumulative coverages > min.cum.cov in any plot" )
			av.r2[as.character(j)]<-NA
			beta.mean.groups[as.character(j)]<-NA
			next
		}

		#stop loop if prevalence of group (= no. of occurences) < min.prev
		if(sum(y[,1]>0)/nrow(env)< min.prev){                        
			#print(".......prevalence of group is too low" )
			av.r2[as.character(j)]<-NA
			beta.mean.groups[as.character(j)]<-NA
			next
		}

		
		##############################################
		#calculate UNIVARIATE MODELS, loop over all explanatory variables (v)
		##############################################
		#print(paste(".....Group ",j," - Univariate models",sep=""))

		for (v in 1:ncol(env)){					#3o (v, variable loop)         
			r2<-"xxx"
			shape<-"xxx"	
			lr<-"xxx"
			x<-data[,colnames(env)[v]]                            
			name.x<-colnames(env)[v]                        
		    #print(paste(".....Variable ",name.x,sep=""))
			#sigmoidal response
		    #---------------------------------
			m.s<-glm(y~x,family=binomial)
			coeff.p.s<-summary(m.s)$coefficients[2,4]
			b1.s<-m.s$coefficients[2]
			r2.s<-(cor(predict(m.s,type="response"),y[,1]/100,method="pearson")^2)
			if(is.na(r2.s)) r2.s<-0
			lr.s<-lrtest(m.s)$"Pr(>Chisq)"[2]

		    #unimodal response	
		    #-----------------------------------		
			m.u<-glm(y~x+I(x^2),family=binomial)
			#auto.u<-bgtest(y~x+I(x^2), order = 1, type = "Chisq")$p.value
			if(length(summary(m.u)$coefficients[,4])<3) {		#4o
			   r2.u<-0
			   coeff.p.u<-1
			   b1.u<-0
			   b2.u<-0
			} else {								#4c/4o
			   coeff.p.u<-summary(m.u)$coefficients[2:3,4]
			   b1.u<-m.u$coefficients[2]
			   b2.u<-m.u$coefficients[3]
			}									#4c
			r2.u<-(cor(predict(m.u,type="response"),y[,1]/100)^2)
			if(is.na(r2.u)) r2.u<-0
			lr.u<-lrtest(m.u)$"Pr(>Chisq)"[2]	

		   #compare unimodal and sigmoidal response
		   #---------------------------------------
		   #case1: UNIMOLDAL response, hump shaped, 
		   #x with positive coefficient, x^2 with negative coefficient OR
		   #x with negative coefficient, x^2 with negative coefficient  
  			if(all(lr.u<=0.05,r2.u>r2.s,max(coeff.p.u)<=0.05,b1.u>0,b2.u<0)|| #4o
			   all(lr.u<=0.05,r2.u>r2.s,max(coeff.p.u)<=0.05,b1.u<0,b2.u<0)){				   
			   form<-"unimodal"  
			   shape<-"hump"  
			   lr<-lr.u	
			   r2<-r2.u			   
	   		   out.group.coeff[name.x,as.character(j)]<-m.u$coefficients[2]
			   out.group.coeff[paste("I(",name.x,"^2)",sep=""),as.character(j)]<-m.u$coefficients[3]
			 }											#4c
		   #case2: SIGMOIDAL response
			if(lr.s<=0.05&&r2.u<r2.s&&coeff.p.s<=0.05){				#4o
			   form<-"sigmoidal"
			   out.group.coeff[name.x,as.character(j)]<-m.s$coefficients[2]
			   if(b1.s<0) shape<-"decreasing"
			   if(b1.s>0) shape<-"increasing"
			   r2<-r2.s
			   lr<-lr.s
			 }											#4c
		   #case3: r2 of unimodal response larger than of sigmoidal response,
		   #but coefficients of unimodal response not significant
		   #--> sigmoidal response
		 	if(lr.s<=0.05&&r2.u>r2.s&&coeff.p.s<=0.05&&max(coeff.p.u)>0.05){	#4o
			   form<-"sigmoidal"
			   out.group.coeff[name.x,as.character(j)]<-m.s$coefficients[2]
			   if(b1.s<0) shape<-"decreasing"
			   if(b1.s>0) shape<-"increasing"
			   r2<-r2.s
			   lr<-lr.s
			 }											#4c
		   #case4: r2 of unimodal response larger than of sigmoidal response,
		   #but bowl shaped
		   #--> sigmoidal response, if coefficients are significant
			#if(lr.u<=0.05&&r2.u>r2.s&&max(coeff.p.u)<=0.05&&b1.u<0&&b2.u>0)
			if(all(lr.u<=0.05,r2.u>r2.s,max(coeff.p.u)<=0.05,b1.u<0,b2.u>0,
			       coeff.p.s<=0.05,lr.s<=0.05) ||
			       all(lr.u<=0.05,r2.u>r2.s,max(coeff.p.u)<=0.05,b1.u>0,b2.u>0,
				 coeff.p.s<=0.05,lr.s<=0.05)){
			   form<-"sigmoidal"
			   out.group.coeff[name.x,as.character(j)]<-m.s$coefficients[2]
			   if(b1.s<0) shape<-"decreasing"
			   if(b1.s>0) shape<-"increasing"		  
 			   r2<-r2.s
			   lr<-lr.s
			 }
		   #case5: univariate model not significant
			 if(coeff.p.s>0.05&&max(coeff.p.u)>0.05) r2<-"n.s."
			 if(lr.s>0.05&&lr.u>0.05) r2<-"n.s."
   
		   #write results to output matrices	  	    
			 out.group.r2[v,as.character(j)]<-r2
			 out.group.shape[v,as.character(j)]<-shape
			 out.group.lr[v,as.character(j)]<-lr
		}											#3c (v, variable loop)

		#excude significant variables with r2 below threshold
		out.group.r2[out.group.r2<r2.min]<-"xxx"
		out.group.shape[out.group.r2<r2.min]<-"xxx"
		out.group.lr[out.group.r2<r2.min]<-"xxx"											
		
		#catch the case of all univariate models not being significant
		if(sum(out.group.r2[,as.character(j)]!="n.s."&out.group.r2[,as.character(j)]!="xxx")==0){
			av.r2[as.character(j)]<-0
			next
		}		

		##########################################
		#caclulate MULTIPLE MODELS for each group
		##########################################
		#print(paste(".....Group ",j," - Multiple models",sep=""))
		
		#extract significant variables for each group, extract variable names
            n.var<-sum(out.group.shape[,as.character(j)]!="xxx") #number of sig. variables
		names.var<-rownames(out.group.shape[out.group.shape[,as.character(j)]!="xxx",,drop=FALSE])
		names.var.1<-names.var #variable names without ^2
		
		#names of squared variables
		names.var.2<-vector(mode="character", length=length(names.var))
		for (v in 1:length(names.var.1)){						#3o
			names.var.2[v]<-paste(c("I(",names.var.1[v],"^2)"),collapse="")
		}											#3c
		form.var<-out.group.shape[out.group.shape[,as.character(j)]!="xxx",as.character(j)]
		for(v in 1:length(names.var)){						#3o #variables that need to be squared
			if(form.var[v]=="hump") names.var[v]<-paste(names.var[v],"+I(",names.var[v],"^2)",sep="")
		}											#3c
			
   
		#function to create all subsets with combinations of 1 to max.var variables
		subsets<-function(n, k, set = names.var){					#3o
			if(k <= 0) NULL
			else if(k >= n) set
			else rbind(cbind(set[1], Recall(n-1, k-1, set[-1])), Recall(n-1, k, set[-1]))
		}											#3c

		#create list of variable combinations
		var.combinations<-vector("list")
		for(combi in 1:max.var){							#3o
			temp<-subsets(n.var,max.var-combi+1)
			if(is.vector(temp)==T) temp<-t(as.matrix(temp))
			temp<-apply(temp,1,as.list)
			var.combinations<-c(var.combinations,temp)		
		}											#3c

		#create formulas
		formula<-vector(length=length(var.combinations),mode="character")
		for(combi in 1:length(var.combinations)){					#3o
			temp<-unlist(var.combinations[combi],use.names=F)
			formula[combi]<-paste(temp,collapse="+")
		}											#3c
		formula<-paste("y~",formula)

		out.mult.mod<-matrix(ncol=7+2*length(names.var.1),nrow=length(formula))
		colnames(out.mult.mod)<-c("formula","r2","aicc","max.p.coeff","good.mod","lr",
			"(Intercept)",names.var.1,names.var.2)
    
    
		#calculate model for each formula	
		################################################
		for(f in 1:length(formula)){							#3o (f, formula loop)
			m<-glm(formula=as.formula(formula[f]),family=binomial)
			p.coeff<-max(summary(m)$coefficients[,4][-1])
			
			#calculate aicc
			k<-length(m$coefficients)
			n<-length(m$fitted.values)
			aicc<-m$aic + ((2*k*(k+1))/(n-k-1))
			r2<-(cor(predict(m,type="response"),y[,1]/100,method="pearson")^2)

			#do LR-Test for all models with one variable less
			if(length(unlist(var.combinations[f]))>1) {			#4o
			   for (v in 1:length(unlist(var.combinations[f]))){		#5o
			   	temp<-unlist(var.combinations[f])[-v]
			   	formula.temp<-paste(temp,collapse="+")
			   	formula.temp<-paste("y~",formula.temp)
			   	m.temp<-glm(formula=as.formula(formula.temp),family=binomial)
			   	lr<-lrtest(m,m.temp)[2,5]
			   	if(lr>0.05) out.mult.mod[f,6]<-"n.s."
			   }										#5c
			}										#4c
      
			#write results to output matrix
			out.mult.mod[f,1]<-formula[f]
			out.mult.mod[f,2]<-r2
			out.mult.mod[f,3]<-aicc
			out.mult.mod[f,4]<-p.coeff
			
			#check model significance
			if(p.coeff<=0.05 && lr!="n.s.") out.mult.mod[f,5]<-1
			for (koeff in 1:length(m$coefficients)){				#4o
			   out.mult.mod[f,names(m$coefficients[koeff])]<-m$coefficients[koeff]
			}										#4c
			#check form of relationships (should be the same as in univariate models)
			for(koeff in 2:length(m$coefficients)){

			   if(all(m$coefficients[koeff]<0, 
				is.na(charmatch("I(",names(m$coefficients)[koeff]))==F,
				out.group.coeff[names(m$coefficients)[koeff],as.character(j)]>0))			
			   	out.mult.mod[f,5]<-"n.s."
			   if(all(m$coefficients[koeff]>0, 
				out.group.coeff[names(m$coefficients)[koeff],as.character(j)]<0,
				is.na(charmatch("I(",names(m$coefficients)[koeff])==T)))
				out.mult.mod[f,5]<-"n.s."
			   if(all(m$coefficients[koeff]<0, 
				out.group.coeff[names(m$coefficients)[koeff],as.character(j)]>0,
				is.na(charmatch("I(",names(m$coefficients)[koeff])==T)))
				out.mult.mod[f,5]<-"n.s."
			} 	
		}											#3c (f, formula loop)			
		
		mult.mod.sig<-out.mult.mod[which(out.mult.mod[,5]==1),] #significant multiple models

		#catch the case of only 1 significant model
		if(is.vector(mult.mod.sig)==T) mult.mod.sig<-t(as.matrix(mult.mod.sig))
  
		#########################################
		#MODEL AVERAGING for significant models (following the method of Burnham)
		#########################################
		#print(paste(".....Group ",j," - Model averaging",sep=""))
		#sort models according to their aicc
		############################################# ERROR found and corrected ##########################################
    #if(nrow(mult.mod.sig)>1) sorted.aicc<-mult.mod.sig[order(mult.mod.sig[,3],decreasing=F),]
		if(nrow(mult.mod.sig)>1) sorted.aicc<-mult.mod.sig[order(as.numeric(mult.mod.sig[,3],decreasing=F)),]
    #################################################################################################################
    if(nrow(mult.mod.sig)==1) sorted.aicc<-mult.mod.sig
		#calculate delta.i between models and aicc best model
		delta.i <- vector(mode="numeric",length=nrow(sorted.aicc)) 
		if(nrow(sorted.aicc)>1){							#3o
			for (m in 2:nrow(sorted.aicc)){					#4o
			   delta.i[m]<-as.numeric(sorted.aicc[m,3]) - as.numeric(sorted.aicc[1,3])
			}										#4c
		}											#3c
		if(nrow(sorted.aicc)==1) delta.i[1]<-0

		likelihood<-exp(-0.5*delta.i) #calculate, from delta.i, likelihood for the give data(Burnham p. 74)
		wi<-likelihood/sum(likelihood) #and from the likelihoods the Akaike weights (Burnham p. 75)

		#calculate model averaged coefficients (beta.mean)
		beta.mean<-vector(mode="numeric") 
		weight.var<-vector(mode="numeric") 
		coeffs<-sorted.aicc[,6:ncol(sorted.aicc)]
		mode(coeffs)<-"numeric"
		if(nrow(sorted.aicc)>1) coeffs<-coeffs[,which(colSums(coeffs,na.rm=T,dims=1)!=0)]
		if(nrow(sorted.aicc)==1) coeffs<-as.matrix(t(coeffs[which(coeffs!="NA")]))

		if(nrow(sorted.aicc)>1){							#3o
			for (c in 1:ncol(coeffs)){						#4o
			   temp<-na.exclude(data.frame(wi,as.numeric(coeffs[,c])*wi))
			   if(mean(temp[,1])!= "NaN") beta.mean[c]<-sum(temp[,2])#/nrow(sorted.aicc)
			   weight.var[c]<-sum(temp$wi)
			}										#4c
		}											#3c
		if(nrow(sorted.aicc)==1) {							#3o
			weight.var<-rep(1,times=length(coeffs))
			beta.mean<-round(as.vector(coeffs),digits=4)
		}											#3c
     		#calculate predicted values for averaged model
		if(nrow(sorted.aicc)>1) {							#3o
			products<-matrix(nrow=nrow(data),ncol=ncol(coeffs)-1) 
			colnames(products)<-colnames(coeffs)[-1]
		}											#3c
		if(nrow(sorted.aicc)==1) {							#3o
			products<-matrix(nrow=nrow(data),ncol=length(coeffs)-1)
			colnames(products)<-rownames(as.matrix(coeffs))[-1]
		}											#3c		
		data.1<-data[,names.var.1]	#data without ^2
		data.2<-cbind(data.1,data.1^2)#data ^2
		colnames(data.2)<-c(names.var.1,names.var.2)
		data.2<-data.2[,colnames(coeffs)[-1],drop=FALSE]
		weight.var[-1]<-round(weight.var[-1]*100/sum(weight.var[-1]),digits=4)
		beta.mean<-cbind(beta.mean,weight.var)
		colnames(beta.mean)<-c("beta.mean","weight.var[%]")
		rownames(beta.mean)<-c("Intercept",colnames(data.2))

		for (c in 1:ncol(products)){							#3o
			products[,c]<-data.2[,c]*beta.mean[c+1,1]
		}											#3c
		sum.products<-rowSums(products,na.rm=T,dims=1)+beta.mean[1,1]
		pred<-100/(1+exp(-1*sum.products)) 
		obs<-y[,1]
		r2<-(cor(pred,obs,method="pearson"))^2	#r2 of averaged model
		av.r2[as.character(j)]<-r2
		beta.mean.groups[as.character(j)]<-list(beta.mean)
		
	}												#2c (j, group loop,
	}												#and if-condition for stable groups only)

	#write results to output list
	mod.clust$final.combi[[i]]<-final.combinations[i]
	mod.clust$r2.groups[[i]]<-av.r2	
	mod.clust$r2.av[[i]]<-mean(av.r2,na.rm=T)		#maybe consider using a weighted mean here?
									#with weight=no.spec.modelled in the respective group
	mod.clust$min.r2[[i]]<-min(av.r2,na.rm=T)
	mod.clust$beta.mean[[i]]<-beta.mean.groups
	mod.clust$spec.modelled[[i]]<-as.list(spec.modelled)
	mod.clust$no.modelled.group[[i]]<-no.spec.modelled
	mod.clust$no.modelled.clustering[[i]]<-sum(no.spec.modelled,na.rm=T)
}													#1c (i, cluster loop)

#write objects that are needed for further analysis into file
dump(list=c("mod.clust","spec.plots.t","env"),file="modelling.output.txt")    
}													#0c (missing trait info)



if(stop==1) warning("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Trait values are missing for at least 1 species:
		Names of the missing species were exported to the file 'missing.txt'. 
		Missing information needs to be completed before you can proceed.
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

#clear workspace
#rm(list = ls())








