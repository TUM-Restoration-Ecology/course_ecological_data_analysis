####################################################
### R script (Version 2.9.0) to run multivariate analyses of   
### species and community functional responses to environmental gradients
### See paper Michael Kleyer, Francesco de Bello, Stéphane Dray, Jan Lepš, 
### Tonia Meier, Robin J. Pakeman, Barbara Strauss, Wilfried Thuiller, Sandra Lavorel (20XX):
### Assessing species and community functional responses to environmental gradients: which multivariate methods?
###################################################
### Version 17/06/2009
### Graphics produced by this script may deviate slightly from those in the tutorial (Appendix S2)
#######################################################
### chunk number 2: pre1
###################################################
traits <- read.table(file="data/Species_traits.txt", sep="\t")
spe<- read.table(file="data/Site_species.txt", sep="\t")
env<-read.table("data/Site_env.txt", sep="\t")

oldpar<-par()   #used to reset graphical parameters 

###################################################
### chunk number 3: pre2
###################################################
source("scripts/Inference_modelset.r")
source("scripts/Inference_compute.r")
source("scripts/corratio.R")
source("scripts/calinski.R")
source("scripts/VarScoreOMI.r")
source("scripts/doublerda.R")


###################################################
### chunk number 4: pre3
###################################################
library(ade4)
library(MASS)
library(vegan)
library(ecodist)
library(maptools)
library(rpart)
library(splines)
library(gam)         
library(pgirmess)     
library(utils)
library(combinat)
library(mvpart)
library(cluster)
library(fpc)
library(clusterSim)
library(lmtest)
library(Hmisc)
library(gplots)


###################################################
### chunk number 5: cca1
###################################################
coa1 <- dudi.coa(spe, scannf = F)
cca1 <- pcaiv(coa1, env, scannf = F)


###################################################
### chunk number 6: cca2
###################################################
100 * sum(cca1$eig) / sum(coa1$eig)


###################################################
### chunk number 7: cca3
###################################################
s.label(cca1$c1, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(cca1$c1,row.names(cca1$c1), cex=0.7)
s.arrow(cca1$cor[-1,], add.plot=TRUE)


###################################################
### chunk number 8: cwm1
###################################################
cwm.tab <- prop.table(as.matrix(spe),1)%*%as.matrix(scale(traits))


###################################################
### chunk number 9: cwm2
###################################################
pca.cwm <- dudi.pca(cwm.tab,scannf=FALSE)
rda.cwm <- pcaiv(pca.cwm,env, scannf=FALSE)


###################################################
### chunk number 10: cwm3
###################################################
100 * sum(rda.cwm$eig) / sum(pca.cwm$eig)


###################################################
### chunk number 11: cwm4
###################################################
pca.cwm <- dudi.pca(cwm.tab,scannf=FALSE)
100 * rda.cwm$eig / sum(rda.cwm$eig)


###################################################
### chunk number 12: cwm5
###################################################
s.arrow(rda.cwm$c1, xlim=c(-1,1), boxes = FALSE)
s.label(rda.cwm$cor[-1,], add.plot=T, clab=1.5)


###################################################
### chunk number 13: clusmod1
###################################################
max.traits<-6		
min.cophenetic.corr<-0.7
max.no.groups<-10
noisecut.value<-5
min.prop<-0.9
no.boot<-200
mean.jac<-0.7
name.output.table<-"result.cluster.boot.txt"
traits.to.consider<-"all"


###################################################
### chunk number 14: clusmod2
###################################################
cor(traits)


###################################################
### chunk number 15: clusmod3
###################################################
source("scripts/script1_cluster_analysis_report.r")


###################################################
### chunk number 16: clusmod3b
###################################################
clus.boot<-read.table("result.cluster.boot.txt",sep="\t",header=TRUE)
clus.boot[,c(1:6,13:14)]


###################################################
### chunk number 17: clusmod4
###################################################
max.var<-3 
r2.min<-0	
min.cum.cov<-0.5
min.prev<-0.1		


###################################################
### chunk number 18: clusmod5
###################################################
source("scripts/script2_group_models_report.r")


###################################################
### chunk number 19: clusmod6
###################################################
name.outputfile.1<-"modelling.output.groupwise.txt"
name.outputfile.2<-"modelling.output.clusterwise.txt"
source("scripts/script3_output_tables_report.r")
mod.groupwise<-read.table(name.outputfile.1,sep="\t",header=TRUE)
mod.clusterwise<-read.table(name.outputfile.2,sep="\t",header=TRUE)
mod.clusterwise[order(-mod.clusterwise$no.clusters,-mod.clusterwise$r2.av),c(1:9,17)]
mod.groupwise[,1:3]


###################################################
### chunk number 20: clusmod7
###################################################
combis.to.plot<-c(25)
min.weight<-10
mycol<-palette()
plot.points<-0
name.table<-"spec.groups.txt"		
source("scripts/script4_group_plots_report.r")
par(oldpar)

###################################################
### chunk number 21: clusmod8
###################################################
name.table<-"spec.groups.txt"		#name of output file
#combis.to.compare<-as.vector(output.table.2[,"no.combi"])
combis.to.compare<-c(25,24,44,26,12,17)
source("scripts/script5_group_assignements_report.r")
spec.in.groups<-read.table("spec.groups.txt",sep="\t",header=TRUE)
spec.in.groups


###################################################
### chunk number 22: rdasrta1
###################################################
rda.phosp <- rda(X=spe, Y=env$SOIL.P, Z=env$dist.int, scale = TRUE)
rda.dist <- rda(X=spe, Y=env$dist.int, Z=env$SOIL.P, scale = TRUE)
rda.both <- rda(X=spe, Y=env[,c("dist.int","SOIL.P")], scale = TRUE)


###################################################
### chunk number 23: rdasrta2
###################################################
100 * rda.phosp$CCA$tot.chi/rda.phosp$tot.chi
100 * rda.dist$CCA$tot.chi/rda.dist$tot.chi


###################################################
### chunk number 24: rdasrta2b
###################################################
cor(env$SOIL.P,env$dist.int)
cor(rda.phosp$CCA$v,rda.dist$CCA$v)


###################################################
### chunk number 25: rdasrta2c
###################################################
rda.phosp$CCA$biplot


###################################################
### chunk number 26: rdasrta3
###################################################
df.phosp <- cbind(RDA1=rda.phosp$CCA$v[,1],traits)
df.dist <- cbind(RDA1=rda.dist$CCA$v[,1],traits)
rta.dist<-rpart::rpart(RDA1~.,data=df.dist, xval = 100, minbucket = 3)
rta.phosp<-rpart::rpart(RDA1~.,data=df.phosp, xval = 100, minbucket = 3)


###################################################
### chunk number 27: rdasrta4
###################################################
rta.phosp


###################################################
### chunk number 28: rdasrta5
###################################################
plotcp(rta.phosp)


###################################################
### chunk number 29: rdasrta6
###################################################
idx.min <- which.min(rta.phosp$cptable[,4])
rta.phosp.pruned <- prune(rta.phosp, cp = rta.dist$cptable[,1][idx.min])
plot(rta.phosp.pruned, ylim = c(0.28,1.1))
text(rta.phosp.pruned, use.n=T)
arrows(0.9, 0.28, 2, 0.28, length = 0.1)
text(1.5, 0.26, "Species score on RDA axis (phosphorous)", cex = 1)


###################################################
### chunk number 30: rdasrta7
###################################################
table(rta.phosp.pruned$where)


###################################################
### chunk number 31: rdasrta8
###################################################
100 * (1-rta.phosp$cptable[2,3])
100 * (1-rta.phosp$cptable[3,3])


###################################################
### chunk number 32: rdasrta9
###################################################
rda.dist$CCA$biplot


###################################################
### chunk number 33: rdasrta10
###################################################
rta.dist


###################################################
### chunk number 34: rdasrta11
###################################################
plotcp(rta.dist)
cptab <- rta.dist$cptable


###################################################
### chunk number 35: rdasrta12
###################################################
idx.min <- which.min(cptab[,4])
rta.dist.pruned <- prune(rta.dist, cp = cptab[,1][idx.min])
plot(rta.dist.pruned, ylim = c(0.33,1.1))
text(rta.dist.pruned, use.n=T)
arrows(0.9, 0.33, 5, 0.33, length = 0.1)
text(3, 0.31, "Species score on RDA axis (grazing)", cex = 1)


###################################################
### chunk number 36: rdasrta13
###################################################
idx.min2 <- min((1:nrow(cptab))[cptab[,4]<min(cptab[,4]+cptab[,5])])
rta.dist.pruned2 <- prune(rta.dist, cp = cptab[,1][idx.min2])
plot(rta.dist.pruned2, ylim = c(0.33,1.1))
text(rta.dist.pruned2, use.n=T)
arrows(0.9, 0.33, 2, 0.33, length = 0.1)
text(1.5, 0.31, "Species score on RDA axis (grazing)", cex = 1)


###################################################
### chunk number 37: rdasrta14
###################################################
100 * (1-cptab[idx.min,3])
100 * (1-cptab[idx.min2,3])


###################################################
### chunk number 38: rdamrta1
###################################################
rda.both <- rda(X=spe, Y=env[,c("dist.int","SOIL.P")], scale = TRUE)


###################################################
### chunk number 39: rdamrta2
###################################################
100 * rda.both$CCA$tot.chi/rda.both$tot.chi


###################################################
### chunk number 40: rdamrta3
###################################################
100 * rda.both$CCA$eig/rda.both$CCA$tot.chi


###################################################
### chunk number 41: rdamrta4
###################################################
plot(rda.both, display=c("bp","sp"))


###################################################
### chunk number 42: rdamrta5
###################################################
df.both <- cbind(RDA1=rda.both$CCA$v.eig[,1],RDA2=rda.both$CCA$v.eig[,2],traits)
rta.both<-mvpart(data.matrix(df.both[,1:2])~Polycarpic+Cnratio+seed.mass.log
                 +SLA+height+Onset.flower,data=df.both, xval = 100, xv="1se")


###################################################
### chunk number 43: omi1
###################################################
pca.env<-dudi.pca(env, scannf=FALSE)
scatter(pca.env)


###################################################
### chunk number 44: omi2
###################################################
100 * pca.env$eig/sum(pca.env$eig)


###################################################
### chunk number 45: omi3
###################################################
pca.env$co


###################################################
### chunk number 46: omi4
###################################################
omi1<-niche(pca.env, spe, scannf=FALSE)
plot(omi1)


###################################################
### chunk number 47: omi5
###################################################
100 * omi1$eig/sum(omi1$eig)


###################################################
### chunk number 48: omi6
###################################################
s.arrow(omi1$c1, clab = 0.8, xlim=c(-2.5,2.5))
s.label(omi1$li, xax = 1,yax = 2, clabel=0,add.plot = TRUE)	
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(omi1$li, rownames(omi1$li), cex=0.7)


###################################################
### chunk number 49: omi7
###################################################
par(mfrow=c(1,2))
sco.distri(omi1$ls[,1],spe,clab=0.7)
sco.distri(omi1$ls[,2],spe,clab=0.7)

par(mar = c(5, 4, 4, 2) + 0.1)
par(oldpar)
###################################################
### chunk number 50: omi8
###################################################
traits[,1]<-as.factor(traits[,1])
modelset<-Inference_modelset(Explanatory=traits)
inf.axis1 <- Inference_compute(Fam="gaussian",  combin=modelset[[1]], Mat=modelset[[2]], 
  Response=omi1$li[,1], Explanatory=traits, Average = TRUE)
inf.axis2 <- Inference_compute(Fam="gaussian",  combin=modelset[[1]], Mat=modelset[[2]], 
  Response=omi1$li[,2], Explanatory=traits, Average = TRUE)


###################################################
### chunk number 51: omi9
###################################################
dd.names <- c('Poly- carpic','CN ratio', 'Log(SM)',
 'Height','Onset flowering', 'SLA')
dd.names.2 <- sapply(dd.names, function(x) gsub("\\s", "\\\n", x))

barplot(inf.axis1$Var.importance[,1], names.arg=dd.names.2)


###################################################
### chunk number 52: omi10
###################################################
Limits <- apply(inf.axis1$Plot.response[, seq(2, 12, by=2)], 2, range)
lim <- c(min(Limits[1,]), max(Limits[2,]))
par(mfrow=c(3,2))
plot(as.factor(inf.axis1$Plot.response[,1]), inf.axis1$Plot.response[,2],
     ylim=lim, type="l", xlab="Policarpic", ylab="Species position axis 1") 
plot(inf.axis1$Plot.response[,3], inf.axis1$Plot.response[,4], ylim=lim,  
     type="l", xlab="CN ratio", ylab="Species position on OMI axis 1")
plot(inf.axis1$Plot.response[,5], inf.axis1$Plot.response[,6], ylim=lim,  
     type="l", xlab="Log (Seed mass)", ylab="Species position on OMI axis 1")
plot(inf.axis1$Plot.response[,9], inf.axis1$Plot.response[,10], ylim=lim, 
     type="l", xlab="Height", ylab="Species position on OMI axis 1")
plot(inf.axis1$Plot.response[,11], inf.axis1$Plot.response[,12], ylim=lim, 
     type="l", xlab="Onset of flowering", ylab="Species position on OMI axis 1")
plot(inf.axis1$Plot.response[,7], inf.axis1$Plot.response[,8], ylim=lim, 
     type="l",  xlab="SLA", ylab="Species position on OMI axis 1")
par(oldpar)

###################################################
### chunk number 53: omi11
###################################################
Averaged.Pred.1.2<-cbind(inf.axis1$Averaged.Pred, inf.axis2$Averaged.Pred)
hc1 <- hclust(dist(Averaged.Pred.1.2), method = "ward")
plot(hc1) 


###################################################
### chunk number 54: omi12
###################################################
ntest <- 6
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc1, k = i)
  res[i-1] <- calinski(tab=Averaged.Pred.1.2, fac = fac)[1]
}
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")

par(oldpar)

###################################################
### chunk number 55: omi13
###################################################
nbgroup <- 3
spe.group <- as.factor(cutree(hc1, k = nbgroup))
spe.group <- as.factor(spe.group)
s.class(Averaged.Pred.1.2, spe.group, col= 1:nlevels(spe.group)) 
s.arrow(omi1$c1, xax=1, yax=2, csub = 1, clab = 0.8, add.plot=T)


###################################################
### chunk number 56: omi14
###################################################
eta2 <- cor.ratio(traits[,-1], data.frame(spe.group), weights = rep(1, length(spe.group)))
par(mfrow=n2mfrow(ncol(traits)))
plot(table(spe.group,traits[,1]), main =names(traits)[1])
for(i in 2:ncol(traits)){
  label <- paste(names(traits)[i], "(cor.ratio =", round(eta2[i-1],3), ")")
  plot(traits[,i]~spe.group, main = label, border = 1:nlevels(spe.group))
}
par(oldpar)

###################################################
### chunk number 57: omi15
###################################################
traits[,1]=as.numeric(traits[,1])


###################################################
### chunk number 58: rlq0
###################################################
pca.traits <- dudi.pca(traits, row.w = coa1$cw, scannf = FALSE)
pca.env <- dudi.pca(env, row.w = coa1$lw, scannf = FALSE)


###################################################
### chunk number 59: rlq1
###################################################
rlq1 <- rlq(pca.env, coa1, pca.traits, scannf = FALSE)
summary(rlq1)


###################################################
### chunk number 60: rlq1b
###################################################
plot(rlq1)


###################################################
### chunk number 61: rlq2
###################################################
## Percentage of co-Inertia for each axis
100*rlq1$eig/sum(rlq1$eig)


###################################################
### chunk number 62: rlq3
###################################################
## weighted correlations axes / env.
t(pca.env$tab)%*%(diag(pca.env$lw))%*%as.matrix(rlq1$mR)
## weighted correlations axes / traits.
t(pca.traits$tab)%*%(diag(pca.traits$lw))%*%as.matrix(rlq1$mQ)
## correlations traits / env.
rlq1$tab


###################################################
### chunk number 63: rlq4
###################################################
s.arrow(rlq1$c1, xlim=c(-1,1), boxes = FALSE)
s.label(rlq1$li, add.plot=T, clab=1.5)


###################################################
### chunk number 64: rlq4b
###################################################
s.label(rlq1$lQ, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(rlq1$lQ,row.names(rlq1$lQ), cex=0.7)


###################################################
### chunk number 65: rlq5
###################################################
hc2 <- hclust(dist(rlq1$lQ), method = "ward")
plot(hc2)


###################################################
### chunk number 66: rlq6
###################################################
ntest <- 6
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- calinski(tab=rlq1$lQ, fac = fac)[1]
}

par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")

par(oldpar)
###################################################
### chunk number 67: rlq7
###################################################
spe.group2 <- as.factor(cutree(hc2, k = which.max(res) +1))
levels(spe.group2) <- c("C","B","D","A")
spe.group2 <- factor(spe.group2, levels=c("A","B","C","D"))
s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))
s.arrow(rlq1$c1, add.plot = T,clab=0.8)


###################################################
### chunk number 68: rlq8
###################################################
eta2 <- cor.ratio(traits[,-1], data.frame(spe.group2), weights = rep(1, length(spe.group2)))
par(mfrow=n2mfrow(ncol(traits)))
plot(table(spe.group2,traits[,1]), main =names(traits)[1])
for(i in 2:ncol(traits)){
  label <- paste(names(traits)[i], "(cor.ratio =", round(eta2[i-1],3), ")")
  plot(traits[,i]~spe.group2, main = label, border = 1:nlevels(spe.group2))
}

par(oldpar)
###################################################
### chunk number 69: rlq9
###################################################
table(spe.group,spe.group2)


###################################################
### chunk number 70: dcca1
###################################################
dbcca1 <- dbrda(coa1,env, traits, scannf = FALSE)


###################################################
### chunk number 71: dcca2
###################################################
## percentage of explained variation by the environment
sum(cca1$eig)/sum(coa1$eig)*100
## percentage of explained variation by both traits and env.
sum(dbcca1$eig)/sum(coa1$eig)*100


###################################################
### chunk number 72: dcca3
###################################################
## Percentage of variation explained by each axis
100*dbcca1$eig/sum(dbcca1$eig)


###################################################
### chunk number 73: dcca4
###################################################
s.arrow(dbcca1$corZ[-1,], xlim=c(-1.2,1.2), boxes = FALSE)
s.label(dbcca1$corX[-1,], add.plot=T, clab=1.5)


###################################################
### chunk number 74: dcca5
###################################################
s.label(dbcca1$co, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(dbcca1$co,row.names(dbcca1$co), cex=0.7)


###################################################
### chunk number 75: dcca6
###################################################
hc3 <- hclust(dist(dbcca1$co), method = "ward")
plot(hc3)


###################################################
### chunk number 76: dcca7
###################################################
ntest <- 6
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc3, k = i)
  res[i-1] <- calinski(tab=dbcca1$co, fac = fac)[1]
}
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")

par(oldpar)

###################################################
### chunk number 77: dcca8
###################################################
nbgroup <- ifelse((which.max(res) + 1) == ntest, nlevels(spe.group2), which.max(res) + 1)
spe.group3 <- as.factor(cutree(hc3, k = nbgroup))
levels(spe.group3) <- c("B","C","D","A")
spe.group3 <- factor(spe.group3, levels=c("A","B","C","D"))
s.class(dbcca1$co, spe.group3, col = 1:nlevels(spe.group3))
s.arrow(dbcca1$corZ[-1,],add.plot=T,clab=0.8)


###################################################
### chunk number 78: dcca9
###################################################
eta2 <- cor.ratio(traits[,-1], data.frame(spe.group3), weights = rep(1, length(spe.group3)))
par(mfrow=n2mfrow(ncol(traits)))
plot(table(spe.group3,traits[,1]), main =names(traits)[1])
for(i in 2:ncol(traits)){
  label <- paste(names(traits)[i], "(cor.ratio =", round(eta2[i-1],3), ")")
  plot(traits[,i]~spe.group3, main = label, border = 1:nlevels(spe.group3))
}


###################################################
### chunk number 79: dcca10
###################################################
table(spe.group3,spe.group2)


