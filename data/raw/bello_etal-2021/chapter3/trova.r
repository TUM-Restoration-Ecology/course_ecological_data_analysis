
#########################################################################################################
# The function trova computes pairwise functional distances between a group of species, taking
# into account not only the interspecific, but also the intraspecific variability of the species.
#
# This intraspecific variability can be expressed in 4 ways:
#   1- Quantitative traits:   Providing the mean and sd for each species.
#   2- Quantitative traits:   Providing a number of observations (measurements) for each species.
#   3- Qualitative traits:    Providing the proportion of teh individuals of each species that belong
#                             to each of the categories in that trait.
#   4-Period of flowering:    Providing the day of the beginning and the day of the ending of the
#                             flowering period of each species.
#
# When the intraspecific variability is taken into account, the distance between any pair of species
# is calculated as the proportion of its range in which those species are not overlapping.
#
# Besides these 4 types of data, the function can also deal with pure categorical data (no
# intraspecific variability in the categories) and quantitative traits for which there is no
# knowledge of the intraspecific variability (only the mean is provided).
#
# INPUTS:
# species.names:     A vector containing the names of the species used.
# gaussian.data:     A matrix or data.frame containing the mean (even numbered columns) and
#               standard deviation (uneven columns) of a series of quantitative traits. The number of rows
#               of this matrix must be equal to the number of species given by species.names. The names
#               of the species are given in rownames(gaussian.data), but its order does not need
#               to be the same as that of species.names. The name of the corresponding
#               trait must be provided as the name of its column of means. NA values or standard
#               deviations less or equal than zero are not allowed in this matrix.
# kernel.data:       A matrix (NOT a data.frame) containing a number of observations of quantitative
#               traits for each species. The names of the species must be given in rownames(kernel.data),
#               and the names have to be the same as in species.names. Note that rownames(kernel.data)
#               must have repeated names (one name per observation of each species).
#               NA values are allowed in this matrix provided that there is at least a minimun number of
#               observations (defined by min.observations) for each trait on each species. The name of
#               the corresponding trait must be provided as the name of its column.
# mean.data:        A matrix or data.frame containing the mean of a series of quantitative
#               traits for each species. The number of rows of this matrix must be equal to
#               the number of species given by species.names. The names of the species are given
#               in rownames(mean.data), but its order does not need to be the same as that of
#               species.names. The name of the corresponding trait must be provided as the name of its
#               columns. NA values are not allowed in this matrix.
# multiple.category: A list containing one matrix or data frame for each of the traits that belong to
#               this category. The name of each x element of this list (multiple.category$x) must be
#               the name of the trait in question. On each matrix, the number of rows
#               must be equal to the number of species given by species.names. The names of the
#               species must be given in rownames(multiple.category$x), but its order does
#               not need to be the same as that of species.names. On each matrix, there
#               has to be a column for each possible category that the trait x can have. The value
#               expected on each cell ij of the matrix is the proportion (0-1) of individuals of the
#               species i that belong to the category j of the trait x. NA values or values above 1
#               or under 0 are not allowed. The sum of all the elements of a row must be 1.
# single.category:   A matrix or data.frame containing a unique for each species for a series of
#               cualitative traits. The number of rows of this matrix must be equal to the number
#               of species given by species.names. The names of the species must be given in
#               rownames(single.category), but its order does not need to be the same as
#               that of species.names. The name of the corresponding trait must be provided
#               as the name of its column. NA values are not allowed in this matrix
# flowering:         A matrix or data.frame containing two columns, the first with the start of the
#               flowering period of each species and the second with the end. The number of rows of
#               this matrix must be equal to the number of species given by species.names. The names
#               of the species must be given in rownames(single.category), but its order does
#               not need to be the same as that of species.names. Start and end of the flowering
#               period must be given as days (1-365), not months. NA values are not allowed
#               in this matrix.
# min.observations:  Minimun number of observations of a given trait for each species for kernel.data
#
# OUTPUTS
# The results are given in a list, that have a component for each of the traits. That component contains
# the pairwise distances between all the species given in species.names (sorted alphabetically).
# Besides, the last two components give the average of the pairwise distances, taking all traits
# into account ($gower) and the sqared root of the sum of the squared values standardized
# to be between 0 and 1 ($euclid)
#########################################################################################################

trova <- function ( species.names , gaussian.data=NULL , kernel.data=NULL , mean.data=NULL ,
                    multiple.category=NULL , single.category=NULL , flowering=NULL ,
                     min.observations=4, missing=FALSE) {
#Check wheter a minimum of inputs have been provided:
require(gtools)
species.names<-as.character(species.names)
species.names<-mixedsort(species.names)
if (is.null(species.names)){
  stop("No species names provided.")}
if (length(unique(species.names)) != length(species.names)){
  stop("Repeated names are not allowed in the species.names vector.")}
if (is.null(c(gaussian.data, kernel.data, mean.data, multiple.category, single.category, flowering))){
  stop("No data provided for the calculation of species distances.")}

n.traits <- 0
names.traits <- numeric()
trait.type <- numeric()

#Check wheter all the matrices fulfil the requirements:
#gaussian.data
if (!is.null(gaussian.data)){ #Start checking conditions for gaussian.data
  gaussian.data<-as.matrix(gaussian.data)
  gaussian.data <- gaussian.data [mixedorder(rownames(gaussian.data)),]
  if(any(rownames(gaussian.data) != species.names)){
    stop("The rownames of gaussian.data are not the same as those in the species.names vector.")}
  if(length(species.names)!= dim(gaussian.data)[1]){
    stop("For the matrix gaussian.data, only one observation per species is allowed.")}
  if(dim(gaussian.data)[2]%%2!=0) {
    stop("The matrix gaussian.data must have an even number of columns (two columns per trait)." )}
  if(any(gaussian.data[,which(1:dim(gaussian.data)[2]%%2==0)]==0,na.rm=T) & missing==FALSE){
    stop("In the matrix gaussian.data, sd cannot be 0 when the argument missing=FALSE.")}
  if(any(is.na(gaussian.data) & missing==FALSE)){
    stop("NA values are not allowed in the gaussian.data matrix when the argument missing=FALSE.")}
  if(any(is.na(gaussian.data[,which(1:dim(gaussian.data)[2]%%2!=0)]))){
    stop("NA values are not allowed in the mean values of traits in the gaussian.data matrix.")}
  n.traits <- n.traits + dim(gaussian.data)[2]/2
  names.traits <- c(names.traits, colnames(gaussian.data)[which(1:dim(gaussian.data)[2]%%2!=0)])
  trait.type <- c(trait.type, rep("gaussian", dim(gaussian.data)[2]/2))
        } #End checking conditions for gaussian.data
#kernel.data
if (!is.null(kernel.data)){ #Start checking conditions for kernel.data
  kernel.data<-as.matrix(kernel.data)
  kernel.data <- as.matrix(kernel.data [mixedorder(rownames(kernel.data)),])
  if(any(unique(rownames(kernel.data)) != species.names)){
    stop("The rownames of kernel.data are not the same as those in the species.names vector.")}
  n.traits <- n.traits + dim(kernel.data)[2]
  names.traits <- c(names.traits, colnames(kernel.data))
  trait.type <- c(trait.type, rep("kernel", dim(kernel.data)[2]))

      } #End checking conditions for kernel.data
#mean.data
if (!is.null(mean.data)){ #Start checking conditions for mean.data
  mean.data<-as.matrix(mean.data)
  mean.data <- as.matrix(mean.data [mixedorder(rownames(mean.data)),])
  if(any(rownames(mean.data) != species.names)){
    stop("The rownames of mean.data are not the same as those in the species.names vector.")}
  if(length(species.names)!= dim(as.matrix(mean.data))[1]){
    stop("For the matrix mean.data, only one observation per species is allowed.")}
  if(any(is.na(mean.data))){
    stop("NA values are not allowed in the mean.data matrix.")}
  n.traits <- n.traits + dim(mean.data)[2]
  names.traits <- c(names.traits, colnames(mean.data))
  trait.type <- c(trait.type, rep("mean", dim(mean.data)[2]))

      } #End checking conditions for mean.data
#multiple.category
if (!is.null(multiple.category)){ #Start checking conditions for multiple.category
  for (i in 1:length(multiple.category)) { #Start checking for every trait in multiple.category
    multiple.category[[i]]<-as.matrix(multiple.category[[i]])
    multiple.category[[i]] <- multiple.category[[i]] [mixedorder(rownames(multiple.category[[i]])),]
    if(any(rownames(multiple.category[[i]]) != species.names)){
      stop(paste("The rownames of the trait -", names(multiple.category)[i],
          "- in the \n list multiple.category are not the same as those in the species.names vector."))}
    if(length(species.names)!= dim(multiple.category[[i]])[1]){
      stop(paste("In the list multiple.category, only one observation per species is allowed.",
          "\n Check the species names in the trait -", names(multiple.category)[i],"-"))}
    if(any(is.na(multiple.category[[i]]))){
      stop(paste("NA values are not allowed in the multiple.category list.",
          "\n Check the values in the trait -", names(multiple.category)[i],"-"))}
    if(any(multiple.category[[i]]<0) | any(multiple.category[[i]]>1)){
       stop(paste("Values in the list multiple.category must be between 0 and 1.",
          "\n Check the values in the trait -", names(multiple.category)[i],"-"))}
    if(any(apply(multiple.category[[i]],1,sum) != 1)){
      stop(paste("In the multiple.category list, the sum of all the elements of each row must be 1.",
          "\n Check the values in the trait -", names(multiple.category)[i],"-"))}

           } #End checking for every trait in multiple.category
  n.traits <- n.traits + length(multiple.category)
  names.traits <- c(names.traits, names(multiple.category))
  trait.type <- c(trait.type, rep("multiple", length(multiple.category)))
  multiple.category[[i]] <- multiple.category[[i]] [mixedorder(rownames(multiple.category[[i]])),]
    } #End checking conditions for multiple.category
#single.category
if (!is.null(single.category)){ #Start checking conditions for single.category
  single.category<-as.matrix(single.category)
  single.category <-as.matrix( single.category [mixedorder(rownames(single.category)),] )
  if(any(rownames(single.category) != species.names)){
    stop("The rownames of single.category are not the same as those in the species.names vector.")}
  if(length(species.names)!= dim(single.category)[1]){
    stop("For the matrix single.category, only one observation per species is allowed.")}
  if(any(is.na(single.category))){
    stop("NA values are not allowed in the single.category matrix.")}
  n.traits <- n.traits + dim(single.category)[2]
  names.traits <- c(names.traits, colnames(single.category))
  trait.type <- c(trait.type, rep("single", dim(single.category)[2]))
    } #End checking conditions for single.category
#flowering
if (!is.null(flowering)){ #Start checking conditions for flowering
  flowering<-as.matrix(flowering)
  flowering <- flowering [mixedorder(rownames(flowering)),]
  if(any(rownames(flowering) != species.names)){
    stop("The rownames of flowering are not the same as those in the species.names vector.")}
  if(length(species.names)!= dim(flowering)[1]){
    stop("For the matrix flowering, only one observation per species is allowed.")}
  if(dim(flowering)[2] != 2){
    stop("The matrix flowering, must have 2 columns (Start and End of the flowering period).")}
  if(any(is.na(flowering))){
    stop("NA values are not allowed in the flowering matrix.")}
  n.traits <- n.traits + 1
  names.traits <- c(names.traits, "flowering")
  trait.type <- c(trait.type, "flowering")
        } #End checking conditions for flowering
n.gaussian <- length(which(trait.type=="gaussian"))
n.kernel <- length(which(trait.type=="kernel"))
n.mean <- length(which(trait.type=="mean"))
n.multiple <- length(which(trait.type=="multiple"))
n.single <- length(which(trait.type=="single"))
n.flowering <- length(which(trait.type=="flowering"))
#create a list to store the results
results <- list()
#Calculate the species pairwise distances for each trait
for (i in 1:n.traits){   #For each trait in n.traits
  trait.name.aux <- names.traits [i]
  print(paste("Trait ", i," / ",n.traits, sep=""))
#1-for gaussian.data
  if (trait.type [i] == "gaussian"){# start of trait.type == "gaussian"
    gaussian.dist<-function(sp1,sp2) {
        range_sp<-range(c(sp1[1]+5*sp1[2],sp1[1]-5*sp1[2],sp2[1]+5*sp2[2],sp2[1]-5*sp2[2]))
        length.grad<-seq(range_sp[1],range_sp[2],length.out=128)
        pdf_sp1<-dnorm(length.grad,sp1[1],sp1[2])
        pdf_sp2<-dnorm(length.grad,sp2[1],sp2[2])
        overlapped<-apply(cbind(pdf_sp1,pdf_sp2), 1, min)
        overlapped.curve<-splinefun(length.grad,overlapped)
        similarity<-integrate(overlapped.curve,range_sp[1],range_sp[2],subdivisions=128)$value
        gaussian.results<-1-similarity
        return(gaussian.results)
        if(gaussian.results>1){gaussian.results<-1}
        if(gaussian.results<0){gaussian.results<-0}
        }
    mean.trait <- gaussian.data[ , which(colnames(gaussian.data)==trait.name.aux)[1]]
    sd.trait <- gaussian.data[ ,1 + which(colnames(gaussian.data)==trait.name.aux)[1]]
    if((any(is.na(sd.trait)) & missing) | (any(sd.trait==0) & missing)){
      index.missing<-which(is.na(sd.trait) | (sd.trait==0))
      sd.trait[index.missing]<-mean(sd.trait,na.rm=TRUE)
    }
    distances.aux<-numeric(length(species.names)^2)  ##create a matrix to store the results
    p<-1
    for (j in 1:length(species.names)){ #for each species j
      sp1<-c(mean.trait[j], sd.trait[j])
      for (k in 1:length(species.names)){ #for each species k
        sp2<-c(mean.trait[k], sd.trait[k])
        if (k > j){
            distances.aux[p]<-gaussian.dist(sp1,sp2)
            } #end k<j
        p<-p+1
        } #end species k
    } #end species j
    matrix.dist<-matrix(distances.aux,nrow=length(species.names),ncol=length(species.names),dimnames=list(species.names,species.names))
    results[[i]] <- as.dist(matrix.dist)
    names(results)[i] <- trait.name.aux
  } #end of trait.type[i] == "gaussian"
#2-for kernel.data
  if (trait.type [i] == "kernel"){  # start of trait.type == "kernel"
    kernel.dist<-function(sp1,sp2,type.sp1,type.sp2)  {
        if (type.sp1=="kern" & type.sp2=="kern"){ #If there are enough observations of both species
        n.obs<-c(length(sp1),length(sp2))
        sdevs<-c(sd(sp1,na.rm=TRUE),sd(sp2,na.rm=TRUE))
        means<-c(mean(sp1,na.rm=TRUE),mean(sp2,na.rm=TRUE))
        h<-1.06*sdevs*(n.obs^-0.2)
        range_sp<-range(c(means[1]+5*sdevs[1],means[1]-5*sdevs[1],means[2]+5*sdevs[2],means[2]-5*sdevs[2]))
        length.grad<-seq(range_sp[1],range_sp[2],length.out=128)
        pdf_x<-density(sp1, bw=h[1],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$x
        pdf_sp1<-density(sp1, bw=h[1],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$y
        pdf_sp2<-density(sp2, bw=h[2],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$y
        }
        if (type.sp1=="gauss" & type.sp2=="gauss"){
          range_sp<-range(c(sp1[1]+5*sp1[2],sp1[1]-5*sp1[2],sp2[1]+5*sp2[2],sp2[1]-5*sp2[2]))
          length.grad<-seq(range_sp[1],range_sp[2],length.out=128)
          pdf_sp1<-dnorm(length.grad,sp1[1],sp1[2])
          pdf_sp2<-dnorm(length.grad,sp2[1],sp2[2])
          }
        if (type.sp1=="gauss" & type.sp2=="kern"){
          n.obs<-length(sp2)
          sdevs<-c(sp1[2],sd(sp2,na.rm=TRUE))
          means<-c(sp1[1],mean(sp2,na.rm=TRUE))
          h<-1.06*sdevs*(n.obs^-0.2)
          range_sp<-range(c(means[1]+5*sdevs[1],means[1]-5*sdevs[1],means[2]+5*sdevs[2],means[2]-5*sdevs[2]))
          length.grad<-seq(range_sp[1],range_sp[2],length.out=128)
          pdf_x<-density(sp2, bw=h[1],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$x
          pdf_sp1<-dnorm(length.grad,sp1[1],sp1[2])
          pdf_sp2<-density(sp2, bw=h[2],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$y
          }
        if (type.sp1=="kern" & type.sp2=="gauss"){
          n.obs<-length(sp1)
          sdevs<-c(sd(sp1,na.rm=TRUE),sp2[2])
          means<-c(mean(sp1,na.rm=TRUE),sp2[1])
          h<-1.06*sdevs*(n.obs^-0.2)
          range_sp<-range(c(means[1]+5*sdevs[1],means[1]-5*sdevs[1],means[2]+5*sdevs[2],means[2]-5*sdevs[2]))
          length.grad<-seq(range_sp[1],range_sp[2],length.out=128)
          pdf_x<-density(sp1, bw=h[1],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$x
          pdf_sp1<-density(sp1, bw=h[2],kernel="gaussian",from=range_sp[1],to=range_sp[2],n=128,na.rm=TRUE)$y
          pdf_sp2<-dnorm(length.grad,sp2[1],sp2[2])
          }
        overlapped<-apply(cbind(pdf_sp1,pdf_sp2), 1, min)
        overlapped.curve<-splinefun(x=length.grad, y=overlapped)
        similarity<-integrate(overlapped.curve,range_sp[1],range_sp[2],subdivisions=128)$value
        kernel.results<-1-similarity
        if(kernel.results>1){kernel.results<-1}
        if(kernel.results<0){kernel.results<-0}
        return(kernel.results)
        }
    distances.aux<-numeric(length(species.names)^2)  ##create a matrix to store the results
    p<-1
    if(missing){community.sd<-mean(tapply(kernel.data[,(i - n.gaussian)],
                                rownames(kernel.data),sd,na.rm=T))}
    for (j in 1:length(species.names)){ #for each species j
      index.sp1<-which(rownames(kernel.data) == species.names[j])
      sp1<-kernel.data [index.sp1 , (i - n.gaussian)]
      for (k in 1:length(species.names)){ #for each species k
        index.sp2<-which(rownames(kernel.data) == species.names[k])
        sp2<-kernel.data [index.sp2 , (i - n.gaussian)]
        if((length(sp1)-sum(is.na(sp1))) < min.observations & missing==FALSE)
        {stop(paste("When type is kernel and missing=FALSE, the minimun number of observations is",
                  min.observations,"\n Only", length(sp1)-sum(is.na(sp1)),"observations for species -",
                  species.names[j],"-, trait -", names.traits[i],"-"))}
        if((length(sp1)-sum(is.na(sp1))) < min.observations & missing){
          type.sp1<-"gauss"
          sp1<-c(mean(sp1,na.rm=T),community.sd)
          } else { type.sp1<-"kern"}
        if((length(sp2)-sum(is.na(sp2))) < min.observations & missing==FALSE)
        {stop(paste("When type is kernel and missing=FALSE, the minimun number of observations is",
                  min.observations,"\n Only", length(sp2)-sum(is.na(sp2)),"observations for species -",
                  species.names[k],"-, trait -", names.traits[i],"-"))}
        if((length(sp2)-sum(is.na(sp2))) < min.observations & missing){
          type.sp2<-"gauss"
          sp2<-c(mean(sp2,na.rm=T),community.sd)
          } else { type.sp2<-"kern"}
        if (k > j){
            distances.aux[p]<-kernel.dist(sp1,sp2,type.sp1,type.sp2)
        } #end k<j
        p<-p+1
      } #end species k
    } #end species j
  matrix.dist<-matrix(distances.aux,nrow=length(species.names),ncol=length(species.names),dimnames=list(species.names,species.names))
  results[[i]] <- as.dist(matrix.dist)
  names(results)[i] <- trait.name.aux
  } #end of trait.type == "kernel"
#3-for mean.data
  if (trait.type [i] == "mean"){# start of trait.type == "mean"
    mean.dist<-function(sp1,sp2) {
        distance<-abs(sp1-sp2)
        mean.results<-list()
        mean.results$distance<-distance
        return(mean.results)
        }
    distances.aux<-matrix(0,nrow=length(species.names),ncol=length(species.names),
                      dimnames=list(species.names,species.names))  ##create a matrix to store the results
    length.gradient<-abs(max(mean.data[, i-(n.gaussian + n.kernel)])-
                     min(mean.data[, i-(n.gaussian + n.kernel)]))
    for (j in 1:length(species.names)){ #for each species j
      index.sp1<-which(rownames(mean.data) == species.names[j])
      sp1<-mean.data[index.sp1, i-(n.gaussian + n.kernel)]
      for (k in 1:length(species.names)){ #for each species k
        index.sp2<-which(rownames(mean.data) == species.names[k])
        sp2<-mean.data[index.sp2, i-(n.gaussian + n.kernel)]
        if (k < j){
            aux<-mean.dist(sp1,sp2)
            distances.aux[j,k]<-aux$distance/length.gradient
            } #end k<j
        } #end species k
    } #end species j
  results[[i]] <- as.dist(distances.aux)
  names(results)[i] <- trait.name.aux
  } #end of trait.type == "mean"
#4-for multiple.category
  if (trait.type [i] == "multiple"){  # start of trait.type == "multiple"
    multiple.dist<-function(sp1,sp2)  {
        similarity<-sum(pmin(sp1,sp2))
        distance<-1-similarity
        multiple.results<-list()
        multiple.results$distance<-distance
        return(multiple.results)
        }
    distances.aux<-matrix(0,nrow=length(species.names),ncol=length(species.names),
                      dimnames=list(species.names,species.names))  ##create a matrix to store the results
    selected.matrix <- multiple.category[[i- (n.gaussian + n.kernel + n.mean)]]
    for (j in 1:length(species.names)){ #for each species j
      index.sp1<-which(rownames(selected.matrix) == species.names[j])
      sp1<-selected.matrix[index.sp1 ,]
      for (k in 1:length(species.names)){ #for each species k
        index.sp2<-which(rownames(selected.matrix) == species.names[k])
        sp2<-selected.matrix[index.sp2 ,]
        if (k < j){
            aux<-multiple.dist(sp1,sp2)
            distances.aux[j,k]<-aux$distance
            } #end k<j
        } #end species k
    } #end species j
  results[[i]] <- as.dist(distances.aux)
  names(results)[i] <- trait.name.aux
  } #end of trait.type == "multiple"
#4-for single.category
  if (trait.type [i] == "single"){  # start of trait.type == "single"
    single.dist<-function(sp1,sp2)  {
        distance<-ifelse (sp1==sp2, 0, 1)
        single.results<-list()
        single.results$distance<-distance
        return(single.results)
        }
    distances.aux<-matrix(0,nrow=length(species.names),ncol=length(species.names),
                      dimnames=list(species.names,species.names))  ##create a matrix to store the results
    for (j in 1:length(species.names)){ #for each species j
      index.sp1<-which(rownames(single.category) == species.names[j])
      sp1<-single.category[index.sp1, i-(n.gaussian + n.kernel + n.mean + n.multiple)]
      for (k in 1:length(species.names)){ #for each species k
        index.sp2<-which(rownames(single.category) == species.names[k])
        sp2<-single.category[index.sp2, i-(n.gaussian + n.kernel + n.mean + n.multiple)]
        if (k < j){
            aux<-single.dist(sp1,sp2)
            distances.aux[j,k]<-aux$distance
            } #end k<j
        } #end species k
    } #end species j
  results[[i]] <- as.dist(distances.aux)
  names(results)[i] <- trait.name.aux
  } #end of trait.type == "single"
#4-for flowering
  if (trait.type [i] == "flowering"){  # start of trait.type == "flowering"
    flowering.dist<-function(sp1,sp2)  {
        if (sp1[1] <= sp1[2]){#si la floracion empieza antes (el el ano natural) y acaba despues eg: febrero-mayo
          flower1 <- sp1[1]:sp1[2]}
          else{flower1 <- c(1:sp1[2],sp1[1]:365)}
        if (sp2[1] <= sp2[2]){#si la floracion empieza antes (el el ano natural) y acaba despues eg: febrero-mayo
          flower2 <- sp2[1]:sp2[2]}
          else{flower2 <- c(1:sp2[2],sp2[1]:365)}
        similarity<-sum(flower1 %in% flower2)/(min (length(flower1),length(flower2)))
        distance<-1-similarity
        flowering.results<-list()
        flowering.results$distance<-distance
        return(flowering.results)
        }
    distances.aux<-matrix(0,nrow=length(species.names),ncol=length(species.names),
                      dimnames=list(species.names,species.names))  ##create a matrix to store the results
    for (j in 1:length(species.names)){ #for each species j
      index.sp1<-which(rownames(flowering) == species.names[j])
      sp1<-flowering[index.sp1, ]
      for (k in 1:length(species.names)){ #for each species k
        index.sp2<-which(rownames(flowering) == species.names[k])
        sp2<-flowering[index.sp2, ]
        if (k < j){
            aux<-flowering.dist(sp1,sp2)
            distances.aux[j,k]<-aux$distance
            } #end k<j
        } #end species k
    } #end species j
  results[[i]] <- as.dist(distances.aux)
  names(results)[i] <- trait.name.aux
  } #end of trait.type == "flowering"
}#end of trait[i]
#Calculte the distance for all the traits:
sumdist<-mindist<-maxdist<-(results[[1]])
sumdist2<-(results[[1]]^2)
if(n.traits>1){

  for(i in 2:n.traits){
    sumdist<-sumdist+(results[[i]])
    sumdist2<-sumdist2+(results[[i]]^2)
    mindist<-pmin(mindist,results[[i]])
    maxdist<-pmax(maxdist,results[[i]])
  }
  results$gower<-sumdist/n.traits
  results$euclid<-sqrt(sumdist2)/sqrt(n.traits)     
  results$mindist<-mindist
  results$maxdist<-maxdist
  }
out <- results
return(out)
} #end of trova



####################################################################################################################
####################################################################################################################
### Creation of example data:
###Gaussian data (mean and sd)
sp.names<- c( "spec_1" ,"spec_3","spec_2","spec_4")

spec_1<-c(22, 2, 10, 2)
spec_2<-c(8, 3, 8, 3)
spec_3<-c(12, 5, 12, 5)
spec_4<-c(20, 7, 5, 2)
  
datos.gaussian<-data.frame(rbind(spec_4,spec_3,spec_2,spec_1))
colnames(datos.gaussian)<-c("trait1.gaussian","sd.trait1.gaussian",
                            "trait2.gaussian","sd.trait2.gaussian")
                              
###Kernel data (a number of observations for each species)
observations<-60
trait1.kernel<-c(rnorm(observations/4, spec_1[1],spec_1[2]),
                 rnorm(observations/4, spec_2[1],spec_2[2]),
                 rnorm(observations/4, spec_3[1],spec_3[2]),
                 rnorm(observations/4, spec_4[1],spec_4[2]))
                 
trait2.kernel<-c(rnorm(observations/4, spec_1[3],spec_1[4]),
                 rnorm(observations/4, spec_2[3],spec_2[4]),
                 rnorm(observations/4, spec_3[3],spec_3[4]),
                 rnorm(observations/4, spec_4[3],spec_4[4]))
                 
datos.kernel<-as.matrix(cbind(trait1.kernel,trait2.kernel))
species.kernel<-c(rep("spec_1",observations/4),rep("spec_2",observations/4),
                          rep("spec_3",observations/4),rep("spec_4",observations/4))
rownames(datos.kernel)<-species.kernel
colnames(datos.kernel)<-c("trait.kernel1","trait.kernel2")
    
###Mean data (quantitative traits for which we only know the mean value for each species)
spec_1<-c(22, 10)
spec_2<-c(8, 8)
spec_3<-c(12, 12)
spec_4<-c(20, 5)
datos.mean<-rbind(spec_4,spec_3,spec_2,spec_1)
colnames(datos.mean)<-c("trait1.mean","trait2.mean")
  
###Multiple category data (proportion of individuals belonging to each value of the traits)
#Trait 1  (3 categories)
T1.spec_1<-c(0.5, 0.25, 0.25)
T1.spec_2<-c(0.5, 0.5, 0)
T1.spec_3<-c(0.5, 0.25, 0.25)
T1.spec_4<-c(0, 0, 1)
#Trait 2  (4 categories)
T2.spec_1<-c(0, 0, 0.5, 0.5)
T2.spec_2<-c(0.5, 0.5, 0, 0)
T2.spec_3<-c(0.25, 0.25, 0.25, 0.25)
T2.spec_4<-c(1, 0, 0, 0)
    
datos.multiple<-list()
T1.multiple<-rbind(T1.spec_3,T1.spec_2,T1.spec_1,T1.spec_4)
rownames(T1.multiple)<-c("spec_3","spec_2","spec_1","spec_4")
datos.multiple[[1]]<-T1.multiple
names(datos.multiple)[1]<-"trait1.multiple"
  
T2.multiple<-rbind(T2.spec_2,T2.spec_1,T2.spec_4,T2.spec_3)
rownames(T2.multiple)<-c("spec_2","spec_1","spec_4","spec_3")
datos.multiple[[2]]<-T2.multiple
names(datos.multiple)[2]<-"trait2.multiple"
###Single category data (all individuals of each species have to the same value for the trait)
spec_1<-c("annual", "terophyte")
spec_2<-c("perennial", "criptophyte")
spec_3<-c("annual", "terophyte")
spec_4<-c("perennial", "hemicriptophyte")
    
datos.single<-rbind(spec_4,spec_2,spec_3,spec_1)
colnames(datos.single)<-c("life_span","life_form")
###Flowering data (the start and end in Julian days of the flowering period)
spec_1<-c(50, 200)
spec_2<-c(300, 100)  
spec_3<-c(1, 49)
spec_4<-c(1, 150)
datos.flowering<-rbind(spec_1,spec_3,spec_4,spec_2)