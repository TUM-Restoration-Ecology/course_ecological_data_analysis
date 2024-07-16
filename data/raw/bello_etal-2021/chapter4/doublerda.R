dbrda=function (dudi, df1=NULL,df2=NULL, scannf = TRUE, nf = 2,scaleX=TRUE,centerX=TRUE,scaleZ=TRUE,centerZ=TRUE) 
{
    lm.pcaiv <- function(x, df, weights, use) {
        if (!inherits(df, "data.frame")) 
            stop("data.frame expected")
        reponse.generic <- x
        begin <- "reponse.generic ~ "
        fmla <- as.formula(paste(begin, paste(names(df), collapse = "+")))
        df <- cbind.data.frame(reponse.generic, df)
        lm0 <- lm(fmla, data = df, weights = weights)
        if (use == 0) 
            return(predict(lm0))
        else if (use == 1) 
            return(residuals(lm0))
        else if (use == -1) 
            return(lm0)
        else stop("Non convenient use")
    }
    if (!inherits(dudi, "dudi")) 
        stop("dudi is not a 'dudi' object")
    if(!is.null(df1)){    
        df1 <- data.frame(df1)
        if (!inherits(df1, "data.frame")) 
            stop("df1 is not a 'data.frame'")
        if (nrow(df1) != length(dudi$lw)) 
            stop("Non convenient dimensions")
        weights1 <- dudi$lw
        isfactor1 <- unlist(lapply(as.list(df1), is.factor))
        for (i in 1:ncol(df1)) {
            if (!isfactor1[i]) 
                df1[, i] <- scalewt(df1[, i], weights1,scale=scaleX,center=centerX)
            }
        }
    if(!is.null(df2)){
        df2 <- data.frame(df2)
        if (!inherits(df2, "data.frame")) 
            stop("df2 is not a 'data.frame'")
        if (nrow(df2) != length(dudi$cw)) 
            stop("Non convenient dimensions")
        weights2 <- dudi$cw
        isfactor2 <- unlist(lapply(as.list(df2), is.factor))
        for (i in 1:ncol(df2)) {
            if (!isfactor2[i]) 
                df2[, i] <- scalewt(df2[, i], weights2,scale=scaleZ,center=centerZ)
        }
    }
    tab=dudi$tab
    if(!is.null(df1)){
        tab <- data.frame(apply(dudi$tab, 2, lm.pcaiv, df = df1, use = 0, 
        weights = dudi$lw))}
    if(!is.null(df2)){    
        tab <- data.frame(apply(tab, 1, lm.pcaiv, df = df2, use = 0, 
        weights = dudi$cw))
        tab=as.data.frame(t(tab))
        }
    X <- as.dudi(tab, dudi$cw, dudi$lw, scannf = scannf, nf = nf, 
        call = match.call(), type = "double rda")
    if(!is.null(df1)){X$X <- df1}
    if(!is.null(df2)){X$Z <- df2}
    X$Y <- dudi$tab
    if(!is.null(df1)){
      ## anciennement on predisait les l1 --- A checker
      U <-  apply(X$li, 2, function(x) coefficients(lm.pcaiv(x, df1, weights1, -1)))
      ##print(U)
      ##row.names(U)[-1] <- names(X$X)
      X$e <- data.frame(U)
      fmla <- as.formula(paste("~ ", paste(names(df1), collapse = "+")))
      w <- scalewt(model.matrix(fmla, data = df1), weights1) * weights1
      w <- t(w) %*% as.matrix(X$l1)
      w <- data.frame(w)
      X$corX <- w
    }
    if(!is.null(df2)){
      U <- data.frame(apply(tab, 1, function(x) coefficients(lm.pcaiv(x, df2, weights2, use=-1))))
      X$U<-U
      U <- as.matrix(signif(U,15))%*%diag(dudi$lw)%*%as.matrix(X$l1)
      X$t=U
      fmla <- as.formula(paste("~ ", paste(names(df2), collapse = "+")))
      w <- scalewt(model.matrix(fmla, data = df2), weights2) * weights2
      w <- t(w) %*% as.matrix(X$c1)
      w <- data.frame(w)
      X$corZ <- w
    }
    #row.names(U)[-1] <- names(X$Z)
    #names(U) <- names(X$li)
    return(X)
}


dbcentrage=function(Y)
{
row.mean=apply(Y,1,mean)
col.mean=apply(Y,2,mean)
res=sweep(Y,1,row.mean,"-")
res=sweep(res,2,col.mean,"-")
return(res+mean(Y))
}



#data=scan("data.txt")
#data=matrix(data,33,15,byrow=T)
#Y=data[,1:4]
#X=data[,5:15]
#Z=matrix(c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1),4,3,byrow=T)
## dudiY=dudi.pca(Y,scannf=F,nf=2,scale=F,center=F)
## dudiXYZ1=dbrda(dudiY,X,Z,scannf=F,nf=4)
## newY=dbcentrage(Y)
## dudinewY=dudi.pca(newY,scannf=F,nf=2,scale=F,center=F)
## dudiXYZ2=dbrda(dudi = dudinewY, df1 = X, df2 = Z,scannf=F,nf=4)
## profY=sweep(Y,1,apply(Y,1,sum),"/")
## dudinewY2=dudi.pca(profY,scannf=F,nf=2,scale=F)
## dudiXYZ3=dbrda(dudi = dudinewY2, df1 = X, df2 = Z,scannf=F,nf=4)
## dudinewY3=dudi.pca(dbcentrage(profY),scannf=F,nf=2,scale=F,center=F)
## dudiXYZ4=dbrda(dudi = dudinewY3, df1 = X, df2 = Z,scannf=F,nf=4)
