VarScoreOMI<-function (score, df)
{
    nvar <- ncol(df)
    mat<-matrix(0,ncol=2,nrow=nvar)
    sum.col <- apply(df, 2, sum)
    df <- df[, sum.col > 0]
    nvar <- ncol(df)
    sum.col <- apply(df, 2, sum)
    df <- sweep(df, 2, sum.col, "/")
    for (i in 1:nvar) {
        w <- df[, i]
        x.moy <- sum(w * score)
        x.et <- sqrt(sum(w * (score - x.moy)^2))
        mat[i,1]<-x.moy
        mat[i,2]<-(x.et)*(x.et)
    }
    return (mat)
}