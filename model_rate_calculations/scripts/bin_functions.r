bin_df = function(df,col=1,nbins=20,method="equal") {
    p25 = function(x) quantile(x,prob=0.25,na.rm=T)
    p75 = function(x) quantile(x,prob=0.75,na.rm=T)
    p5 = function(x) quantile(x,prob=0.05,na.rm=T)
    p95 = function(x) quantile(x,prob=0.95,na.rm=T)

# bin data 
    if(method=="quantile"){
        probs = c(0:nbins)/nbins # probability value of bins
        bin = cut(df[,col],breaks=quantile(df[,col],probs))
    } else if (method=="equal"){
        bin = cut(df[,col],breaks=nbins)
    } else {
        print("error in bin_df, method should be quantile or equal")
        stop
    }

# means
    df.mean = aggregate(df,by=list(bin),FUN=mean)
    names(df.mean) = paste("mean",names(df.mean),sep=".")
    names(df.mean)[1] = paste(names(df)[col],"bin",sep=".")

# medians
    df.median = aggregate(df,by=list(bin),FUN=median)
    names(df.median) = paste("median",names(df.median),sep=".")
    names(df.median)[1] = paste(names(df)[col],"bin",sep=".")

# standard deviations of the means
    df.sd = aggregate(df,by=list(bin),FUN=sd)[,-1]
    names(df.sd) = paste("sd",names(df.sd),sep=".")
    df.tmp = df
    df.tmp$ones = 1
    df.sum = aggregate(df.tmp,by=list(bin),FUN=sum)
    df.sd$npoints = df.sum$ones

# 5% quantile
    df.p5 = aggregate(df,by=list(bin),FUN=p5)[,-1]
    names(df.p5) = paste("p5",names(df.p5),sep=".")

# 25% quantile
    df.p25 = aggregate(df,by=list(bin),FUN=p25)[,-1]
    names(df.p25) = paste("p25",names(df.p25),sep=".")

# 75% quantile
    df.p75 = aggregate(df,by=list(bin),FUN=p75)[,-1]
    names(df.p75) = paste("p75",names(df.p75),sep=".")

# 95% quantile
    df.p95 = aggregate(df,by=list(bin),FUN=p95)[,-1]
    names(df.p95) = paste("p95",names(df.p95),sep=".")

    df.bin= cbind(df.mean,df.median[,-1],df.sd,df.p25,df.p75,df.p5,df.p95)
    df.bin
}

