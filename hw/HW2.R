rm(list = ls())
ptm<-proc.time()
## library ##
library(pdist)
#library(ForImp)
library(imputation)
library(reshape2)
library(ggplot2)
library(grid)
library(scales)
library(samr)

processFile <- function(filename){ as.matrix(read.table(filename,header=F,sep='\t')) }
nrmse <- function(a, b){ sqrt(mean((b-a)**2)/var(a)) }

opt <- list(
  M = c(1,5,10,15,20),
  Method = c("zero", "rowaverage", "iknn_cor", "sknn_cor", "knn_cor", "LSimpute"),
  outFile = format(Sys.time(), "%Y_%m_%d_%H_%M"),
  K = 15	
)

zero <- function(xmiss){
  xmiss[is.na(xmiss)] <- 0
  xmiss
}

rowaverage <- function(xmiss){
  mis <- is.na(xmiss)
  rowM <- rowMeans(xmiss, na.rm=T)
  rowM <- matrix(rowM, nc=ncol(xmiss), nr=nrow(xmiss), byrow=F)
  
  xmiss[mis] <- rowM[mis]
  xmiss
}

LLS <- function(xmiss, K=15, sim.method="euclidean"){
  
  similarityCal<-function(vec, mat, method="euclidean"){
    methods<-c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    row[row.miss] <- ans<-t(x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}

LS <- function(xmiss, K=15, sim.method="pearson"){
  
  similarityCal<-function(vec, mat, method="euclidean"){
    methods<-c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  
  xcomplete <- xmiss[-miss.row, ]
  xincomplete <- xmiss[miss.row, ]
  
  impute <- function(row) {
    
    row.exp <- which(is.na(row))
    gene <- row[-row.exp]
    cand_x <- xcomplete[, -row.exp, drop=F]
    
    sim <- similarityCal(gene,xcomplete[,-row.exp], method=sim.method)
    row.idx <- order(sim, decreasing=T)[1:K]
    row.r <- sim[row.idx]
    row.cand <- cand_x[row.idx, , drop=F]
    lg <- apply(row.cand, 1, function(x){lm(gene~x)$coefficients})
    row.impcand <- xcomplete[row.idx, row.exp, drop=F]
    
    y <- matrix(0, nc=ncol(row.impcand), nr=nrow(row.impcand))
    w <- (row.r**2/(1-row.r**2+0.000001))**2
    sw <- sum(w)
    w <- w/sw
    
    for(i in 1:nrow(row.impcand)) {
      y[i, ] <- lg[2, i] * row.impcand[i, ] + lg[1, i]
    }
    row[row.exp]<-apply(y, 2, function(x){sum(w*x)})
    return(row)
  }
  
  xmiss[miss.row, ] <- t(apply(xincomplete, 1, impute))
  
  return (xmiss)
}

knn <- function(xmiss, K=15, sim.method="euclidean"){
  
  similarityCal<-function(vec, mat, method="euclidean"){
    methods<-c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xcomplete <- xmiss[-miss.row, ]
  xincomplete <- xmiss[miss.row, ]
  
  xmiss[miss.row, ] <- t(apply(xincomplete, 1, function(row){
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
    id.idx <- order(d, decreasing=T)[1:K]
    id.sel <- d[id.idx]
    const <- sum(1/id.sel)
    w <- 1/(const*id.sel)
    w <- matrix(w, nc=length(w), nr=1)
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    return (row)
  }))
  
  return (xmiss)
}

sknn <- function(xmiss, K=15, sim.method="euclidean"){
  
  similarityCal<-function(vec, mat, method="euclidean"){
    methods<-c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xincomplete <- xmiss[miss.row, ]
  miss.inc <- is.na(xincomplete)
  miss.origin <- order(rowSums(miss.inc))
  xincomplete <- xincomplete[order(rowSums(miss.inc)), ]
  xcomplete <<- xmiss[-miss.row, ]
  xtmp <- matrix(nc=ncol(xincomplete))
  xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, function(row){
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F])
    id.idx <- order(d, decreasing=T)[1:K]
    id.sel <- d[id.idx]
    const <- sum(1/id.sel)
    w <- 1/(const*id.sel)
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    xcomplete <<- rbind(xcomplete, row)
    return (row)
  }))
  
  return (xmiss)
}

iknn <- function(xmiss, K=15, sim.method="euclidean", Niter=3){
  
  similarityCal<-function(vec, mat, method="euclidean"){
    methods<-c("euclidean","pearson","cosine")
    switch(match.arg(method,methods),
           euclidean=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           pearson=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           cosine=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  rowMeanSubstitution <- function(xmiss) {
    rowM <- rowMeans(xmiss, na.rm = T)
    rowM <- matrix(rowM, nc = ncol(xmiss), nr = nrow(xmiss), byrow=F)
    
    E <- xmiss
    E[is.na(xmiss)] <- rowM[is.na(xmiss)]
    
    return(E)
  }
  
  xcomplete <- rowMeanSubstitution(xmiss);
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, ])))
  
  for(h in 1:Niter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- is.na(row)
      row.exp <- which(row.miss)
      d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
      id.idx <- order(d, decreasing=T)[2:(K+1)]
      id.sel <- d[id.idx]
      const <- sum(1/id.sel)
      w <- 1/(const*id.sel)
      row[row.exp] <- w %*% xcomplete[id.idx, row.exp, drop=F]
      return (row)
    }))
  }
  return(xcomplete) 
}


# algorithm list
methodList <- list(
  zero = function(m){
    zero(m)
  },
  rowaverage = function(m){
   rowaverage(m)
  },
  knn_eu = function(m){
    knn(m, opt$K, "euclidean")
  },
  knn_cor = function(m){
    knn(m, opt$K, "pearson")
  },
  knn_angle = function(m){
    KNN(m ,opt$K, "consine")
  },
  sknn_cor = function(m){
    sknn(m, opt$K, "pearson")
  },
  iknn_cor = function(m){
    iknn(m, opt$K, "pearson", 2)
  },
  LSimpute = function(m){
    LS(m, opt$K, "pearson")
  },
  LLSimpute = function(m){
    LLS(m, opt$K, "euclidean")
  }
)
# for cpp
cpp <- function(a, b){
  c <- matrix(0,nrow=10,ncol=10)
  for(i in seq(1:length(a))) c[a[i], b[i]] <- c[a[i], b[i]] + 1
  sum(apply(c, 1, max))/length(a)
}

# for BLCI
findsig <- function(x, genenames){
  x.d<-list(x=x,eigengene.number=10,geneid=genenames,genenames=genenames)
  x.samr<-samr(x.d,resp.type="Pattern discovery",nperms=100)
  x.delta<-samr.compute.delta.table(x.samr)
  x.siggenes.table<-samr.compute.siggenes.table(x.samr,del=0,x.d,x.delta,all.genes=T)
  x.totalgenes<-rbind(x.siggenes.table$genes.lo,x.siggenes.table$genes.up)
  x.totalgenes[as.numeric(x.totalgenes[,7])<10,2]
}
blci <- function(cd, id, total){
  cdc <- setdiff(total, cd)
  idc <- setdiff(total, id)
  length(intersect(cd,id))/length(cd)+length(intersect(cdc,idc))/length(cdc)-1
}
#plotting
myplot <- function(data, name){
  melt <- melt(data, varnames=c('M','Method'), value.name='NRMSE')
  p1 <- ggplot(melt, aes(x=M,y=NRMSE,colour=Method)) + 
    geom_line(aes(group=Method),size=1.5) + 
    geom_point(aes(shape=Method),size=7) +
    xlab("Missing Rate (%)") +
    ylab(name) +
    theme(axis.text=element_text(size=26),
          axis.title=element_text(size=30),
          legend.title=element_text(size=28),
          legend.text=element_text(size=26),
          legend.key.size=unit(1,"cm"),
          legend.background=element_rect(fill=alpha("gray",0.5)),
          legend.margin=unit(2.5,"cm"),
          legend.position=c(.9,.1))
  
  png(filename=paste("/Users/andyisman/Documents/BioInfo/lymphoma/",opt$outFile,'_',name,'.png',sep=""), width=700, height=700, bg='white')
  print(p1)
  dev.off()
  
}

# main
# create empty table
resTable.nrmse <- matrix(0, nrow = length(opt$M), ncol = length(opt$Method), byrow = F, dimname=list(opt$M, opt$Method))
resTable.blci <- resTable.nrmse
resTable.cpp <- resTable.nrmse
resTable.time <- resTable.nrmse
# read ans and calculate ans kmeans and siggene
ans <- processFile('/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt')
totalGene <- as.character(1:nrow(ans))
ans.cluster <- kmeans(ans, 10)$cluster
ans.sig <- findsig(ans, totalGene)
# read incomplete and imputing
for ( mis in opt$M ) { 
  input <- processFile(paste('/Users/andyisman/Documents/BioInfo/lymphoma/m_', mis, '_1.txt', sep=''))
  ans_v <- ans[is.na(input)]
  for( met in opt$Method ) {
    cat(met,"\n")
    ptm <- proc.time()
    output <- methodList[[met]](input)
    resTable.time[toString(mis), met] <- (proc.time() - ptm)[3]
    output.cluster <- kmeans(output, 10)$cluster
    resTable.cpp[toString(mis), met] <- cpp(ans.cluster, output.cluster)
    output.sig <- findsig(output, totalGene)
    resTable.blci[toString(mis), met] <- blci(ans.sig, output.sig, totalGene)
    output_v <- output[is.na(input)]
    resTable.nrmse[toString(mis), met] <- nrmse(ans_v, output_v)
  }
}


# plot region
stopCluster(c1)
myplot(resTable.nrmse,"NRMSE")
#myplot(resTable.blci,"BLCI")
#myplot(resTable.cpp,"CPP")
myplot(resTable.time, "TIME")
running.time<-proc.time()-ptm
print(running.time)

# vi:nu:nowrap:ts=4:st=4