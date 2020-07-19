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

LLS <- function(xmiss, K=15, sim.method="EuDist"){
  
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
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
    row[row.miss] <- t(x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}

LS <- function(xmiss, K=15, sim.method="cor"){
  
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
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
    
    sim <- similarityCal(gene,cand_x, method=sim.method)
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

knn <- function(xmiss, K=15, sim.method="EuDist"){
  
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
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
    const <- sum(id.sel)
    w <- 1/const*id.sel
    w <- matrix(w, nc=length(w), nr=1)
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    return (row)
  }))
  
  return (xmiss)
}

sknn <- function(xmiss, K=15, sim.method="EuDist"){
  
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xincomplete <- xmiss[miss.row, ]
  miss.inc <- is.na(xincomplete)
  miss.origin <- order(rowSums(miss.inc))
  xincomplete <- xincomplete[order(rowSums(miss.inc)), ]
  xcomplete <- xmiss[-miss.row, ]
  xtmp <- matrix(nc=ncol(xincomplete))
  xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, function(row){
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
    id.idx <- order(d, decreasing=T)[1:K]
    id.sel <- d[id.idx]
    const <- sum(id.sel)
    w <- 1/const*id.sel
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    xcomplete <<- rbind(xcomplete, row)
    return (row)
  }))
  
  return (xmiss)
}

iknn <- function(xmiss, K=15, sim.method="EuDist", Niter=3){
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
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
      const <- 1/sum(id.sel)
      w <- const * id.sel
      row[row.exp] <- w %*% xcomplete[id.idx, row.exp, drop=F]
      return (row)
    }))
  }
  return(xcomplete) 
}