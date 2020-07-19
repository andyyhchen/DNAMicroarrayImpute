
similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}


sknn <- function(xmiss, K) {
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xincomplete <- xmiss[miss.row, ]
  miss.inc <- is.na(xincomplete)
  miss.origin <- order(rowSums(miss.inc))
  xincomplete <- xincomplete[miss.origin, ]
  xcomplete <- xmiss[-miss.row, ]
  xtmp <- matrix(nc=ncol(xincomplete))
  impute <- function(row) {  
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    gene <- row[-row.exp]
    cand_x <- xcomplete[, -row.exp]
    d <- sqrt(rowSums((cand_x - matrix(gene, nc = length(gene), nr = nrow(cand_x), byrow=T))^2))
    
    id.idx <- order(d)[1:K]
    id.sel <- d[id.idx]
    
    const <- sum(1/id.sel)
    w <- 1/(const*id.sel)
    
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    xcomplete <<- rbind(xcomplete, row)
    
    return (row)
  } 
  
  xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, impute))
    
  return (xmiss)
}
sknn_1 <- function(xmiss, K=15, sim.method="EuDist"){
  
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

SKNNimpute<-function(x, k=10, sim.method="EuDist"){
  rowNum<-nrow(x)
  colNum<-ncol(x)
  missIdx<-is.na(x)
  miss.RowIdx<-which(rowSums(missIdx)!=0)
  x.completeRows<-x[-miss.RowIdx,]
  x.missingRows<-x[miss.RowIdx,]
  miss.list<-order(rowSums(is.na(x.missingRows)))
  for (i in seq(nrow(x.missingRows))) {
    target<-x.missingRows[miss.list[i],]
    dist.list<-similarityCal(target, x.completeRows, sim.method)
    neighborsIdx<-order(abs(dist.list),decreasing=T)[1:k]
    missColIdx<-which(is.na(target))
    estimation<-sapply(missColIdx, function(j){
      weight<-dist.list[neighborsIdx]
      weightedAvg<-weight %*% x.completeRows[neighborsIdx,j] / sum(weight)
      return(weightedAvg)
    })
    target[missColIdx]<-estimation
    x.missingRows[miss.list[i],]<-target
    x.completeRows<-rbind(x.completeRows, target)
  }
  x[miss.RowIdx,]<-x.missingRows
  return(x)
}
nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_5_1.txt", sep=""), sep="\t")
ptc <- proc.time()
cat(nrmse(missing, sknn_1(missing, 15), answer),"\n")
