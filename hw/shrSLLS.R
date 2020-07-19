require(MASS)

shrSLLS <- function(xmiss, K=15, sim.method="EuDist"){
  
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
  miss.rate <- rowSums(is.na(x.incomplete))
  miss.rate.order <- order(miss.rate)
  therhold <- sum(miss.rate)/length(miss.rate)
  x.incomplete <- x.incomplete[miss.rate.order, ]
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    x.tmp <- ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    x.shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
    row[row.miss] <- t(x.complete[sim.id, row.miss, drop=FALSE]) %*% x.shrink
    if(length(row.miss) <= therhold){
      x.complete <<- rbind(x.complete, row)
    }
    return(row)
  }))
  
  xmiss[miss.row[miss.rate.order], ] <- x.imputed
  
  return(xmiss)
}

nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_15_1.txt", sep=""), sep="\t")

cat(nrmse(missing, shrSLLS(missing), answer),"\n", proc.time() - ptm, "\n")