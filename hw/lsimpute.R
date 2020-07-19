
LsImpute <- function(xmiss, K=15, sim.method="pearson"){
  
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


nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_1_1.txt", sep=""), sep="\t")

cat(nrmse(missing, LsImpute(missing, 15), answer))