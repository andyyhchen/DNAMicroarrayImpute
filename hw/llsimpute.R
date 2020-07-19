require(MASS)
LLSimpute <- function(xmiss, K=15, sim.method="euclidean"){
  
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


nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_15_1.txt", sep=""), sep="\t")

cat(nrmse(missing, LLSimpute(missing, 15), answer))