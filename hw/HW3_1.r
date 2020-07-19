rm(list=ls())
require(grid)
library(ggplot2)
library(reshape2)

#output 
opt <- list(
  rate=c(1, 5, 10, 15, 20), # missing rate
  k=c(5, 10, 12, 14, 16, 18, 20, 30, 50, 80, 100), #k value
  outfmt="hw3_diff_kNN.png" #output file format
)

#knn_impute

simplyfied_knn_impute <- function(matrix, k) {
  
  #classfied
  nomissing <- matrix[which(rowSums(is.na(matrix)) == 0), ]
  missing <- matrix[which(rowSums(is.na(matrix)) != 0), ]
  
  #euclid distance
  distance <- function(x, y) {
    common <- !is.na(x)&!is.na(y)
    sqrt(sum((x[common]-y[common])^2))
  }
  
  #impute
  for(i in 1:nrow(missing)){
    dist <- NULL
    for(j in 1:nrow(nomissing)){
      dist <- c(dist, distance(missing[i,], nomissing[j,]))
    }
    
    ksd <- order(dist)[1:k] #k_smallest distance
    
    for(id in which(is.na(missing[i,]))) {
      missing[i, id] <- sum(1/dist[ksd]*nomissing[ksd, id])/sum(1/dist[ksd])
    }
  }
  
  matrix[which(rowSums(is.na(matrix)) != 0), ] <- missing
  
  matrix
}

#nrmse

nrmse <- function(missing, guess, answer) { 
  i = is.na(missing)
  sqrt(mean((guess[i]-answer[i])^2)/var(answer[i])) 
}

#---

result = matrix(nr=length(opt$rate), nc=length(opt$k), dimnames=list(opt$rate, opt$k))
answer = as.matrix(read.table("ans.txt", sep="\t"))

for(i in opt$rate) {
  missing <- as.matrix(read.table(paste("m_", i, "_1.txt", sep="")), sep="\t")
  for(j in opt$k){
    #result[toString(i), toString(j)] <- nrmse(missing, simplyfied_knn_impute(missing, j), answer)
    result[toString(i), toString(j)] <- nrmse(missing, kNNImpute(missing, j)$x, answer)
    cat(i,"%_k=",j,"\n")
  }
}
rm(missing, i, j)

result

#--

pic <- ggplot(melt(result, varnames=c('RATE', 'K'), value.name='NRMSE'), aes(x=RATE, y=NRMSE, colour=factor(K))) + 
  geom_line(aes(group=factor(K)),size=1.5) + 
  xlab("Missing Rate (%)") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.title=element_text(size=14),
        legend.text=element_text(size=10),
        legend.key.size=unit(1,"cm"),
        legend.background=element_rect(colour="gray"),
        legend.margin=unit(2.5,"cm"),
        legend.position=c(1,.40))
pic

png(filename=opt$outfmt, width=1500, height=700, bg='white')
print(pic)
dev.off()

