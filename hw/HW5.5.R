blci <- function(cd, id, total){
  cdc <- setdiff(total, cd)
  idc <- setdiff(total, id)
  length(intersect(cd,id))/length(cd)+length(intersect(cdc,idc))/length(cdc)-1
}

cpp <- function(a, b){
  c <- matrix(0,nrow=10,ncol=10)
  for(i in seq(1:length(a))) c[a[i], b[i]] <- c[a[i], b[i]] + 1
  sum(apply(c, 1, max))/length(a)
}

load("/Users/andyisman/Documents/BioInfo/lymphoma/cpp.RData") 
load("/Users/andyisman/Documents/BioInfo/lymphoma/blci.RData")

cat("cpp: ", cpp(ans.clust, imp.clust), "\n")
cat("blci: ", blci(ans.siggene.list, imp.siggene.list, union(ans.nonsiggene.list, ans.siggene.list)))