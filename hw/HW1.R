rawData <- as.matrix(read.table("/Users/andyisman/Downloads/lymphoma.pcl", na.strings="", sep="\t"))

mx <- matrix(as.numeric(rawData[-1,-1]), nr=nrow(rawData[-1,-1]))


zero_impute <- function(rawData, nu_mat){
  nu_mat[is.na(nu_mat)] <- 0
  rawData[-1,-1] <- nu_mat
  return (rawData)
}

write.table(zero_impute(rawData, mx), "/Users/andyisman/Downloads/i54016212_zeroimpute_output", sep="\t", quote=F, row.names=F, col.names=F)

rowaverage_impute <- function(rawData, nu_mat) {
  impute <- function(x){
    x[is.na(x)] <- mean(x, na.rm=T)
    return (x)
  }
  rawData[-1,-1] <- t(apply(nu_mat, 1, impute))
  return (rawData)
}

write.table(rowaverage_impute(rawData, mx), "/Users/andyisman/Downloads/i54016212_rowaverageimpute_output", sep="\t", quote=F, row.names=F, col.names=F)
