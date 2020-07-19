


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