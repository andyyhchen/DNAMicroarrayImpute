iknn = function(xmiss, Niter, K=NULL, e=1e-3) {
  
  # KNN impute Algorithm
  # 
  # Some definitions:
  # Target gene: gene with missing entries.
  # Candidate genes: genes that are candidates to be the k-nearest neighbors
  # of the target gene.
  # 
  # Routine Details:
  # In a 1st round, the missing entries are estimated by row (gene) average.
  # Next, Niter rounds of KNN imputations are performed.
  # In each cycle, the set of candidate genes is composed by the last imputed
  # matrix. So the candidate set as the maximum possible dimension. 
  # 
  # If K is not given, before the KNN imputation, the optimum number of neighbors 
  # is determined automatically using estimateK
  # 
  # Additionally, all the missing entries of a given target gene are imputed
  # at once. Therefore, for each target gene, its missing experiments are
  # eliminated in the target and in the set of candidate genes. So, the
  # target gene and the set of candidate genes share the same dimension (i.e.
  # they have the same n of experiments).
  # Because of this the simple Euclidean distance can be used.
  # 
  # After selecting the k-nearest neighbors, the missing value is imputed
  # using a weighted average between these genes. The weights for each
  # K-nearest neighbor are given by:
  # wi = (1/di) / (sum of all di, with i=1,...,K)   
  # so each missing entry is determined by:
  #
  # ximputed = sum of wi * xi, with i=1,...,K  
  #  
  # INPUTS:
  # xmiss - matrix with missing values (NA values): p genes x n arrays
  # !! GENES IN ROWS  !!!
  # Niter - number of iterations to perform
  # K - (optional input) number of neighbors
  # e - convergence tolerance (default: 1e-3)
  #
  # OUTPUTS:
  # xcomplete -  complete matrix with missing entries estimated by IKNN (p genes x n arrays)
  # Kopt - number of neighbors (defined by the user or estimated automatically)
  # 
  # LPB Jan 2005
  
  similarityCal<-function(vec, mat, method="EuDist"){
    methods<-c("EuDist","cor","Angle")
    switch(match.arg(method,methods),
           EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
           cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
           Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
    )
  }
  # 1) Impute missing values by row average:
  Eold <- rowMeanSubstitution(xmiss) 
  
  if (is.null(K)) Kopt <- estimateK(xmiss, Eold) else Kopt <- K
  
  # 1st Cycle: Impute missing values by KNN imputation:
  
  # determine the columns (genes) with missing values (miss_gene):
  R <- !is.na(xmiss) # missing indicator matrix (1 if entry is observed, 0 if it is missing)
  
  miss_gene <- which(rowSums(!R)!=0)
  miss_idx <- which(!R)
  
  cat("\n\n")
  cat("----------------------------------------------\n")
  cat("KNNimputation/Round 1 \n")
  cat("----------------------------------------------\n")
  cat("\n\n")
  
  
  xcomplete <- Eold 
  
  # arrays with missing entries for each target gene i
  miss_exp <- lapply(miss_gene, function(i) which(is.na(xmiss[i,])))
  
  E <- xmiss
  
  for (j in 1:length(miss_gene)){
    
    gene <- xmiss[miss_gene[j],] # target gene
    exp <- miss_exp[[j]] # row indexes (experiments) with missing entries within the target gene  
    gene <- gene[-exp] # remove the missing entries
    
    # Remove the experiments corresponding to the missing entries of the target gene and delete the target gene from the candidate set
    cand_x <- xcomplete[-miss_gene[j], -exp]
    # ... compute the Euclidean distance between each candidate gene and
    # the target gene j:
    #d <- sqrt(rowSums((cand_x - matrix(gene, ncol=length(gene),nrow=nrow(cand_x), byrow=TRUE))^2))
    
    d <- similarityCal(gene, cand_x, "cor")
    # determine the candidate genes corresponding to the Kopt "smallest"
    # distances:
    ds <- d[order(d, decreasing=T)[1:Kopt]]
    idx <- order(d, decreasing=T)[1:Kopt]
    idx_gene <- c(1:nrow(xmiss))[-miss_gene[j]]
    idx_gene <- idx_gene[idx] # k similar genes
    y <- xcomplete[idx_gene,, drop=FALSE] 
    
    const <- 1/sum(ds) # normalising weight constant 
    w <- const * ds #weights for each gene
    
    E[miss_gene[j], exp] <- w%*%y[, exp, drop=FALSE] #weighted average of the k-nearest neighbors
  }
  error <- sum((E[miss_idx]-Eold[miss_idx])^2)
  cat(sprintf("Sum of squares between the last 2 imputation cycles = %g \n",error))
  
  
  # 2nd and next cycles: Repeat the imputations using the last imputed matrix
  for (h in 2:Niter){
    
    cat("\n\n")
    cat("----------------------------------------------\n")
    cat(sprintf("KNNimputation/Round %s\n", h))
    cat("----------------------------------------------\n")
    cat("\n\n")
    
    Eold <- E
    error_old <- error
    xcomplete <- Eold
    
    for (j in 1:length(miss_gene)){
      gene <- xmiss[miss_gene[j],] # target gene
      exp <- miss_exp[[j]] # row indexes (experiments) with missing entries within the target gene  
      gene = gene[-exp] # remove the missing entries
      
      cand_x <- xcomplete[-miss_gene[j], -exp]
      
      #d <- sqrt(rowSums((cand_x - matrix(gene, ncol=length(gene),nrow=nrow(cand_x), byrow=TRUE))^2))
      d <- similarityCal(gene, cand_x, "cor")
      # determine the candidate genes corresponding to the Kopt "smallest"
      # distances:
      ds <- d[order(d, decreasing=T)[1:Kopt]]
      idx <- order(d, decreasing=T)[1:Kopt]
      idx_gene <- c(1:nrow(xmiss))[-miss_gene[j]]
      idx_gene <- idx_gene[idx] # k similar genes
      y <- xcomplete[idx_gene,, drop=FALSE] 
      
      const <- 1/sum(ds) # normalising weight constant 
      w <- const*ds #weights for each gene
      
      E[miss_gene[j], exp] <- w%*%y[, exp, drop=FALSE] #weighted average of the k-nearest neighbors
    }
    
    error <- sum((E[miss_idx]-Eold[miss_idx])^2) 
    cat(sprintf("Sum of squares between the last 2 imputation cycles = %g \n",error))
    cat(sprintf("Ratio of the sum of squares between the last 2 imputation cycles = %g \n",error_old/error))
    
    if (( error_old/error < 4 ) | (error < e ) ) break 
  }
  
  
  cat(sprintf("A total of %g rounds were performed \n",h))
  
  return(out = list(xcomplete=E, kest = Kopt))
}



## =======================================
rowMeanSubstitution = function(xmiss){
  
  # Impute missing entries by row mean
  # genes in rows
  # 
  # E = Rmeanimpute(data);
  # 
  # INPUTS:
  # xmiss - matrix with missing entries (genes x arrays)
  # 
  # OUTPUTS:
  # E - complete data matrix with missing entries imputed by row (gene) mean
  # 
  # LPB Fev 2005
  
  rowM = rowMeans(xmiss, na.rm=TRUE)
  rowM = matrix(rowM, ncol=ncol(xmiss), nrow=nrow(xmiss), byrow=FALSE)
  
  E = xmiss
  E[is.na(xmiss)] = rowM[is.na(xmiss)] 
  
  return(E)
}
## =======================================
nrmse <- function(missing, guess, answer) { 
  x = is.na(missing)
  sqrt(mean((guess[x]-answer[x])^2)/var(answer[x])) 
}
answer = as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/ans.txt", sep="\t"))
missing <- as.matrix(read.table("/Users/andyisman/Documents/BioInfo/lymphoma/m_15_1.txt", sep=""), sep="\t")

cat(nrmse(missing, iknn(missing, 2, 15)$x, answer),"\n")