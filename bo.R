# Function for BO algorithm
# Eric Yanchenko
# March 6, 2026


A <- B <- list()

for(i in 1:5){
  A[[i]] <- matrix(runif(5*2, 0, 10), nrow=5)
  B[[i]] <- matrix(runif(5*2, 0, 10), nrow=5)
}

# Computes the kernel matrix for EMD
# A is a list where each entry in the list is the location
# of the k sites
# Row of A is the sites included in the set
# Number of rows is number of sets
computeK <- function(A, weights=c("equal", "loc")){
  
  n = length(A) 
  
  # Compute weights
  if(weights=="equal"){
    wfunc <- function(aa){
      x = nrow(aa)
      return(rep(1/x, x))
    }
  }else if(weights=="loc"){
    wfunc <- function(aa){
      KA <- exp(-as.matrix(rdist(aa))) 
      KAinv <- solve(KA) 
      return(as.numeric(colSums(KAinv) / sum(KAinv) ))
    }
  }

  W <- lapply(A, wfunc)
  
  K <- matrix(0, n, n)
  diag(K) <- 1

  for(i in 2:n){
    for(j in 1:(i-1)){
      K[i,j] <- emdw(A[[i]], W[[i]], A[[j]], W[[j]])
      K[j,i] <- K[i,j]
    }
  }
  
  return(K)
  
}