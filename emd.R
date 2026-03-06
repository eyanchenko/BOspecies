# Earth mover's distance as kernel function
# Small EMD means the two sets are close, i.e., highly correlated
library(emdist)
library(rdist)


n = 5
# First column is weight
A <- matrix(runif(5*2, 0, 10), nrow=n)
B <- matrix(runif(5*2, 0, 10), nrow=n)

KA <- exp(-as.matrix(rdist(A))) 
KB <- exp(-as.matrix(rdist(B)))
KAinv <- solve(KA) 
KBinv <- solve(KB)

wA <- colSums(KAinv) / sum(KAinv) 
wB <- colSums(KBinv) / sum(KBinv) 
emdw(A, wA, B, wB)

# This weighting scheme ensures that if two points are close,
# their total weight is like if there was just one point
# Also compare performance with uniform weights

# Next steps
# Step 1: Convert code to R
# Step 2: Add EMD kernel function
  # Choose weights as in Garnett et al. (2010) paper
# Step 3: Compare EMD with Hamming distance