# Function for BO algorithm
# Eric Yanchenko
# March 6, 2026
library(extraDistr)

# Computes the kernel matrix
# Uses exponential kernel where distance is EMD
# X is a list where each entry in the list is the location
# of the k sites (allows for other covariates with the locations here)
# Number of entries in X is number of sets
computeK <- function(X, weights=c("equal", "loc"), lambda=1, rho=1){
  
  n = length(X) 
  
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

  W <- lapply(X, wfunc)
  
  K <- matrix(0, n, n)
  diag(K) <- 1

  for(i in 2:n){
    for(j in 1:(i-1)){
      d <- emdw(X[[i]], W[[i]], X[[j]], W[[j]])
      K[i,j] <- lambda*exp(-d/rho)
      K[j,i] <- K[i,j]
    }
  }
  
  # Return covariance matrix and weights
  return(list(K=K, W=W))
  
}

# MCMC sampler for beta0 and sigma2
mcmc <- function(Y, K, mc_burn = 1000, mc_iters = 1000){
  N         <- length(Y)
  Vinv  <- solve(K + diag(N))
  Vsum <- sum(Vinv)
  Vy    <- sum(Vinv%*%Y)
  
  # Initialize MCMC
  beta0 = 0
  sigma2 = rinvgamma(1,1,1)
  
  post <- matrix(0, nrow=mc_iters, ncol=2)
  
  for(i in 1:(mc_burn+mc_iters)){
    # Sample from beta0
    vvv   <- 1/(Vsum/sigma2 + 1/100)
    mmm   <- Vy/sigma2

    beta0 <- rnorm(1, mmm*vvv, sqrt(vvv))
    
    # Sample from sigma2
    sigma2 <- rinvgamma(1, N/2-1, 0.5*t(Y-beta0)%*%Vinv%*%(Y-beta0))
    
    if(i > mc_burn){
      post[i-mc_burn,] <- c(beta0, sigma2)
    }
    
  }
  
  return(apply(post, 2, median))
}

# Calculate covariance of new point with existing points
new_cov <- function(xstar, X, K, weights=c("equal", "loc"), W, lambda, rho){
  N <- ncol(K)
  Vinv <- solve(K + diag(N))
  
  # Calculate covariance of new point with existing points
  # First, get weights
  if(weights=="equal"){
    Wstar <- rep(1/nrow(xstar), nrow(xstar))
  }else if(weights=="loc"){
    Kstar <- exp(-as.matrix(rdist(xstar))) 
    Kstarinv <- solve(Kstar) 
    Wstar <- as.numeric(colSums(Kstarinv) / sum(Kstarinv) )
  }
  # Then calculate covariance
  kstar <- unlist(lapply(1:N, function(i){
    lambda*exp(-emdw(xstar, Wstar, X[[i]], W[[i]])/rho)
  }))
  return(
    list(kstar=kstar, Wstar=Wstar)
  )
}

# Convert mean and cov for new point
mean_cov <- function(xstar, X, Y, beta0, sigma2, K, weights=c("equal", "loc"), W, lambda=1, rho=1){
  N <- ncol(K)
  Vinv <- solve(K + diag(N))
  
  #  Calculate covariance
  kstar <- new_cov(xstar, X, K, weights, W, lambda, rho)$kstar
  
  mu = as.numeric(beta0 + kstar%*%Vinv%*%(Y-beta0))
  sd = as.numeric(sqrt(sigma2*(2-kstar%*%Vinv%*%kstar)))
  
  return(list(mu=mu, sd=sd))
}

# Calculate Expected improvement
AEI <- function(xstar, X, Y, ystar, beta0, sigma2, K, weights=c("equal", "loc"), W, lambda=1, rho=1){
  
  out <- mean_cov(xstar, X, Y, beta0, sigma2, K, weights, W, lambda, rho)
  mu = out$mu
  sd = out$sd
  
  Delta = mu - ystar
  
  EI = max(Delta, 0) + sd * dnorm(Delta/sd) - abs(Delta)*pnorm(Delta/sd)

  return(EI*(1-sqrt(sigma2)/sqrt(sd^2+sigma2)))
  
}

# Greedy function to maximize AEI
greedy <- function(locs, X, Y, V, k, ystar, beta0, sigma2, K, weights, W, lambda, rho){
  
  
  # True if candidate set, x, has already been evaluated. False otherwise
  # Could probably figure out a way to speed this up
  # Like have a separate data structure which just holds site numbers of sets
  check_fun <- function(x){
    x <- sort(x)
    
    apply_fun <- function(Vj){
      return(sum(x==Vj)==k)
    }
    
    sum(unlist(apply(V, 1, apply_fun))) > 0 
    
  }
  
  # Initialize set of sites
  L <- nrow(locs)
  xstar_id <- sample(1:L, k)
  xstar <- locs[xstar_id, ]
  
  # Find current EI
  aei_cur <- AEI(xstar, X, Y, ystar, beta0, sigma2, K, weights, W, lambda, rho)
  ind = TRUE
  
  while(ind){
    ind = FALSE
    
    # Swap ith site with all other sites and keep largest
    for(i in 1:k){
      # Randomize order
      idx <- sample((1:L)[-xstar_id] )
      
      for(j in idx){
        xtry <- xstar
        xtry_id <- xstar_id
        
        # Swap current ith site with all other sites
        xtry[i,] <- locs[j, ]
        xtry_id[i] <- j
        
        
        # Check if we have already sampled this set
        ifelse(check_fun(xtry_id), aei_try <- -Inf, aei_try <- AEI(xtry, X, Y, ystar, beta0, sigma2, K, weights, W, lambda, rho))
        
        if(aei_try > aei_cur){
          aei_cur = aei_try
          xstar = xtry
          xstar_id[i] <- j 
          ind = TRUE
        }
      }
    }
  }
  
  return(list(xstar=xstar, id=sort(xstar_id)))
}

# Main algorithm
# f is objective function. Make sure inputs for it work lol
# This seems to work!
bayesopt <- function(f, locs, k, weights=c("equal", "loc"), N0=5, B=20, lambda=1, rho=1, mc_burn=1000, mc_iters=5000){
  
  L = nrow(locs) # total number of locations
  p = ncol(locs) # number of covariates at each site
  
  # Initial random sample
  X <- list() # keeps track of covariates of sampled places
  V <- matrix(0, nrow=N0+B, ncol = k)# keeps track of IDs of sampled places
  for(i in 1:N0){
    ids <- sort(sample(1:L, k))
    X[[i]] <- locs[ids, ]
    V[i, ] <- ids
  }

  # Evaluate objective function at each set
  Y <- unlist(lapply(X, f))
  meany = mean(Y)
  Y <- Y - meany
  
  # Construct covariance matrix
  out <- computeK(X, weights, lambda, rho)
  K = out$K
  W = out$W
  
  # Repeat while budget remaining
  for(b in 1:B){
    out <- mcmc(Y, K, mc_burn, mc_iters)
    b0 = out[1]
    sigma2 = out[2]
    
    # Find current optimum as largest conditional mean (removes uncertainty of evaluations)
    ystar <- max(unlist(lapply(X, function(xxx){mean_cov(xxx, X, Y, b0, sigma2, K, weights, W, lambda, rho)$mu})))
    
    # Find set of locations which maximize AEI
    out <- greedy(locs, X, Y, V, k, ystar, b0, sigma2, K, weights, W, lambda, rho)
    xstar = out$xstar
    xstar_id = out$id
    
    # Append xstar and evaluate objective function
    Y <- c(Y, f(xstar)-meany)
    X[[N0+b]] <- xstar
    V[N0+b, ] <- xstar_id
    
    # Update K and W
    out <- new_cov(xstar, X, K, weights, W, lambda, rho)
    kstar <- out$kstar
    Wstar <- out$Wstar
    K <- rbind(K, kstar)
    K <- cbind(K, c(kstar,1))
    rownames(K) <- NULL
    
    W[[N0+b]] <- Wstar
    print(b)
  }
  
  # Select largest value
  yvec  <- unlist(lapply(X, function(xxx){mean_cov(xxx, X, Y, b0, sigma2, K, weights, W, lambda, rho)$mu}))
  idmx  <- which.max(yvec)
  xstar <- X[[idmx]]
  
  # Return trajectory
  ytraj <- numeric(B+1)
  ytraj[1] <- max(yvec[1:N0])
  for(b in 1:B){ytraj[b+1] <- max(yvec[1:(N0+b)])}
  
  return(list(xstar=xstar, y=ytraj))
}




# All possible locations to choose from
L = 500
locs <- matrix(runif(2*L, -1, 1), nrow=L, ncol=2)

# Crazy objective function
f <- function(x){
  -sin(12*x[1,1])*x[1,1]*x[2,1] + 0.5*x[1,1]^3*x[2,2] + sin(6*x[1,2]*x[2,2]) - x[2,1]*x[1,2]^4 +
    cos(x[3,1]*x[3,2])*x[1,2] + rnorm(1, 0, 0.01)
}



f(locs[sample(1:L, 3), ])

k = 3 # number of sites considered in each set


out <- bayesopt(f, locs, 3, "loc", N0=5, B=25)
f(out$xstar)
plot(out$y, type="l")


x <- numeric(10); for(i in 1:10){x[i] <- f(bayesopt(f, locs, 3, "loc", N0=5, B=10)$xstar)}; mean(x) # Optimal input
z <- numeric(100); for(i in 1:100){z[i] <- f(locs[sample(1:L, k),])}; mean(z) # Average over random inputs

# So we are doing significantly better than random guessing. Good sanity check.



