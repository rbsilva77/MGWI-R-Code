#-------------------------------------------------------
# R code for generating a non-stationary MGWI trajectory
#-------------------------------------------------------

# Generating count responses based on non-stationary MGWI processes

# Instructions

# Set values for beta0, beta1, beta2, gamma0 and gamma1

# Set the sample size n

  r.mgwi_nst <- function(size, par){
    
    n <- size; b0 <- par[1]; b1 <- par[2]; b2 <- par[3]; g0 <- par[4]; g1<- par[5] 
    
    #----------------------------------------------------------------
    # Zero-modified geometric (ZMG) distribution
    
    dzmg <- function(x, prob, prob0){
      
      ifelse(x == 0, prob0 + (1 - prob0) * dgeom(0, prob = prob),
             (1 - prob0) * dgeom(x, prob = prob))
    }
    
    rzmg <- function(n, prob, prob0){
      
      B <- rbinom(n, size = 1, prob = 1 - prob0)
      
      Y <- rgeom(n, prob = prob)
      
      return(B * Y)
    }
    #----------------------------------------------------------------
    
    t <- seq(1:n)
    
    beta <- c(b0, b1, b2)
    w <- rbind(array(1, c(1, n)), t / n, cos(2 * pi * t / 12)) # covariates for mu
    eta <- t(w) %*% beta # linear predictor for mu
    mu <- exp(eta) # theoretical structure for mu
    
    gamma <- c(g0, g1)
    v <- rbind(array(1, c(1, n)), t / n) # covariates for alpha
    nu <- t(v) %*% gamma # linear predictor for alpha
    alpha <- exp(nu) # theoretical structure for alpha
    
    # Generation of the MGWI process
    # Obs. # Xalpha ~ Geom(1 / (1 + alpha)); epsilon ~ ZMG(1 / ( 1 + mu), alpha / (1 + mu + alpha))
    
    x <- array(0,c(n,1))
    
    x[1] <- rgeom(1, 1 / (1 + mu[1]))
    for (t in 2:n){
      Xalpha <- rgeom(1, 1 / (1 + alpha[t]))
      x[t] <- min(x[t-1], Xalpha) + rzmg(1, prob = 1 / (1 + mu[t]), prob0 = alpha[t]/(1 + mu[t-1] + alpha[t]))
    }
    
    return(x)
  } 
  
  # Example
  
  n <- 100 # sample size
  
  b0 <- 2; b1 <- 1; b2 <- .7; g0 <- 2; g1 <- 1 # parameter values
  
  x <- r.mgwi_nst(size = n, par = c(b0, b1, b2, g0, g1))
  
  x