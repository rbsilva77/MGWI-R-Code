#---------------------------------------------------
# R code for generating a stationary MGWI trajectory
#---------------------------------------------------

# Generating count responses based on a stationary MGWI process

# Instructions

# Set values for mu and alpha

# Set the sample size n

r.mgwi_st <- function(size, par){
  
  n <- size
  mu <- par[1]
  alpha <- par[2]
  
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
  
  x <- NULL
  x[1] <- rgeom(1,1 / (1 + mu))
  for (i in 2:n){
    Xalpha <- rgeom(1, 1 / (1 + alpha)) #  X ~ Geo(alpha)
    x[i] <- min(x[i-1], Xalpha) + rzmg(1,1 / (1 + mu), alpha / (1 + mu + alpha))
  }
  return(x)
}

# Example

n <- 100 # sample size

mu <- 2; alpha <- 1 # parameter values

x <- r.mgwi_st(size = n, par = c(mu, alpha))

x