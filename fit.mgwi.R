# R code for fitting a MGWI trajectory

fit.mgwi <- function(x, st){
  
  pacman::p_load(VGAM, dplyr)
  
  # ZMG distribution
  
  dzmg <- function(x, prob, prob0){
    
    ifelse(x == 0, prob0 + (1 - prob0) * dgeom(0, prob = prob),
           (1 - prob0) * dgeom(x, prob = prob))
  }
  
  rzmg <- function(n, prob, prob0){
    
    B <- rbinom(n, size = 1, prob = 1 - prob0)
    
    Y <- rgeom(n, prob = prob)
    
    return(B * Y)
  }
  
  # Fitting MGWI model for stationary cases (without covariates)
  
  fit.st <- function(x){
    
    n <- length(x)
    
    #-------------------------------------
    # CLS estimation for initial guesses #
    #-------------------------------------
    
    # Q-function
    
    f_cls <- function(z){
      
      mu <- z[1]
      alpha <- z[2]
      
      # q-function
      
      sum((x[2:n] - alpha * (1 - (alpha / (1 + alpha))^(x[1:(n-1)]))
           - mu * (1 + mu) / (1 + mu + alpha))^2)
    }
    
    # Gradient vector
    g_cls <- function(z){
      
      mu <- z[1]
      alpha <- z[2]
      
      # grad functions
      
      grr1 <- function(par){
        
        x <- par[1:n]
        mu <- par[n+1]
        alpha <- par[n+2]
        
        -2 * (1 - alpha * (1 + alpha) / (1 + mu + alpha)^2) *
          sum(x[2:n] - alpha * (1 - (alpha / (1 + alpha))^x[1:(n-1)]) - mu * (1 + mu) / (1 + mu + alpha))
      }
      
      grr2 <- function(par){
        
        x <- par[1:n]
        mu <- par[n+1]
        
        alpha <- par[n+2]
        -2 * sum((x[2:n] - alpha * (1 - (alpha / (1 + alpha))^x[1:(n-1)]) -
                    mu * (1 + mu)/(1 + mu + alpha)) * (1 - (alpha / (1 + alpha))^x[1:(n-1)] * (1 + x[1:(n-1)] / (1 + alpha))
                                                       - mu * (1 + mu) / (1 + mu + alpha)^2))
      }
      
      c(grr1(c(x, mu, alpha)), grr2(c(x, mu, alpha)))
    }
    
    fit_cls <- optim(par = c(1,1), fn = f_cls, gr = g_cls, method = "BFGS", hessian = FALSE)
    
    est_cls <- fit_cls$par
    
    #============================================
    # Conditional maximum likelihood estimation #
    #============================================
    
    loglik <- function(z){
      
      mu <- z[1]
      alpha <- z[2]
      
      m <- cbind(x[1:(n-1)], x[2:n])
      
      # Transition probabilities
      
      tp <- function(z){
        
        x <- z[1]
        y <- z[2]
        
        if(x == 0){
          res <- dzmg(y, 1/(1 + mu), alpha/(1 + mu + alpha))
          return(res)
        }
        else{
          
          # P(alpha min_operator x = z)
          
          P <- function(x, w = 0:y){
            
            case_when(
              x < w ~ 0,
              x == w ~ (alpha / (1 + alpha))^w,
              x > w ~ alpha^w / (1 + alpha)^(w + 1)
            )
          }
          
          res <- P(x, 0:y) * dzmg(y - 0:y, 1 / (1 + mu), alpha/(1 + mu + alpha))
          
          return(sum(res))
          
        }
      }
      
      # likelihood function
      return(sum(log(apply(m, 1, tp))))
    }
    
    par <- c(1,1)
    
    fit_cml <- optim(par = par, fn = loglik, gr = NULL, method = "Nelder-Mead",
                     control = list(fnscale = -1), hessian = TRUE)
    
    est_cml <- fit_cml$par
    
    return(list("Coef_cml" = est_cml, "Conv_cml" = fit_cml$convergence,
                "Coef_cls" = est_cls, "Conv_cls" = fit_cls$convergence))
  }
  
  # Fitting MGWI model for non stationary cases (with covariates)
  
  fit.nst <- function(x){
    
    n <- length(x)
    
    # CLS estimation
    
    # Dimensions of p and q
    p <- 2
    q <- 2
    
    # Sn function 
    f_cls <- function(z){
      
      beta <- z[1:p]
      gamma <- z[(p+1):(p+q)]
      
      t <- seq(1:n)
      
      w <- rbind(array(1, c(1,n)), t/n) # covariates for mu
      eta <- t(w) %*% beta # linear predictor for mu
      mu <- exp(eta) # structure for mu
      
      v <- rbind(array(1, c(1,n)), t/n) # covariates for alpha
      nu <- t(v) %*% gamma # linear predictor for alpha
      alpha <- exp(nu) # structure for alpha
      
      # s-function
      Sn <- sum((x[2:n] - alpha[2:n] * (1 - (alpha[2:n]/(1 + alpha[2:n]))^(x[1:(n-1)])) 
                 - mu[2:n] * (1 + mu[2:n])/(1 + mu[2:n] + alpha[2:n]))^2)
      
      return(Sn)
    }
    
    # Gradient vector
    g_cls <- function(z){
      
      beta <- z[1:p]
      gamma <- z[(p+1):(p+q)]
      
      t <- seq(1:n)
      
      w <- rbind(array(1, c(1,n)), t/n) # covariate values for mu
      eta <- t(w) %*% beta # linear predictor for mu
      mu <- exp(eta) # structure for mu
      
      v <- rbind(array(1, c(1,n)), t/n) # covariate values for alpha
      nu <- t(v) %*% gamma # linear predictor for alpha
      alpha <- exp(nu) # structure for alpha
      
      # gradient functions
      
      grad1 <- function(par){
        
        x <- par[1:n]
        beta <- par[(n+1):(n+p)]
        gamma <- par[(n+p+1):(n+p+q)]
        
        -2 * sum((mu[2:n] * (1 + alpha[2:n] + mu[2:n] * (2 + 2 * alpha[2:n] + mu[2:n])) / (1 + mu[2:n] + alpha[2:n])^2) * 
                   (x[2:n] - alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - mu[2:n] * (1 + mu[2:n])/(1 + mu[2:n] + alpha[2:n])))
      }
      
      grad2 <- function(par){
        
        x <- par[1:n]
        beta <- par[(n+1):(n+p)]
        gamma <- par[(n+p+1):(n+p+q)]
        
        -2 * sum(((t[2:n] / n) * mu[2:n] * (1 + alpha[2:n] + mu[2:n] * (2 + 2 * alpha[2:n] + mu[2:n])) / (1 + mu[2:n] + alpha[2:n])^2) * 
                   (x[2:n] - alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - mu[2:n] * (1 + mu[2:n])/(1 + mu[2:n] + alpha[2:n])))
      }
      
      grad3 <- function(par){
        
        x <- par[1:n]
        beta <- par[(n+1):(n+p)]
        gamma <- par[(n+p+1):(n+p+q)]
        
        -2 * sum((alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - alpha[2:n] * mu[2:n] * (1 + mu[2:n]) / (1 + mu[2:n] + alpha[2:n])^2 -
                    x[1:(n-1)] * (alpha[2:n] / (1 + alpha[2:n]))^(x[1:(n-1)] + 1)) * 
                   (x[2:n] - alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - mu[2:n] * (1 + mu[2:n])/(1 + mu[2:n] + alpha[2:n])))
      }
      
      grad4 <- function(par){
        
        x <- par[1:n]
        beta <- par[(n+1):(n+p)]
        gamma <- par[(n+p+1):(n+p+q)]
        
        -2 * sum(((t[2:n] / n) * alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - (t[2:n] / n) * alpha[2:n] * mu[2:n] * 
                    (1 + mu[2:n]) / (1 + mu[2:n] + alpha[2:n])^2 -
                    (t[2:n] / n) * x[1:(n-1)] * (alpha[2:n] / (1 + alpha[2:n]))^(x[1:(n-1)] + 1)) * 
                   (x[2:n] - alpha[2:n] * (1 - (alpha[2:n] / (1 + alpha[2:n]))^x[1:(n-1)]) - mu[2:n] * (1 + mu[2:n])/(1 + mu[2:n] + alpha[2:n])))
      }
      
      c(grad1(c(x, beta, alpha)), grad2(c(x, beta, alpha)), grad3(c(x, beta, alpha)), 
        grad4(c(x, beta, alpha)))
    }
    
    par <- rep(1, 4)
    
    fit_cls <- optim(par = par, fn = f_cls, gr = g_cls, method = "Nelder-Mead", hessian = FALSE)
    
    est_cls <- fit_cls$par
    
    # ML estimation via transition probabilities
    loglik <- function(z){
      
      beta <- z[1:p]
      gamma <- z[(p+1):(p+q)]
      
      t <- seq(1:n)
      
      w <- rbind(array(1, c(1,n)), t/n) # covariates for mu
      eta <- t(w) %*% beta # linear predictor for mu
      mu <- exp(eta) # structure for mu
      
      v <- rbind(array(1, c(1,n)), t/n) # covariates for alpha
      nu <- t(v) %*% gamma # linear predictor for alpha
      alpha <- exp(nu) # structure for alpha
      
      # Transition probabilities
      tp <- function(z){
        
        x <- z[1] # x_(t-1)
        y <- z[2] # x_t
        w <- z[3] # mu
        u <- z[4] # alpha
        
        if(x == 0){
          
          res <- dzigeom(y, 1/(1 + w), u/(1 + w + u))
          return(res)
        }
        
        else{
          
          P <- function(x, w = 0:y, u){  # P(alpha min_operator x = z)
            
            case_when(
              x < w ~ 0,
              x == w ~ (u / (1 + u))^w,
              x > w ~ dgeom(w, 1 / (1 + u))
            )
          }
          
          res <- P(x, 0:y, u) * dzigeom(y - (0:y), 1/(1 + w), u/(1 + w + u))
          
          return(sum(res))
        }
      }
      
      m <- cbind(x[1:(n-1)], x[2:n], mu[2:n], alpha[2:n])
      
      # likelihood function
      return(sum(log(apply(m, 1, tp))))
    }
    
    par <- rep(1, p+q) # initial values
    
    fit_cml <- optim(par = par, fn = loglik, gr = NULL, method = "BFGS",
                     control = list(fnscale = -1),
                     hessian = TRUE)
    
    est_cml <- fit_cml$par
    
    # Output
    return(list("Coef_cml" = est_cml, "Conv_cml" = fit_cml$convergence,
                "Coef_cls" = est_cls, "Conv_cls" = fit_cls$convergence))
  }
  
  # Model selection
  
  if(st == TRUE) return(fit.st(x))
  
  if(st == FALSE) return(fit.nst(x))
  
}