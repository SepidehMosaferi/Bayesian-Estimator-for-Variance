library(MASS)
library(MCMCpack)
library(fda)

# P-spline-based models 
# survey weight is incorporated as normalized weight for likelihood
HB_var <- function(Y, x, Pi, p=2, L=7, mc=5000, bn=2000, verbose=F){
  ## preparation 
  x <- scale(x) # scaled version of x will be used in MCMC
  ww <- (1/Pi)/mean(1/Pi)
  H <- length(Y)
  X <- 1
  for(j in 1:p){
    X <- cbind(X, x^j)
  }
  qq <- dim(X)[2]
  Po <- function(x){ x*(x>0) }
  kappa <- seq(min(x), max(x), length=(L+2))[2:(L+1)]
  Z <- matrix(NA, H, L)
  for(l in 1:L){
    Z[,l] <- ( Po(x-kappa[l]) )^l
  }
  XX <- cbind(X, Z)
  
  ## objects to store posterior samples
  Beta_pos <- matrix(NA, mc, qq+L)
  Gam_pos <- matrix(NA, mc, qq+L)
  
  ## initial values 
  Beta <- rep(0, qq+L)
  Gam <- rep(0, qq+L)
  Sig2 <- as.vector( exp(XX%*%Gam) ) 
  Tau_beta <- Tau_gam <- 1
  
  ## MCMC 
  for(k in 1:mc){
    # Beta (mean)
    mat <- diag(c(rep(0.01, qq), rep(Tau_beta, L)))
    A <- solve(mat + t(XX)%*%(ww*XX/Sig2) + 10^(-5)*diag(p+L+1))
    A <- (A + t(A))/2    # to make A symmetric (although A is theoretically symmetric)
    B <- t(XX)%*%(ww*Y/Sig2)
    Beta <- as.vector( mvrnorm(1, A%*%B, A) )
    Mu <- as.vector(XX%*%Beta)
    Beta_pos[k,] <- Beta
    
    # Tau_beta
    Tau_beta <- rgamma(1, 1+L/2, 1+sum(Beta[-(1:qq)]^2)/2)
    
    # Gam (variance)
    Gam_new <- Gam + 0.1*rnorm(qq+L)
    Sig2_new <- as.vector( exp(XX%*%Gam_new) )
    like <- sum( ww*dnorm(Y, Mu, sqrt(Sig2), log=T)) + sum( dnorm(Gam[-(1:qq)], 0, 1/sqrt(Tau_gam), log=T) ) + sum( dnorm(Gam[1:qq], 0, 1, log=T) )
    like_new <- sum( ww*dnorm(Y, Mu, sqrt(Sig2_new), log=T)) + sum( dnorm(Gam_new[-(1:qq)], 0, 1/sqrt(Tau_gam), log=T) ) + + sum( dnorm(Gam_new[1:qq], 0, 1, log=T) )
    pp <- exp(like_new - like)
    if(runif(1)<pp){ 
      Gam <- Gam_new
      Sig2 <- Sig2_new
    }
    Gam_pos[k,] <- Gam
    
    # Tau_gamma
    Tau_gam <- rgamma(1, 10+L/2, 1+sum(Gam[-(1:qq)]^2)/2)
    
    # print 
    if(verbose & round(10*k/mc)==(10*k/mc)){
      print( paste0("MCMC ", 100*k/mc, "% completed.") )
      print( Gam )
    }
  }
  
  # omit burn-in samples
  Beta_pos <- Beta_pos[-(1:bn),]
  Mu_pos <- t( XX%*%t(Beta_pos) )
  Gam_pos <- Gam_pos[-(1:bn),]
  Sig2_pos <- t( exp(XX%*%t(Gam_pos)) )
  Sig2_pos[Sig2_pos>var(Y)] <- var(Y)   # adjustment 
  
  return(list(Beta=Beta_pos, Gam=Gam_pos, Mu=Mu_pos, S2=Sig2_pos))
}
