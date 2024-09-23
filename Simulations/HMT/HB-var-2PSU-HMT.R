library(MASS)
library(MCMCpack)
library(fda)

# P-spline-based models 
HB_PSU2 <- function(Y, x, Pi, ID, p=2, L=7, mc=5000, bn=2000){
  ## preparation 
  x <- scale(x)   # scaled version of x will be used in MCMC
  N <- length(ID)
  H <- max(ID)
  ww <- (1/Pi)/mean(1/Pi)
  W <- diag(ww)
  X <- 1
  for(j in 1:p){
    X <- cbind(X, x^j)
  }
  qq <- dim(X)[2]
  Po <- function(x){ x*(x>0) }
  kappa <- seq(min(x), max(x), length=(L+2))[2:(L+1)]
  Z <- matrix(NA, N, L)
  for(l in 1:L){
    Z[,l] <- ( Po(x-kappa[l]) )^l
  }
  XX <- cbind(X, Z)
  
  ID_first <- rep(NA, H)  # index of the first unit in each stratum
  for(i in 1:H){
    ID_first[i] <- (1:N)[ID==i][1]
  }
  
  cor_mat <- function(rho){
    C <- diag(rep(1-rho, N))
    for(i in 1:H){
      C[ID==i, ID==i] <- C[ID==i, ID==i] + rho 
    }
    return(C)
  } 
  
  ## objects to store posterior samples
  Beta_pos <- matrix(NA, mc, qq+L)
  Gam_pos <- matrix(NA, mc, qq+L)
  rho_pos <- rep(NA, mc)
  
  ## initial values 
  Beta <- rep(0, qq+L)
  Gam <- rep(0, qq+L)
  Sig2 <- as.vector( exp(XX%*%Gam) ) 
  Tau_beta <- Tau_gam <- 1
  rho <- 0.5   # correlation parameter 
  
  ## MCMC 
  for(k in 1:mc){
    # Beta (mean)
    mat <- diag(c(rep(0.01, qq), rep(Tau_beta, L)))
    V_mat <- cor_mat(rho) 
    for(i in 1:H){
      V_mat[ID==i, ID==i] <- (Sig2[ID==i][1])*V_mat[ID==i, ID==i]
    }
    invS <- solve(V_mat + 10^(-2)*diag(N))
    A <- solve(mat + t(sqrt(ww)*XX)%*%invS%*%(sqrt(ww)*XX) + 10^(-5)*diag(p+L+1))
    A <- (A + t(A))/2    # to make A symmetric (although A is theoretically symmetric)
    B <- t(XX)%*%invS%*%(ww*Y)
    #Beta <- as.vector( A%*%B + t(chol(A))%*%rnorm(p+L+1) )
    Beta <- as.vector( mvrnorm(1, A%*%B, A) )
    Mu <- as.vector(XX%*%Beta)
    Beta_pos[k,] <- Beta
    
    # Tau_beta
    Tau_beta <- rgamma(1, 1+L/2, 1+sum(Beta[-(1:qq)]^2)/2)
    
    # Gam (variance) & rho 
    Gam_new <- Gam + 0.1*rnorm(qq+L)
    rho_new <- rho + 0.1*rnorm(1)
    if(rho_new<0){ rho_new <- 0 }
    if(rho_new>0.99){ rho_new <- 0.99 }
    Sig2_new <- as.vector( exp(XX%*%Gam_new) )
    V_mat_new <- cor_mat(rho_new) 
    for(i in 1:H){
      V_mat_new[ID==i, ID==i] <- (Sig2_new[ID==i][1])*V_mat_new[ID==i, ID==i]
    }
    invS_new <- solve(V_mat_new + 10^(-2)*diag(N))
    resid <- Y - Mu
    like <- like_new <- c()
    for(i in 1:H){
      like[i] <- (-0.5)*t(resid[ID==i])%*%invS[ID==i, ID==i]%*%(resid[ID==i]) + 0.5*log(det(invS[ID==i, ID==i]))
      like_new[i] <- (-0.5)*t(resid[ID==i])%*%invS_new[ID==i, ID==i]%*%(resid[ID==i]) + 0.5*log(det(invS_new[ID==i, ID==i]))
    }
    val <- sum(ww*like) + sum( dnorm(Gam[-(1:qq)], 0, 1/sqrt(Tau_gam), log=T) ) + sum( dnorm(Gam[1:qq], 0, 1, log=T) )
    val_new <- sum(ww*like_new) + sum( dnorm(Gam_new[-(1:qq)], 0, 1/sqrt(Tau_gam), log=T) ) + + sum( dnorm(Gam_new[1:qq], 0, 1, log=T) )
    pp <- exp(val_new - val)
    if(runif(1)<pp){ 
      Gam <- Gam_new
      Sig2 <- Sig2_new
      rho <- rho_new
    }
    Gam_pos[k,] <- Gam
    rho_pos[k] <- rho
    
    # Tau_gamma
    Tau_gam <- rgamma(1, 10+L/2, 1+sum(Gam[-(1:qq)]^2)/2)
  }
  
  # omit burn-in samples
  Beta_pos <- Beta_pos[-(1:bn),]
  Gam_pos <- Gam_pos[-(1:bn),]
  rho_pos <- rho_pos[-(1:bn)]
  Mu_pos <- t(XX%*%t(Beta_pos))[,ID_first]
  Sig2_pos <- t(exp(XX%*%t(Gam_pos)))[,ID_first]
  Sig2_pos[Sig2_pos>var(Y)] <- var(Y)   # adjustment 
  
  return(list(Beta=Beta_pos, Gam=Gam_pos, Mu=Mu_pos, S2=Sig2_pos, rho=rho_pos))
}
