# ==============================================================================
# R Code for the Simulation Studies (HMT case SYSTEMATIC with 2 PSUs)
# Date: Feb 2024 
# ==============================================================================

# Required Libraries
require(PracTools); require(sampling); require(LaplacesDemon)
require(pscl); require(MCMCpack); require(extraDistr); require(jipApprox)

setwd("/Users/sepidehmosaferi/Documents/GitHub/Bayesian-Estimator-for-Variance/Simulations/HMT")
rm(list = ls(all = TRUE))
source("HB-var-2PSU-HMT.R")

sigma2 <- 0.0625 #from HMT (1983) paper
set.seed(2012)
DATA <- as.data.frame(HMT(N=2000,H=20))
DATA$x <- (DATA$x-min(DATA$x))/(max(DATA$x)-min(DATA$x))
DATA <- DATA[-1,]

H <- 20
Ni <- as.vector(table(DATA$strat))
N <- sum(Ni)

# first-order inclusion probabilities (for 2 units sample selection):
pik <- sapply(1:H,function(i){inclusionprobabilities(DATA$x[DATA$strat==i],2)}) 
DATA$pik <- unlist(pik)

# second-order (joint) inclusion probabilities:
pik2 <- sapply(1:H,function(i){jip_approx(DATA$pik[DATA$strat==i],method='Hajek')})

# Variance for the purpose of bias and MSE:
S2_h <- sapply(1:H,function(i){HTvar(DATA$y[DATA$strat==i],pik2[[i]],sample=FALSE,method="HT")})
true_Var <- sum(S2_h)/N^2 

# Sample Sizes (selecting PSUs)
ni_PSU <- rep(2,H)  # vector of stratum sample sizes 

# Target Stratum Mean 
Meanstratum <- tapply(DATA$y,INDEX =DATA$strat,FUN=mean)

# Stratified Mean PSU design
overall_mean <- sum((Ni/N)*as.vector(Meanstratum))

## Intro. Repetition Loop
one.rep <- function(){
  cat("start rep",date(),"\n")
  
  # Sample selection process (selecting PSUs)
  SMP.IDs_PSU <- strata(data=DATA,stratanames="strat",size=ni_PSU,method="systematic",pik=DATA$pik)
  SAMPLE_PSU <- getdata(DATA,SMP.IDs_PSU)  # Output of selected sample
  
  ### Collapsed Variance estimator
  ng <- 2*ni_PSU[1]
  CollIndex <- rep(1:(H/2),each=ng)
  SAMPLE_PSU$CollIndex <- CollIndex 
  
  tHT_PSU <- rep(0,H)
  for (i in 1:H)
  {
    tHT_PSU[i] <- sum(SAMPLE_PSU$y[SAMPLE_PSU$Stratum==i]/
                        SAMPLE_PSU$Prob[SAMPLE_PSU$Stratum==i])
  }
  VarColl_each_PSU <- rep(0,H/2)
  for(i in 1:(H/2)){
    VarColl_each_PSU[i] <- 0.5*sum(tHT_PSU[2*i-1]^2,tHT_PSU[2*i]^2)
  }
  VarColl_PSU <- sum(VarColl_each_PSU)/(N^2)
  
  ### H-T mean 
  mean_PSU <- sum(tHT_PSU)/N 
  
  ### Kernel-Based Variance estimator
  h <- 0.06 # h should be in the range of [1/H,2/H].
  
  Kernel_PSU <- function(x,xi,h){u<-(x-xi)/h; k<-0.75*(1-u^2); k*(k>0)}
  
  x_ker <- sapply(1:H,function(i){sum(SAMPLE_PSU$x[SAMPLE_PSU$Stratum==i])}) 
  
  d_ji_PSU <- matrix(0,nrow=H,ncol=H) #column is stratum
  for(i in 1:H){
    d_ji_PSU[i,] <- Kernel_PSU(x_ker,x_ker[i],h)
    d_ji_PSU[i,] <- d_ji_PSU[i,]/sum(d_ji_PSU[i,]) 
  }
  
  Cd_PSU <-mean(rep(1,H)-2*diag(d_ji_PSU)+apply(d_ji_PSU^2,MAR=1,FUN="sum"))
  
  VarKer_each_PSU <- rep(0,H)
  for(i in 1:H){
    VarKer_each_PSU[i] <- (tHT_PSU[i]-sum(d_ji_PSU[i,]*tHT_PSU))^2
  }
  VarKer_PSU <- sum(VarKer_each_PSU/Cd_PSU)/(N^2)
  VarKer_PSU <- ifelse(Cd_PSU==0,NA,VarKer_PSU)
  
  ### Hierarchical Bayesian Based Variance estimator 
  ni <- ni_PSU[1]  
  Y <- SAMPLE_PSU$y
  x <- SAMPLE_PSU$x  
  Pi <- SAMPLE_PSU$Prob
  ID <- SAMPLE_PSU$Stratum
  mcmc <- HB_PSU2(Y, x, Pi, ID, L=7, mc=10000, bn=3000)
  
  S2i_HB <- apply(mcmc$S2, 2, mean)   # strata-specific variance 
  Mu_HB <- apply(mcmc$Mu, 2, mean)
  rho_HB <- mean(mcmc$rho)
  S2i_HB_new <- array(NA,dim=c(2,2,H))
  for(h in 1:H){
    S2i_HB_new[,,h] <- matrix(c(1,rho_HB,rho_HB,1),ncol=2,nrow=2)*S2i_HB[h]
  }
  
  wij_final <- array(NA,dim=c(1,2,H))
  for(h in 1:H){
    wij_final[,,h] <- c((1/SAMPLE_PSU$Prob[2*h-1]),(1/SAMPLE_PSU$Prob[2*h]))
  }
  S2i_HB_final <- rep(NA,H)
  for(h in 1:H){
    S2i_HB_final[h] <- matrix(wij_final[,,h],1,2)%*%S2i_HB_new[,,h]%*%matrix(wij_final[,,h],2,1) 
  }
  VarHB_PSU <- sum(S2i_HB_final)/(N^2) ## final HB var
  
  # CI 
  ind <- c()
  ind[1] <- ifelse(mean_PSU-1.96*sqrt(VarColl_PSU)<overall_mean & mean_PSU+1.96*sqrt(VarColl_PSU)>overall_mean, 1, 0)
  ind[2] <- ifelse(mean_PSU-1.96*sqrt(VarKer_PSU)<overall_mean & mean_PSU+1.96*sqrt(VarKer_PSU)>overall_mean, 1, 0)
  ind[3] <- ifelse(mean_PSU-1.96*sqrt(VarHB_PSU)<overall_mean & mean_PSU+1.96*sqrt(VarHB_PSU)>overall_mean, 1, 0)
  
  cat("End rep", date(),"\n")
  return( c(VarColl_PSU,VarKer_PSU,VarHB_PSU,mean_PSU,ind) )  #Final Output of Loop
  
}

## Replicate 1000 sample: 
R <- 1000

many.reps <- replicate(n=R,one.rep()); dim(many.reps)  
many.reps <- t(many.reps)

colnames(many.reps) <- c("Collapsed","Kernel","Bayes","Mean","ind_coll","ind_ker","ind_Bayes")

## criteria of evaluations:
Bias_VarColl_PSU <- mean(abs(many.reps[,1]-true_Var))
Bias_VarKer_PSU <- mean(abs(many.reps[,2]-true_Var))
Bias_VarHB_PSU <- mean(abs(many.reps[,3]-true_Var))

RMSE_VarColl_PSU <- sqrt(mean((many.reps[,1]-true_Var)^2))
RMSE_VarKer_PSU <- sqrt(mean((many.reps[,2]-true_Var)^2))
RMSE_VarHB_PSU <- sqrt(mean((many.reps[,3]-true_Var)^2))

# boxplot
Box <- data.frame(cbind(many.reps[,1],many.reps[,2],many.reps[,3]))
colnames(Box) <- c("varColl","varKer","varHB")
boxplot(Box[,c(1,2,3)],main="Comparison")
abline(h=true_Var)

# results 
c(Bias_VarColl_PSU,Bias_VarKer_PSU,Bias_VarHB_PSU)  # Bias
c(RMSE_VarColl_PSU,RMSE_VarKer_PSU,RMSE_VarHB_PSU)  # RMSE
apply(many.reps[,5:7], 2, mean)     # CP

