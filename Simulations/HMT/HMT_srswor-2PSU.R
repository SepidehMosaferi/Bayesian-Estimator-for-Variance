# ==============================================================================
# R Code for the Simulation Studies (HMT case SRSWOR with 2 PSUs)
# Date: Feb 2024 
# ==============================================================================

# Required Libraries
require(PracTools); require(sampling); require(LaplacesDemon)
require(pscl); require(MCMCpack); require(extraDistr)

setwd("/Users/sepidehmosaferi/Documents/GitHub/Fine_Stratification_Variance/Simulations/HMT")
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

# Sample Sizes (selecting PSUs)
ni_PSU <- rep(2,H)  # vector of stratum sample sizes 

# Target Stratum Mean
Meanstratum <- tapply(DATA$y,INDEX =DATA$strat,FUN=mean)
SDstratum <- tapply(DATA$y,INDEX =DATA$strat,FUN=sd)

# Stratified Mean PSU design
overall_mean <- sum((Ni/N)*as.vector(Meanstratum))

# Variance for the purpose of bias and MSE:
true_Var <- sum((Ni/N)^2*((Ni-ni_PSU)/Ni)*(SDstratum^2/ni_PSU))

## Intro. Repetition Loop
one.rep <- function(){
  cat("start rep",date(),"\n")
  
  # Sample selection process (selecting PSUs)
  SMP.IDs_PSU <- strata(data=DATA,stratanames="strat",size=ni_PSU,method="srswor")
  SAMPLE_PSU <- getdata(DATA,SMP.IDs_PSU)  # Output of selected sample
  
  ### Collapsed Variance estimator 
  ng <- 2*ni_PSU[1]
  Ng <- rep(NA,H/2)
  for(i in 1:(H/2)){
    Ng[i] <- Ni[2*i-1]+Ni[2*i]
  }
  N <- sum(Ni)
  fg <- (Ng-ng)/Ng
  
  CollIndex <- rep(1:(H/2),each=ng)
  SAMPLE_PSU$CollIndex <- CollIndex 
  
  s2g_PSU <- rep(NA,H/2)
  for(i in 1:(H/2)){
    s2g_PSU[i] <- (sd(SAMPLE_PSU$y[SAMPLE_PSU$CollIndex==i]))^2
  }
  
  VarColl_PSU <- sum(s2g_PSU*((Ng-ng)/Ng)*(Ng^2/N^2)*(1/ng))
  
  ### H-T mean
  tHT_PSU <- rep(0,H)
  for (i in 1:H)
  {
    tHT_PSU[i] <- sum(SAMPLE_PSU$y[SAMPLE_PSU$Stratum==i]/
                        SAMPLE_PSU$Prob[SAMPLE_PSU$Stratum==i])
  }
  mean_PSU <- sum(tHT_PSU)/sum(Ni)
  
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
  VarHB_PSU <- sum(S2i_HB*((Ni-ni)/Ni)*(Ni^2/N^2)*(1/ni)) 
  
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
write.table(many.reps, sep = "\t",col.names = TRUE,row.names = FALSE,
            file = "/Users/sepidehmosaferi/Library/CloudStorage/Dropbox/Fine Stratification Variance/PSU.txt")


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

