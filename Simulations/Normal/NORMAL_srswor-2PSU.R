# ==============================================================================
# R Code for the Simulation Studies (Normal case SRSWOR with 2 PSUs)
# Date: Feb 2024 
# ==============================================================================

# Required Libraries
require(PracTools); require(sampling); require(LaplacesDemon); require(MASS)
require(pscl); require(MCMCpack); require(extraDistr); require(statmod)

setwd("/Users/sepidehmosaferi/Documents/GitHub/Bayesian-Estimator-for-Variance/Simulations/Normal")
rm(list = ls(all = TRUE))
source("HB-var-2PSU.R")

# values of parameters:
# H={50,100,200}; 1/H < h < 2/H; sigma={0.25,0.5}

# settings
Ni <- 60  # Size of each stratum
H <- 50  # Number of stratum
h <- 0.03  # Bandwidth
i <- 1:H
x_i <- i/H


sigma <- 5

m_starx <- 1+2*(x_i-0.5)
m_starx <- 2*(m_starx-min(m_starx))/(max(m_starx)-min(m_starx)) # Rescaled regression function
plot(x_i, m_starx)

# data generation 
set.seed(2012)
y <- matrix(0,ncol=H,nrow=Ni) # each col represents a stratum
for(i in 1:H){
  y[,i] <- m_starx[i]+rnorm(Ni,0,sigma)
}


# Creating data frame
STRATUM <- rep(1:H,each=Ni)
X <- rep(x_i,each=Ni) # Collapsing Index
Y <- as.vector(y)
DATAframe <- data.frame(Y,STRATUM,X)   # full dataset


# Sample Sizes (selecting PSUs {1,2,3,4,5})
ni_PSU <- rep(2,H)  # vector of stratum sample sizes 

# Target Stratum Mean
Meanstratum <- tapply(DATAframe$Y,INDEX =DATAframe$STRATUM,FUN=mean)
SDstratum <- tapply(DATAframe$Y,INDEX =DATAframe$STRATUM,FUN=sd)

# Stratified Mean PSU design
overall_mean <- sum((Ni/(Ni*H))*as.vector(Meanstratum))   # population mean under design

# Variance for the purpose of bias and MSE:
true_Var <- sum((as.vector(table(DATAframe$STRATUM)/sum(table(DATAframe$STRATUM))))^2
                *((as.vector(table(DATAframe$STRATUM))-ni_PSU)/as.vector(table(DATAframe$STRATUM)))
                *(SDstratum^2/ni_PSU))

## Intro. Repetition Loop
one.rep <- function(){
  cat("start rep",date(),"\n")
  
  # Sample selection process (selecting PSUs)
  SMP.IDs_PSU <- strata(data=DATAframe,stratanames="STRATUM",size=ni_PSU,method="srswor")
  SAMPLE_PSU <- getdata(DATAframe,SMP.IDs_PSU)  # Output of selected sample
  
  ## H-T mean
  sum_PSU <- rep(NA,H)
  for(i in 1:H){
    sum_PSU[i] <- sum(SAMPLE_PSU$Y[SAMPLE_PSU$STRATUM==i]/SAMPLE_PSU$Prob[SAMPLE_PSU$STRATUM==i])
  }
  mean_PSU <- sum(sum_PSU)/(Ni*H)
  
  ### Collapsed Variance estimator 
  ng <- 2*ni_PSU[1]
  Ng <- 2*Ni
  N <- H*Ni
  fg <- (Ng-ng)/Ng
  
  CollIndex <- rep(1:(H/2),each=ng)
  SAMPLE_PSU$CollIndex <- CollIndex 
  
  s2g_PSU <- rep(NA,H/2)
  for(i in 1:(H/2)){
    s2g_PSU[i] <- (sd(SAMPLE_PSU$Y[SAMPLE_PSU$CollIndex==i]))^2
  }
  
  VarColl_PSU <- sum(s2g_PSU)*((Ng-ng)/Ng)*(Ng^2/ng)/(N^2)
  
  ### Kernel-Based Variance estimator
  tHT_PSU <- rep(0,H)
  for (i in 1:H)
  {
    tHT_PSU[i] <- sum(SAMPLE_PSU$Y[SAMPLE_PSU$STRATUM==i]/
                         SAMPLE_PSU$Prob[SAMPLE_PSU$STRATUM==i])
  }
  
  Kernel_PSU <- function(x,xi,h){u<-(x-xi)/h; k<-0.75*(1-u^2); k*(k>0)}
  
  d_ji_PSU <- matrix(0,nrow=H,ncol=H) #column is stratum
  for(i in 1:H){
    d_ji_PSU[i,] <- Kernel_PSU(x_i,x_i[i],h)
    d_ji_PSU[i,] <- d_ji_PSU[i,]/sum(d_ji_PSU[i,]) 
  }
  
  Cd_PSU <-mean(rep(1,H)-2*diag(d_ji_PSU)+apply(d_ji_PSU^2,MAR=1,FUN="sum"))
  
  VarKer_each_PSU <- rep(0,H)
  for(i in 1:H){
    VarKer_each_PSU[i] <- (tHT_PSU[i]-sum(d_ji_PSU[i,]*tHT_PSU))^2
  }
  VarKer_PSU <- (sum(VarKer_each_PSU)/Cd_PSU)/(N^2)
  
  ### Hierarchical Bayesian Based Variance estimator 
  ni <- ni_PSU[1]  
  Y <- SAMPLE_PSU$Y   # strata-specific mean 
  x <- SAMPLE_PSU$X
  Pi <- SAMPLE_PSU$Prob
  ID <- SAMPLE_PSU$STRATUM
  mcmc <- HB_PSU2(Y, x, Pi, ID, p=2,L=7, mc=10000, bn=3000)
  
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
  VarHB_PSU <- sum(S2i_HB_final*((Ni-ni)/Ni)*(Ni^2/N^2)*(1/ni)) 
  
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



