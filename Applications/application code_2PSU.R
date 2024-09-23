## application code for the variance of fine stratification
library(dplyr); library(ggplot2) 
## 2 PSUs case:
load(file="Data_2PSU.RData")
View(Data_2PSU)

## recoding strata
Data_2PSU$SEST[Data_2PSU$SEST=="305"] <- 1
Data_2PSU$SEST[Data_2PSU$SEST=="309"] <- 2
Data_2PSU$SEST[Data_2PSU$SEST=="311"] <- 3
Data_2PSU$SEST[Data_2PSU$SEST=="312"] <- 4
Data_2PSU$SEST[Data_2PSU$SEST=="316"] <- 5
Data_2PSU$SEST[Data_2PSU$SEST=="322"] <- 6
Data_2PSU$SEST[Data_2PSU$SEST=="329"] <- 7
Data_2PSU$SEST[Data_2PSU$SEST=="330"] <- 8
Data_2PSU$SEST[Data_2PSU$SEST=="331"] <- 9
Data_2PSU$SEST[Data_2PSU$SEST=="337"] <- 10
Data_2PSU$SEST[Data_2PSU$SEST=="343"] <- 11
Data_2PSU$SEST[Data_2PSU$SEST=="347"] <- 12
Data_2PSU$SEST[Data_2PSU$SEST=="351"] <- 13
Data_2PSU$SEST[Data_2PSU$SEST=="364"] <- 14
Data_2PSU$SEST[Data_2PSU$SEST=="365"] <- 15
Data_2PSU$SEST[Data_2PSU$SEST=="366"] <- 16
Data_2PSU$SEST[Data_2PSU$SEST=="367"] <- 17
Data_2PSU$SEST[Data_2PSU$SEST=="368"] <- 18

## Y variables:
# 1. AGEPREG (pregnancy age)
# 2. EDUCAT (number of years of schooling)
# 3. NBRNALIV (Number of babies)

## we assume "x" is the poverty level income.

y <- matrix(0,nrow=18,ncol=4)
for(i in 1:18){
  for(j in 1:4){
    y[i,j] <- mean(as.numeric(Data_2PSU$NBRNALIV[Data_2PSU$SEST==i & Data_2PSU$SECU==j]),na.rm=TRUE) 
  }
}

Prob <- matrix(0,nrow=18,ncol=4)
for(i in 1:18){
  for(j in 1:4){
    Prob[i,j] <- 1/mean(as.numeric(Data_2PSU$WGT2015_2017[Data_2PSU$SEST==i & Data_2PSU$SECU==j])/1000,na.rm=TRUE)
  }
}

x <- matrix(0,nrow=18,ncol=4) 
for(i in 1:18){
  for(j in 1:4){
    x[i,j] <- mean(as.numeric(Data_2PSU$POVERTY[Data_2PSU$SEST==i & Data_2PSU$SECU==j]),na.rm=TRUE)
  }
}
x <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

Stratum <- rep(1:18,each=2)
y2 <- as.vector(t(y)) 
y2 <- y2[!is.na(y2)]
x2 <- as.vector(t(x))
x2 <- x2[!is.na(x2)]
Prob2 <- as.vector(t(Prob))
Prob2 <- Prob2[!is.na(Prob2)]

DATA2 <- as.data.frame(cbind(y2,x2,Prob2,Stratum))

### Collapsed Variance estimator
H <- 18
ng <- 4
CollIndex <- rep(1:(H/2),each=ng)
DATA2$CollIndex <- CollIndex 

tHT_PSU <- rep(0,H)
## Note: for all variables, multiply tHT_PSU[i] by 4.
for (i in 1:H)
{
  tHT_PSU[i] <- 4*sum(DATA2$y2[DATA2$Stratum==i]/
                      DATA2$Prob2[DATA2$Stratum==i])
}
VarColl_each_PSU <- rep(0,H/2)
for(i in 1:(H/2)){
  VarColl_each_PSU[i] <- 0.5*sum(tHT_PSU[2*i-1]^2,tHT_PSU[2*i]^2)
}
N <- 2149 ##total number of PSUs
VarColl_PSU <- sum(VarColl_each_PSU)/(N^2)

### H-T mean 
mean_PSU <- sum(tHT_PSU)/N 

### Kernel-Based Variance estimator
h <- (1/H+2/H)/2

Kernel_PSU <- function(x,xi,h){u<-(x-xi)/h; k<-0.75*(1-u^2); k*(k>0)}

x_ker <- sapply(1:H,function(i){sum(DATA2$x2[DATA2$Stratum==i])}) 

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

### Hierarchical Bayesian Based Variance estimator 
source("HB-var-2PSU-HMT.R")
Y <- DATA2$y2
x <- DATA2$x2  
Pi <- DATA2$Prob2
ID <- DATA2$Stratum
mcmc <- HB_PSU2(Y, x, Pi, ID, L=7, mc=10000, bn=3000)

S2i_HB <- apply(mcmc$S2, 2, mean)   # strata-specific variance 
Mu_HB <- apply(mcmc$Mu, 2, mean)
wij_final <- rep(0,H)
for(h in 1:H){
  wij_final[h] <- (1/DATA2$Prob2[2*h-1]^2)+(1/DATA2$Prob2[2*h]^2)
}
S2i_HB <- wij_final*S2i_HB 
VarHB_PSU <- sum(S2i_HB)/(N^2) ## final HB var

low_coll <- mean_PSU-1.96*sqrt(VarColl_PSU)
high_coll <- mean_PSU+1.96*sqrt(VarColl_PSU) 
low_ker <- mean_PSU-1.96*sqrt(VarKer_PSU)
high_ker <- mean_PSU+1.96*sqrt(VarKer_PSU)
low_HB <- mean_PSU-1.96*sqrt(VarHB_PSU)
high_HB <- mean_PSU+1.96*sqrt(VarHB_PSU) 
lowBand <- c(low_coll,low_ker,low_HB)
highBand <- c(high_coll,high_ker,high_HB) 
xAxis <- c(1,2,3)
yAxis <- rep(mean_PSU,3)
method=c("Collpased","Kernel-based","HB")
sample_data <- data.frame(method,xAxis, yAxis, lowBand, highBand) 

ggplot(sample_data, aes(xAxis, yAxis,colour=method)) + geom_point() +  
  geom_errorbar(aes(ymin = lowBand, ymax = highBand))+theme_bw()+
  labs(y="Number of Babies",x="Variance Method Estimators")+
  ggtitle("Confidence Intervals for Number of Babies (2 PSUs case)")+
  theme(legend.position="none")+
  xlim("Collpased Estimator", "Kernel-Based Estimator","Bayesian Estimator")


