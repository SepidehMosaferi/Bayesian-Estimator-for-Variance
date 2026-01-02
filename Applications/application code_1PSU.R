## application code for the variance of fine stratification
library(dplyr); library(ggplot2) 
## 1 PSU case:
load(file="Data_1PSU.RData")
View(Data_1PSU)

## recoding strata
Data_1PSU$SEST[Data_1PSU$SEST=="305"] <- 1
Data_1PSU$SEST[Data_1PSU$SEST=="309"] <- 2
Data_1PSU$SEST[Data_1PSU$SEST=="311"] <- 3
Data_1PSU$SEST[Data_1PSU$SEST=="312"] <- 4
Data_1PSU$SEST[Data_1PSU$SEST=="316"] <- 5
Data_1PSU$SEST[Data_1PSU$SEST=="322"] <- 6
Data_1PSU$SEST[Data_1PSU$SEST=="329"] <- 7
Data_1PSU$SEST[Data_1PSU$SEST=="330"] <- 8
Data_1PSU$SEST[Data_1PSU$SEST=="331"] <- 9
Data_1PSU$SEST[Data_1PSU$SEST=="337"] <- 10
Data_1PSU$SEST[Data_1PSU$SEST=="343"] <- 11
Data_1PSU$SEST[Data_1PSU$SEST=="347"] <- 12
Data_1PSU$SEST[Data_1PSU$SEST=="351"] <- 13
Data_1PSU$SEST[Data_1PSU$SEST=="364"] <- 14
Data_1PSU$SEST[Data_1PSU$SEST=="365"] <- 15
Data_1PSU$SEST[Data_1PSU$SEST=="366"] <- 16
Data_1PSU$SEST[Data_1PSU$SEST=="367"] <- 17
Data_1PSU$SEST[Data_1PSU$SEST=="368"] <- 18

## Y variables:
# 1. AGEPREG (pregnancy age)
# 2. EDUCAT (number of years of schooling)
# 3. NBRNALIV (Number of babies)

## we assume "x" is the poverty level income.

y  <- sapply(1:18,function(i){mean(as.numeric(Data_1PSU$NBRNALIV[Data_1PSU$SEST==i]),na.rm=TRUE)})
Prob <- 1/sapply(1:18,function(i){mean(as.numeric(Data_1PSU$WGT2015_2017[Data_1PSU$SEST==i])/1000,na.rm=TRUE)})
x <- sapply(1:18,function(i){mean(as.numeric(Data_1PSU$POVERTY[Data_1PSU$SEST==i]),na.rm=TRUE)}) 
x <- (x-mean(x))/sd(x)
Stratum <- 1:18

DATA <- as.data.frame(cbind(y,x,Prob,Stratum))

### Collapsed Variance estimator
H <- 18
ng <- 2
CollIndex <- rep(1:(H/2),each=ng)
DATA$CollIndex <- CollIndex 

tHT_PSU <- rep(0,H)

for (i in 1:H)
{
  tHT_PSU[i] <- sum(DATA$y[DATA$Stratum==i]/
                      DATA$Prob[DATA$Stratum==i])
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

x_ker <- DATA$x 

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
source("HB-var-1PSU-HMT.R")
Y <- DATA$y
x <- DATA$x  
Pi <- DATA$Prob
mcmc <- HB_var(Y, x, Pi, L=5, mc=10000, bn=3000)

S2i_HB <- apply(mcmc$S2, 2, mean)   # strata-specific variance 
Mu_HB <- apply(mcmc$Mu, 2, mean)
wi <- 1/(DATA$Prob^2) # weight
S2i_HB <- wi*S2i_HB
VarHB_PSU <- sum(S2i_HB)/(N^2) ## HB final var

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
  ggtitle("Confidence Intervals for Number of Babies (1 PSU case)")+
  theme(legend.position="none")+
  xlim("Collpased Estimator", "Kernel-Based Estimator","Bayesian Estimator")





