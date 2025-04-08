require(MASS)
library(fMultivar)
library(copula)
library(ggplot2)
library(tidyverse)
library(extraDistr)


## Gaussian Copula


bvncop=function(n,size=1000,rho=0,funx=qnorm,funy=qnorm,...)
{
  sigma=matrix(c(1,rho,rho,1),nr=2)
  kt = NULL
  for(m in 1:size){
  data=mvrnorm(n=n,mu=c(0,0),Sigma=sigma)
  u1=pnorm(data[,1])
  u2=pnorm(data[,2])
  z1=funx(u1,...)
  z2=funy(u2)
  kt[m] = as.vector(cor.test(z1,z2,method = 'kendall')$estimate)
  }
  
  kt}
#fgm Copula

#generating data from a bivariate distribution with a specific covariancde structure
#whose marginals are exp(2)and beta(5,2)of 1st kind, using fgm Copula :

#parameters: cov = 0(i.e. under H0),
fgmcop = function(n,size,cov,...)
{ kt = NULL
  for(m in 1:size){
  myCop = fgmCopula(param = cov, dim = 2)
  myMvd = mvdc(copula = myCop, ...)	
  
  x = rMvdc(n, myMvd)
  kt[m] = as.vector(cor.test(x[,1],x[,2],method = 'kendall')$estimate)
  }
  
  kt}
  



##Distribution free: under H0 or rho = 0



## ---------------- n = 5 ----------------


## Gaussian Copula

par(mfrow = c(2,2))

hist( bvncop(n = 5,size = 1000,rho = 0,funx = qbeta,funy = qexp,shape1 = 5,shape2 = 3),breaks = 30,prob = TRUE 
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'beta(5,3) and exp(1)',xlab = 'value')
hist(bvncop(n = 5,size = 1000,rho = 0,funx = qlogis,funy = qnorm,location = 6,scale = 2),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'logistic(6,2) and n(0,1)',xlab = 'value')
hist(bvncop(n = 5,size = 1000,rho = 0,funx = qf,funy = qunif,df1=7,df2 = 4),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'F(7,4) and u(0,1)',xlab = 'value')
hist(bvncop(n = 5,size = 1000,rho = 0,funx = qunif,funy = qnorm,min = -5,max = 10),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'u(-5,10) and n(0,1)',xlab = 'value')

par(mfrow = c(1,1))

##Distribution free: under H0 or rho = 0
## fgm Copula

par(mfrow = c(2,2))

hist( fgmcop(n = 5,size=1000,cov = 0,margins = c("exp", "beta"),
             paramMargins = list(list(rate = 2), 
                                 list(shape1 = 5, shape2 = 2))),breaks = 30,prob = TRUE 
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'exp(rate = 2) and beta(5,2)',xlab = 'value')
hist(fgmcop(n = 5,size=1000,cov = 0,margins = c("logis", "unif"),
            paramMargins = list(list(location = 6,scale = 2), 
                                list(min = -10,max= 20))),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'logistic(6,2) and u(-10,20)',xlab = 'value')
hist(fgmcop(n = 5,size=1000,cov = 0,margins = c("gamma", "t"),
            paramMargins = list(list(shape = 10, rate = 2), 
                                list(df = 9))),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'gamma(10,2) and t(9)',xlab = 'value')
hist( fgmcop(n = 5,size=1000,cov = 0,margins = c("unif", "norm"),
             paramMargins = list(list(min = -5,max = 10), 
                                 list(mean = 3,sd = 3))),breaks = 30,prob = TRUE
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'u(-5,10) and n(3,9)',xlab = 'value')

par(mfrow = c(1,1))



## ---------------- n = 10 ----------------


## Gaussian Copula

par(mfrow = c(2,2))

hist( bvncop(n=10,size = 1000,rho = 0,funx = qbeta,funy = qexp,shape1 = 5,shape2 = 3),breaks = 30,prob = TRUE 
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'beta(5,3) and exp(1)',xlab = 'value')
hist(bvncop(n=10,size = 1000,rho = 0,funx = qlogis,funy = qnorm,location = 6,scale = 2),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'logistic(6,2) and n(0,1)',xlab = 'value')
hist(bvncop(n=10,size = 1000,rho = 0,funx = qf,funy = qunif,df1=7,df2 = 4),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'F(7,4) and u(0,1)',xlab = 'value')
hist(bvncop(n=10,size = 1000,rho = 0,funx = qunif,funy = qnorm,min = -5,max = 10),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'u(-5,10) and n(0,1)',xlab = 'value')

par(mfrow = c(1,1))

##Distribution free: under H0 or rho = 0

## fgm Copula

par(mfrow = c(2,2))

hist( fgmcop(n=10,size=1000,cov = 0,margins = c("exp", "beta"),
             paramMargins = list(list(rate = 2), 
                                 list(shape1 = 5, shape2 = 2))),breaks = 30,prob = TRUE 
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'exp(rate = 2) and beta(5,2)',xlab = 'value')
hist(fgmcop(n=10,size=1000,cov = 0,margins = c("logis", "unif"),
            paramMargins = list(list(location = 6,scale = 2), 
                                list(min = -10,max= 20))),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'logistic(6,2) and u(-10,20)',xlab = 'value')
hist(fgmcop(n=10,size=1000,cov = 0,margins = c("gamma", "t"),
            paramMargins = list(list(shape = 10, rate = 2), 
                                list(df = 9))),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'gamma(10,2) and t(9)',xlab = 'value')
hist( fgmcop(n=10,size=1000,cov = 0,margins = c("unif", "norm"),
             paramMargins = list(list(min = -5,max = 10), 
                                 list(mean = 3,sd = 3))),breaks = 30,prob = TRUE
     ,col = sample(colors()[-c(136:234,151:361)],1),main = 'u(-5,10) and n(3,9)',xlab = 'value')

par(mfrow = c(1,1))




## Large sample Distn: sqrt(n)tao->N(0,4/9)

sam = bvncop(n=100,size = 1000,rho = 0,funx = qbeta,funy = qexp,shape1 = 5,shape2 = 3)
hist( sqrt(900/4)*sam
      ,breaks = 30,prob = TRUE 
      ,col = sample(colors()[-c(136:234,151:361)],1),main = 'beta(5,3) and exp(1)',xlab = 'value')

curve(dnorm(x),add = T,col = sample(colors()[-c(136:234,151:361)],1),lwd = 2)

## qqPlot

qqPlot(sam)


## power comparison from gaussiancopula(unif(-5,10),N(0,1)) for different n   

vec = seq(from = -1, to = 1, length = 100)


#for n = 10


pow = NULL
i = 1

for(j in vec){
  
  pow[i] = sum(replicate(500, (abs(sqrt(10*9/4)*bvncop(n = 10,size = 1,rho = j,funx = qunif,funy = qnorm,min = -5,max = 10)) > 1.95)))/500
  i=i+1
}
plot(vec,pow,type="l",col = 5, lwd = 2,
     main = 'Power function for different sample size', ylab = 'probability',
     xlab = expression(rho),ylim = c(0,1))
abline(h = 0.05, col = 'navyblue', lty = 3, lwd = 1.2)

#for n = 25, 40, 70

l = 6
for(n in c(25, 40, 70)){
  
  pow = NULL
  i=1
  
  for(j in vec){
    
    pow[i] = sum(replicate(500, (abs(sqrt(n*9/4)*bvncop(n = n,size = 1,rho = j,funx = qunif,funy = qnorm,min = -5,max = 10)) > 1.95)))/500
    i=i+1
  }
  lines(vec, pow,col = l,ylab = '',lwd = 2)
  l = l + 1
  
}

legend('bottomright',legend = c(10,25,40,70), title = 'sample size',
       col = 5:8, lty = rep(1,4),merge = TRUE, cex = 0.8)


## violation

#by inbuild R function
ktao.pois = function(n,size = 1000){ # data = cbind(x,y)
  
  
  kt = NULL
  for(m in 1:size){
    
    
    data = rbvpois(n = n, a = 10, b = 10, c = 10)
    x = data[,1]
    y = data[,2]
    
    
    kt[m] = as.vector(cor.test(x,y,method = 'kendall')$estimate)
  }
  
  kt}

hist(ktao.pois(n = 4,size = 1000))




