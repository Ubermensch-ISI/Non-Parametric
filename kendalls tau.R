##libraries

require(car) #for qqPlot()
require(kableExtra)
require(MASS) #for generaing bivariate normal sample
library(fMultivar)
library(copula)
library(tidyverse)

## Code for kendall's tau



ktao = function(n,size = 1000,dist1,dist2){
  
  
  kt = NULL
  for(m in 1:size){
    
    x = dist1(n)
    y = dist2(n)
    
    k = 1
    store = NULL
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        store[k] = sign(x[i]-x[j])*sign(y[i]-y[j])
        k = k+1
      }
    }
  kt[m] = sum(store)/choose(n,2)
    }
    
kt}

#by inbuild R function
ktao = function(n,size = 1000,dist1,dist2){
  
  
  kt = NULL
  for(m in 1:size){
    
    x = dist1(n)
    y = dist2(n)
    
    
    kt[m] = as.vector(cor.test(x,y,method = 'kendall')$estimate)
  }
  
  kt}
   

## for different rho

par(mfrow=c(1,3))
for(r in c(-0.7,0,0.7)){
  data=mvrnorm(n=1000,mu=c(0,0),Sigma = cbind(c(1,r),c(r,1)))
  plot(data[,1],data[,2],col = ifelse(data[,1]*data[,2]>0,'navyblue','tomato')
       ,ylab = '',xlab = paste(expression(rho),'=',r),lwd = 2)
  grid()
  abline(h = 0,v = 0,col = 'green', lwd = 2)
  
  
}


par(mfrow = c(1,1))


##Mean and variance under $H_0$

sam = ktao(n = 10,size = 10000,rnorm,runif)


out = matrix(c(0,(2*(2*10+5)/(9*10*(10-1))),mean(sam),var(sam)),ncol = 2,byrow = F)
rownames(out) = c('mean','variance')
colnames(out) = c('Theoritical','Observed')
kable(out,caption = 'n = 100')







##Distribution free: EDA

par(mfrow = c(2,2))

hist(ktao(n = 10,size = 1000,rnorm,runif),prob = TRUE ,col = 'sienna2')
hist(ktao(n = 10,size = 1000, rnorm,rcauchy),prob = TRUE,col = 'sienna2' )
hist(ktao(n = 10,size = 1000, rexp,runif),prob = TRUE,col = 'sienna2' )
hist(ktao(n = 10,size = 1000, rnorm,rexp),prob = TRUE,col = 'sienna2' )

par(mfrow = c(1,1))

##Distribution free: Testing

chisq.test(ktao(n = 10,size = 1000,rnorm,runif),ktao(n = 10,size = 1000, rnorm,rcauchy))





##Large sample test: sqrt(n-1)*rho => N(0,1) as n -> inf

n = 100
l.sam = sqrt((9*n*(n-1))/(2*(2*n+5)))*ktao(n = 100,size = 1e3,rnorm,rexp)

#Checking for Normality: EDA

hist(l.sam, prob = TRUE,col = 'cadetblue1')
curve(dnorm(x),add = TRUE, col = 'darkorchid')

#comment: Good fit

qqPlot(l.sam, ylab = 'Theoritical Quantiles', pch = 20)


#Testing for Normality:  ks test 



ks.test( l.sam, "pnorm")

#Large sample Power curve for BVN(rho)




##new variation of ktao: it takes two sample x and y with size n


 

ktao.bvn = function(s.size = 1000,r){
    
    
      n = s.size
      sam = mvrnorm(n = s.size, mu = c(0,0), 
                    Sigma = matrix(c(1, r, r, 1), 
                                   ncol = 2))
      
      x = sam[,1]
      y = sam[,2]
      k = 1
      store = NULL
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          store[k] = sign(x[i]-x[j])*sign(y[i]-y[j])
          k = k+1
        }
      }
  sum(store)/choose(n,2)
    }
 
  
vec = seq(from = -1, to = 1, length = 100)

#for n = 10

pow = NULL
i = 1

for(j in vec){
  
  pow[i] = sum(replicate(500, (abs(sqrt((9*10*(10-1))/(2*(2*10+5)))*ktao.bvn(s.size = 10,r = j)) > 1.95)))/500
  i=i+1
}
plot(vec,pow,type="l",col = 5, lwd = 2,
     main = 'Power function for different sample size', ylab = 'probability',
     xlab = expression(rho),ylim = c(0,1))
abline(h = 0.05, col = 'navyblue', lty = 3, lwd = 1.2)

#for n = 50, 100, 200

l = 6
for(n in c(50,100,200)){
  
  pow = NULL
  i=1
  
  for(j in vec){
    
    pow[i] = sum(replicate(100, (abs(sqrt((9*n*(n-1))/(2*(2*n+5)))*ktao.bvn(s.size = n,r = j)) > 1.95)))/100
    i=i+1
  }
  lines(vec, pow,col = l,ylab = '',lwd = 2)
  l = l + 1
  
}

legend('bottomright',legend = c(10,50,100,500), title = 'sample size',
       col = 5:8, lty = rep(1,4),merge = TRUE, cex = 0.8)

## Power curves for two sided t-test in BVN with n: 20, 50, 100, 200

# n: 20

vec = seq(from = -1, to = 1, length = 100)
pow = NULL
i = 1
for(j in vec){

  pow[i] = sum(replicate(500,{sam = mvrnorm(n = 20, mu = c(0,0),
                                            Sigma = matrix(c(1, j, j, 1),
                                                           ncol = 2))
  abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = 18, lower.tail = F)}  ))/500
  i=i+1
}

plot(vec, pow,col = 6,ylab = 'Probability',lwd = 2, type = 'l')
abline(h = 0.05, col = 'grey', lwd = 2)

# n: 50,100,200

for(s in c(50,100,200)){
  vec = seq(from = -1, to = 1, length = 100)
  pow = NULL
  i = 1
  for(j in vec){

    pow[i] = sum(replicate(500,{sam = mvrnorm(n = s, mu = c(0,0),
                                              Sigma = matrix(c(1, j, j, 1),
                                                             ncol = 2))
    abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = (s-2), lower.tail = F)}  ))/500
    i=i+1
  }

  lines(vec, pow,col = s,ylab = '',lwd = 2, type = 'l')

}


legend('bottomright',legend = c(20,50,100,200), title = 'sample size',
       col = c(6,50,100,200), lty = rep(1,4), cex = 0.8)

## Power curves for two sided t-test in BVN with n: 20, 50, 100, 200


# n: 20

vec = seq(from = -1, to = 1, length = 100)
pow = NULL
i = 1
for(j in vec){
  
  pow[i] = sum(replicate(500,{sam = mvrnorm(n = 20, mu = c(0,0),
                                            Sigma = matrix(c(1, j, j, 1),
                                                           ncol = 2))
  abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = 18, lower.tail = F)}  ))/500
  i=i+1
}

plot(vec, pow,col = 6,ylab = 'Probability',lwd = 2, type = 'l')
abline(h = 0.05, col = 'grey', lwd = 2)

# n: 50,100,200

for(s in c(50,100,200)){
  vec = seq(from = -1, to = 1, length = 100)
  pow = NULL
  i = 1
  for(j in vec){
    
    pow[i] = sum(replicate(500,{sam = mvrnorm(n = s, mu = c(0,0),
                                              Sigma = matrix(c(1, j, j, 1),
                                                             ncol = 2))
    abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = (s-2), lower.tail = F)}  ))/500
    i=i+1
  }
  
  lines(vec, pow,col = s,ylab = '',lwd = 2, type = 'l')
  
}


legend('bottomright',legend = c(20,50,100,200), title = 'sample size',
       col = c(6,50,100,200), lty = rep(1,4), cex = 0.8)



## comparison with parametric test

# kendalls tau and t-test for BVN

vec = seq(from = -1, to = 1, length = 100)

#for n = 20

pow1 = NULL
pow2 = NULL
i = 1

for(j in vec){

  pow1[i] = sum(replicate(500, (abs(sqrt((9*20*(20-1))/(2*(2*20+5)))*ktao.bvn(s.size = 20,r = j)) > 1.95)))/500
  pow2[i] = sum(replicate(500,{sam = mvrnorm(n = 20, mu = c(0,0),
                                            Sigma = matrix(c(1, j, j, 1),
                                                           ncol = 2))
  abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = 18, lower.tail = F)}  ))/500



  i=i+1
}
plot(vec,pow1,type="l",col = 5, lwd = 2,
     main = 'Power function for kendalls tau and parametric test ', ylab = 'probability',
     xlab = expression(rho),ylim = c(0,1))
abline(h = 0.05, col = 'navyblue', lty = 3, lwd = 1.2)
lines(vec, pow2,col = 6,ylab = '',lwd = 2)

legend('top',col=c(5,6),lwd = rep(2,2),legend = c('kendalls tau','parametric test'), cex = 0.8)


# comment: t-test perform better

## when sample is not comming from Normal Distribution, kendalls tau will perform better

vec = seq(from = -1, to = 1, length = 100)

#for n = 20

pow1 = NULL
pow2 = NULL
i = 1

for(j in vec){

  pow1[i] = sum(replicate(500, {sam = rt2d(20,rho= j,nu=4)
    abs(sqrt((9*20*(20-1))/(2*(2*20+5)))*as.vector(cor.test(sam[,1],sam[,2],
                                                          method = 'kendall')$estimate)) > 1.95}))/500
  pow2[i] = sum(replicate(500,{sam = rt2d(20,rho= j,nu=4)
  abs(cor.test(sam[,1],sam[,2])$statistic ) > qt(0.025,df = 18, lower.tail = F)}  ))/500



  i=i+1
}
plot(vec,pow1,type="l",col = 5, lwd = 2,
     main = 'Power function for kendalls tau and parametric test ', ylab = 'probability',
     xlab = expression(rho),ylim = c(0,1))
abline(h = 0.05, col = 'navyblue', lty = 3, lwd = 1.2)
lines(vec, pow2,col = 6,ylab = '',lwd = 2)

legend('top',col=c(5,6),lwd = rep(2,2),legend = c('kendalls tau','parametric test'), cex = 0.8)

# comment: kendalls tau performs better







## spearman vs kendall's tau

vec = seq(from = -1, to = 1, length = 100)
#for n = 10

pow = NULL
i = 1
for(j in vec){
  
  pow[i] = sum(replicate(500, (abs(sqrt((9*30*(30-1))/(2*(2*30+5)))*ktao.bvn(s.size = 30,r = j)) > 1.95)))/500
  i=i+1
}
plot(vec,pow,type="l",col = 5, lwd = 2,
     main = '', ylab = 'probability',
     xlab = expression(rho),ylim = c(0,1))
abline(h = 0.05, col = 'navyblue', lty = 3, lwd = 1.2)

#for spearman

pow = NULL
i = 1

for(j in vec){
  
  pow[i] = sum(replicate(500, (abs(sqrt(30 -1)*sprho.bvn(s.size = 30,r = j)) > 1.95)))/500
  i=i+1
}
lines(vec,pow,type="l",col = 6, lwd = 2)
legend('top',col=c(6,5),lwd = rep(2,2),legend = c('spearman rho','kendalls tao'))




                            
   
   