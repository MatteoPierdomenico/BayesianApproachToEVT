## STAN version of the Hierarchical random effects model for daily maximum wind gusts
## varying only by sites.

library(ggplot2)
library(POT)
library(lubridate)
library(evir)
library(coda)
library(nimble)


Ireland <- read.csv("../../../Dataset/Ireland_daily_from1990.csv")
head(Ireland)
unique(Ireland$spot)

df <- data.frame(time=Ireland$year, obs=Ireland$hg, county=Ireland$spot)
n <- length(df$obs)
m <- length(unique(df$county))

spots <- list()
data <- list()
counties <- as.vector(unique(Ireland$spot))
for(i in 1:m){ 
  county <- Ireland[which(Ireland$spot==counties[i]),]
  spots[[i]] <- data.frame(time=county$year, obs=county$hg) # Raccolgo in una lista le osservazioni divise per siti
  data[[i]] <- county$hg
  plot(pacf(spots[[i]]$obs))
  readline(prompt = "press [enter] to continue or ESC to stop")
  graphics.off()
}


thresholds <- c(70,70,68,68,63)  

## The thresholds above, for every site, have been empirically estimated 
## by the following commented lines, which represent the empirical mean residual life plot (mrlplot) 
## and the Threshold Choice Plot (tcplot)


##i=5    # i varie from 1 to 5 for the sites (respectively Clare, Cork, Dublin, Kerry and Westmeath)

# initial data selection
#initial_threshold <- 51
#events0 <- c(spots[[i]][which(spots[[i]]$obs>initial_threshold),])
#events0 <- data.frame(time = events0$time,obs= events0$obs, idx=which(spots[[i]]$obs>initial_threshold))
## threshold selection (try 12)
#par(mfrow = c(2, 2))
#mrlplot(events0[, "obs"]) # should be linear after threshold
#abline(v = initial_threshold, col = "green")
#diplot(events0) # should be close to 1
#abline(v = initial_threshold, col = "green")
#tcplot(events0[, "obs"], which = 1) # should be constant after threshold
#abline(v = initial_threshold, col = "green")
#tcplot(events0[, "obs"], which = 2) # should be constant after threshold
#abline(v = initial_threshold, col = "green")
#graphics.off()

#thresholds[i] <- initial_threshold
#lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>initial_threshold),2])/length(spots[[i]]$obs)

#tim.cond <- 3/365
#
## cluster to avoid short range temporal dependence
#for(i in 1:m){
#spots[[i]] <- clust(spots[[i]], u = thresholds[i], tim.cond = tim.cond, 
#                 clust.max = TRUE, plot = TRUE)
#}

## Computing lambda for every choice of threshold. Lambda is the normalized number of
## peaks over threshold for each site.

lambda <- rep(NA,m)

for(i in 1:m){
  lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>thresholds[i]),2])/length(spots[[i]]$obs)
}


## Collecting data in a matrix in which every column represent a different site

N<-length(data[[1]])
datamatrix<-matrix(NA,N,m)
for(i in 1:N){
  for(j in 1:m){
    datamatrix[i,j]<-spots[[j]][i,2]
  }
}

#FOMC_logit_GPD_joint_density <- nimbleFunction(
#   run = function(x1=double(0),x2=double(0),threshold=double(0),sigma=double(0),xi=double(0),alpha=double(0), lambda=double(0)){
#     returnType(double(0))
#     z1 = 1 + xi*(x1 - threshold)*(1/sigma)
#     z2 = 1 + xi*(x2- threshold)*(1/sigma)
#     inv_xi = 1/xi
#     inv_alpha = 1/alpha 
#     if(x1>threshold & x2>threshold){
#         return(-(z1^inv_xi*1/lambda)^(-inv_alpha)*(z2^inv_xi*1/lambda)^(-inv_alpha)*((z2^inv_xi*1/lambda)^(-inv_alpha) + (z1^inv_xi*1/lambda)^(-inv_alpha))^alpha/(sigma^2*(z1)*(z2)*(((z2)^(1/xi)*1/lambda)^(-inv_alpha) + ((z1)^(inv_xi)*1/lambda)^(-inv_alpha))^2) + ((z1)^(inv_xi)*1/lambda)^(-inv_alpha)*((z2)^(inv_xi)*1/lambda)^(-inv_alpha)*(((z2)^(inv_xi)*1/lambda)^(-inv_alpha) + ((z1)^(inv_xi)*1/lambda)^(-inv_alpha))^alpha/(alpha*sigma^2*(z1)*(z2)*(((z2)^(inv_xi)*1/lambda)^(-inv_alpha) + ((z1)^(inv_xi)*1/lambda)^(-inv_alpha))^2))
#     }
#     if(x1>threshold & x2<=threshold){
#           return(((z1)^(inv_xi)*1/lambda)^(-inv_alpha)*(1/lambda^(-inv_alpha) + ((z1)^(1/xi)*1/lambda)^(-inv_alpha))^alpha/(sigma*(z1)*(1/lambda^(-inv_alpha) + ((z1)^(inv_xi)*1/lambda)^(-inv_alpha))))
#       }
#     if(x1<=threshold & x2>threshold){
#         return(((z2)^(inv_xi)*1/lambda)^(-inv_alpha)*(1/lambda^(-inv_alpha) + ((z2)^(1/xi)*1/lambda)^(-inv_alpha))^alpha/(sigma*(z2)*(1/lambda^(-inv_alpha) + ((z2)^(inv_xi)*1/lambda)^(-inv_alpha))))
#     }
#     else{    
#            return(1-((1/lambda^(-inv_alpha) + 1/lambda^(-inv_alpha))^(alpha)))
#       }
#})


## Modelling First Order Markov Chain dependence

FOMC_logit_GPD_joint_density <- nimbleFunction(
  run = function(x1=double(0),x2=double(0),threshold=double(0),sigma=double(0),xi=double(0),alpha=double(0), lambda=double(0)){
  returnType(double(0))
    if(x1>threshold & x2>threshold){
      return(-((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma^2*(1 + xi*(-threshold + x1)/sigma)*(1 + xi*(-threshold + x2)/sigma)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^2) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(alpha*sigma^2*(1 + xi*(-threshold + x1)/sigma)*(1 + xi*(-threshold + x2)/sigma)*(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^2));
    }
    
    if(x1>threshold & x2<=threshold){
      return(((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma*(1 + xi*(-threshold + x1)/sigma)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x1)/sigma)^(1/xi)/lambda)^(-1/alpha))));
    }
    if(x1<=threshold & x2>threshold){
      return(((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha)*((1 /lambda)^(-1/alpha) + ((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha))^alpha/(sigma*(1 + xi*(-threshold + x2)/sigma)*((1/lambda)^(-1/alpha) + ((1 + xi*(-threshold + x2)/sigma)^(1/xi)/lambda)^(-1/alpha))));
    }
    else{    
      return(1-(((1/lambda)^(-1/alpha) + (1/lambda)^(-1/alpha))^(alpha)));
    }
  })

## Generalized Pareto Density:

GPD_density <- nimbleFunction(
  run=function(x=double(0),threshold=double(0),scale=double(0), xi=double(0), lambda=double(0)){
  returnType(double(0))
  if(x >= threshold){
    if(abs(xi) <= 1e-15)
     {return( lambda/scale*exp(-(x-threshold)/scale))}
    else
     {return((lambda/scale) * (1 + (xi * (x - threshold))/scale)^((-1/xi) - 1))}  
  }
  else{
    return((1-lambda)/threshold)
  }
  
})

## Likelihood model

dmyModel_lpost <- nimbleFunction(
  run=function(x=double(1),threshold=double(0),scale=double(0),xi=double(0), lambda=double(0), alpha=double(0),log = integer(0, default = 0)){
  returnType(double(0))
  loglikelihood <- log(GPD_density(x[1], threshold,scale,xi,lambda))
  num <-0
  denom <-0
  for(i in 1:(length(x)-1)){
    num <- log(FOMC_logit_GPD_joint_density(x[i],x[i+1],  threshold,scale,xi, alpha, lambda)) +num
    denom <- log(GPD_density(x[i],  threshold,scale,xi, lambda)) + denom
  }
  loglikelihood <- loglikelihood + num - denom
  if(log) return(loglikelihood)
  else return(exp(loglikelihood))
}
)

code <- nimbleCode({
  for(i in 1:M){
    y[1:N,i] ~ dmyModel_lpost(threshold[i],sigma[i],xi[i],lambda[i],alpha[i])}
  for(i in 1:M){
    constraint_data[i] ~ dconstraint( xi[i] > -sigma[i]/(max(y[1:N,i])-threshold[i]))
    xi[i] ~ dnorm(a_xi,phi_xi)
    sigma[i] <- exp(logsigma[i])
    logsigma[i] ~ dnorm(a_sigma,phi_sigma)
    alpha[i] ~ dunif(0,1)}
  a_sigma ~ dnorm(b,c)
  a_xi ~ dnorm(b,c)
#  phi_sigma ~ dgamma(f,f)
#  phi_xi ~ dgamma(f,f)
  phi_sigma ~ dunif(0,1)
  phi_xi ~ dunif(0,1)
})


pumpConsts <- list(M= m, N=N,threshold = thresholds, lambda = lambda, b=0, c=10^-2,f=5*10^-2)
pumpData <-list(y = datamatrix)


pumpInits <- list()

## Se si considerano delle gamma per phi_sigma e phi_xi, 
## per ottenere conjugacy, decommentare il comando che segue:

#pumpInits<-list(phi_sigma=1, phi_xi=1)  

Hmodel <- nimbleModel(code=code, name ="Hmodel", constants = pumpConsts, data=pumpData, inits = pumpInits)

configureMCMC(Hmodel)

mcmc.out <- nimbleMCMC(model = Hmodel,
                       niter = 8000, nchains = 3, thin = 10, nburnin = 2000, 
                       monitors = c("a_sigma","a_xi","phi_sigma","phi_xi", "logsigma", "xi", "alpha", "sigma"),
                       summary = TRUE, WAIC = TRUE, samplesAsCodaMCMC  = TRUE, setSeed = TRUE)


coda_chain<- mcmc.out$samples
summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
x11()
plot(coda_chain)
acfplot(coda_chain)
cumuplot(coda_chain)
graphics.off()

geweke.diag(coda_chain)
geweke.plot(coda_chain, frac1 = 0.1, frac2 = 0.5, nbins = 20)

###
#CREATING A LIST OF THE MAIN RESULTS OF THE CHAIN
###

## Change the first four lines to analyze different sites,
## in the following, for example, we are chosing i=4 which is Kerry

sigmacb<-coda_chain$chain2[,"sigma[4]"]
xicb <- coda_chain$chain2[,"xi[4]"]
obs <- datamatrix[,4]
threshold<-thresholds[4]

estim = cbind( sigmacb, xicb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]))
ests2 <- array(0, c(2, 2))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), quantile(estim[, 2], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), quantile(estim[, 2], 0.975))
results <- list(posterior = estim, data = obs, postmean = ests, 
                postmedian = ests1, postCI = ests2, block = 100)
names(results$postmean) <- c("sigma", "xi")
names(results$postmedian) <- c( "sigma", "xi")
dimnames(results$postCI) <- list(c("lower bound", "upper bound"), 
                                 c( "sigma", "xi"))
class(results) <- "gpdp"

###
#RETURN LEVEL COMPUTATION AND PLOT
###


t <- 2  #start ret level
k <- 100   #end ret level
x <- results

sampl <- qgpd(1 - 1/t, x$posterior[, 2], threshold, x$posterior[, 1])
res <- quantile(sampl, 0.5)
ta <- seq(1, k, 1)
n <- length(ta)
li <- array(0, c(n))
ls <- array(0, c(n))
pred <- array(0, c(n))
for (s in 1:n) {
  sampl <- qgpd(1 - 1/s, x$posterior[, 2], threshold, x$posterior[, 1])
  li[s] <- quantile(sampl, 0.025)
  ls[s] <- quantile(sampl, 0.975)
  pred[s] <- quantile(sampl, 0.5)
}
x11()
ggplot(data = data.frame( level= ta, intensity=pred), mapping = aes(x=level,y=intensity)) +
  geom_line( color = "red") +
  xlab("level") + 
  ylab("wind speed (kmh)") + 
  ggtitle("Return level") +
  xlim(0,100)+
  ylim(80,200) +
  geom_line(aes(y = li), color = "black", linetype = "dotted") +
  geom_line(aes(y = ls), color = "black", linetype = "dotted")
graphics.off()



##
#log predictive density
##

data <- x$data[x$data > threshold]
linf = min(data)
lsup = max(data)
dat1 = seq(linf, lsup, (lsup - linf)/300)
n = length(dat1)
int = length(x$posterior[, 1])
res = array(0, c(n))
for (i in 1:n) {
  for (j in 1:int) {
    res[i] = res[i] + (1/int) * dgpd(dat1[i], x$posterior[j, 2], threshold, x$posterior[j, 1])
  }
}
dataf <- data.frame(log_pos = data)
dataf2 <- data.frame(data = dat1, res = res)

ggplot(dataf, aes(x=log_pos)) + 
  geom_histogram(aes(y = ..density.. ), bins = 20, color = "black", fill = "lightblue") +
  xlab("data") + 
  ylab("density") +
  ggtitle("density") +
  geom_line(aes(x=data,y = res), dataf2,color = "red") 

