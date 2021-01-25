library(ggplot2)
library(POT)
library(lubridate)
library(evir)
library(coda)
library(nimble)


Ireland <- read.csv('C:/Users/matte/OneDrive/Documenti/GitHub/Bayesian-Approach-to-EVT/Dataset/Ireland_daily_from1990.csv')
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
  spots[[i]] <- data.frame(time=county$year, obs=county$hg)
  data[[i]] <- county$hg
  plot(pacf(spots[[i]]$obs))
  readline(prompt = "press [enter] to continue or ESC to stop")
  graphics.off()
}


thresholds <- rep(NA,m)
lambda <- rep(NA,m)
thresholds <- c(95,95,94,94,85)  #sono state scelte empiricamente
for(i in 1:m){
  lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>thresholds[i]),2])/length(spots[[i]]$obs)
}




# i=5 #loop manuale (modificare i fino a 5)
# 
# # initial data selection
# initial_threshold <- 85
# events0 <- c(spots[[i]][which(spots[[i]]$obs>initial_threshold),])
# events0 <- data.frame(time = events0$time,obs= events0$obs, idx=which(spots[[i]]$obs>initial_threshold))
# # threshold selection (try 12)
# par(mfrow = c(2, 2))
# mrlplot(events0[, "obs"]) # should be linear after threshold
# abline(v = initial_threshold, col = "green")
# diplot(events0) # should be close to 1
# abline(v = initial_threshold, col = "green")
# tcplot(events0[, "obs"], which = 1) # should be constant after threshold
# abline(v = initial_threshold, col = "green")
# tcplot(events0[, "obs"], which = 2) # should be constant after threshold
# abline(v = initial_threshold, col = "green")
# graphics.off()


# thresholds[i] <- initial_threshold
# lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>initial_threshold),2])/length(spots[[i]]$obs)

N<-length(data[[1]])
datamatrix<-matrix(NA,N,m)
for(i in 1:N){
  for(j in 1:m){
    datamatrix[i,j]<-spots[[j]][i,2]
  }
}


FOMC_logit_GPD_joint_density <- nimbleFunction(
  run = function(x1=double(0),x2=double(0),threshold=double(0),sigma=double(0),xi=double(0),alpha=double(0), lambda=double(0)){
    returnType(double(0))
    z1 <- 1 + xi*(x1 - threshold)/sigma
    z2 <- 1 + xi*(x2- threshold)/sigma
    inv_xi <- 1/xi
    inv_alpha <- 1/alpha 
    inv_lambda <- 1/lambda
    if(x1>threshold & x2>threshold){
     return(-(inv_lambda * z1^inv_xi)^(-inv_alpha)*(inv_lambda * z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
            + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
            + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2) + (inv_lambda*z1^inv_xi)^(-inv_alpha)*(inv_lambda*z2^inv_xi)^(-inv_alpha)*((inv_lambda*z2^inv_xi)^(-inv_alpha) 
            + (inv_lambda*z1^inv_xi)^(-inv_alpha))^alpha/(alpha*sigma^2*z1*z2*((inv_lambda*z2^inv_xi)^(-inv_alpha) + (inv_lambda*z1^inv_xi)^(-inv_alpha))^2))
    
  }
  
    if(x1>threshold & x2<=threshold){
     return((1/(lambda*sigma))*(z1^(inv_xi-1))*((inv_lambda*z1^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z1^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)))
    
  }
    if(x1<=threshold & x2>threshold){
     return((1/(lambda*sigma))*(z2^(inv_xi-1))*((inv_lambda*z2^inv_xi)^(-inv_alpha -1))*(((inv_lambda*z2^inv_xi)^(-inv_alpha) + inv_lambda^(-inv_alpha))^(alpha-1)))
    
  }
    else{    
     return(1-((inv_lambda^(-1/alpha) + inv_lambda^(-1/alpha))^(alpha)))
    
    
  }
  
})

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
    y[1:N,i] ~ dmyModel_lpost(threshold[i],exp(logsigma[i]),xi[i],lambda[i],alpha[i])}
  for(i in 1:M){
    logsigma[i] ~ dnorm(a_sigma,phi_sigma)
    xi[i] ~ dnorm(a_xi,phi_xi)
    constraint_data[i] ~ dconstraint( xi[i] + exp(logsigma[i])/(max(y[1:N,i])-threshold[i])>0 & threshold[i]>50 & threshold[i]<100)
    alpha[i] ~ dunif(0,1)
    lambda[i] ~ dunif(0.001,0.05)
    threshold[i]~dnorm(80,0.1)}
  a_sigma ~ dnorm(b,c)
  a_xi ~ dnorm(b,c)
#  phi_sigma ~ dgamma(f,f)
#  phi_xi ~ dgamma(f,f)
  phi_sigma ~ dunif(d,e)
  phi_xi ~ dunif(f,g)
})


pumpConsts <- list(M= m, N=N, b=0, c=10^-3,d=5*10^-3,e=15*10^-3,f=10^-2,g=15*10^-2)
pumpData <-list(y = datamatrix)
pumpInits<-list(phi_sigma=1, phi_xi=1)  #I valori iniziali servono solo se vengono usate delle gamma

Hmodel <- nimbleModel(code=code, name ="Hmodel", constants = pumpConsts, data=pumpData)


mcmc.out <- nimbleMCMC(model = Hmodel,
                       niter = 6000, nchains = 2, thin = 10, nburnin = 1000, 
                       monitors = c("a_sigma","a_xi","phi_sigma","phi_xi", "logsigma", "xi", "alpha", "threshold", "lambda"),
                       summary = TRUE, WAIC = TRUE, samplesAsCodaMCMC  = TRUE, setSeed = TRUE)
help(nimbleMCMC)


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

sigmacb<-exp(samples1[,4])
xicb <- samples1[,7]

estim = cbind( sigmacb, xicb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]))
ests2 <- array(0, c(2, 2))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), quantile(estim[, 2], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), quantile(estim[, 2], 0.975))
results <- list(posterior = estim, data = y, postmean = ests, 
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

sampl <- qgpd(1 - 1/t, x$posterior[, 2], thresholds[1], x$posterior[, 1])
res <- quantile(sampl, 0.5)
ta <- seq(1, k, 1)
n <- length(ta)
li <- array(0, c(n))
ls <- array(0, c(n))
pred <- array(0, c(n))
for (s in 1:n) {
  sampl <- qgpd(1 - 1/s, x$posterior[, 2], thresholds[1], x$posterior[, 1])
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

