library(ggplot2)
library(POT)
library(lubridate)
library(evir)
library(coda)
library(nimble)
library(tidyr)


Ireland <- read.csv('C:/Users/matte/OneDrive/Documenti/GitHub/Bayesian-Approach-to-EVT/Dataset/Ireland_daily_from1990.csv')
head(Ireland)
unique(Ireland$spot)
Ireland<-separate(Ireland,"date", c("day","month","y"),sep = "-")
unique(Ireland$month)

S <- length(unique(Ireland$spot))
M <- length(unique(Ireland$month))

Ireland_monthly <- list()
counties <- as.vector(unique(Ireland$spot))
months <- as.vector(unique(Ireland$month))
r <- rep(NA,M)

for(i in 1:M){ 
  monthly <- Ireland[which(Ireland$month==months[i]),]
  Ireland_monthly[[i]] <- data.frame(time=monthly$year, obs=monthly$hg, county = monthly$spot)
  r[i] <- length(Ireland_monthly[[i]]$obs)
}



for(i in 1:M){
  data<-matrix(NA,(r[i]/5),S)
  for(j in 1:S){
    county<-Ireland_monthly[[i]][which(Ireland_monthly[[i]]$county==counties[j]),]
    for(n in 1:(r[i]/5)){
      data[n,j]<-county$obs[n]
    }
  }
  Ireland_monthly[[i]]<-data
}


thresholds <- matrix(72,M,S)
lambda <- matrix(NA,M,S)
thresholds[8,5]<-60
for(i in 1:M){
  for(j in 1:S){
    lambda[i,j] <- length(Ireland_monthly[[i]][which(Ireland_monthly[[i]][,j]>thresholds[i,j]),j])/r[i]
  }
}

i=5 #loop manuale (modificare i fino a 5)

# initial data selection
initial_threshold <- 85
events0 <- c(spots[[i]][which(spots[[i]]$obs>initial_threshold),])
events0 <- data.frame(time = events0$time,obs= events0$obs, idx=which(spots[[i]]$obs>initial_threshold))
# threshold selection (try 12)
par(mfrow = c(2, 2))
mrlplot(events0[, "obs"]) # should be linear after threshold
abline(v = initial_threshold, col = "green")
diplot(events0) # should be close to 1
abline(v = initial_threshold, col = "green")
tcplot(events0[, "obs"], which = 1) # should be constant after threshold
abline(v = initial_threshold, col = "green")
tcplot(events0[, "obs"], which = 2) # should be constant after threshold
abline(v = initial_threshold, col = "green")
graphics.off()


thresholds[i] <- initial_threshold
lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>initial_threshold),2])/length(spots[[i]]$obs)


maxim <- matrix(NA,M,S)

for(i in 1:M){
  for(j in 1:S){
    maxim[i,j] <- max(Ireland_monthly[[i]][,j])
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

GPD_density <- nimbleFunction(
  run=function(x=double(0),threshold=double(0),scale=double(0), xi=double(0), lambda=double(0)){
    returnType(double(0))
    if(x>= threshold){
      if (abs(xi) > 1e-15)
      { return((lambda/scale) * (1 + (xi * (x - threshold))/scale)^((-1/xi) - 1))}
      else
       {return( lambda/scale*exp(-(x-threshold)/scale))}  
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
  for(j in 1:S){
    y_dec[1:len[1],j] ~ dmyModel_lpost(threshold[1,j],exp(logsigma[1,j]),xi[1,j],lambda[1,j],alpha[j])
    y_jan[1:len[2],j] ~ dmyModel_lpost(threshold[2,j],exp(logsigma[2,j]),xi[2,j],lambda[2,j],alpha[j])
    y_feb[1:len[3],j] ~ dmyModel_lpost(threshold[3,j],exp(logsigma[3,j]),xi[3,j],lambda[3,j],alpha[j])
    y_mar[1:len[4],j] ~ dmyModel_lpost(threshold[4,j],exp(logsigma[4,j]),xi[4,j],lambda[4,j],alpha[j])
    y_apr[1:len[5],j] ~ dmyModel_lpost(threshold[5,j],exp(logsigma[5,j]),xi[5,j],lambda[5,j],alpha[j])
    y_may[1:len[6],j] ~ dmyModel_lpost(threshold[6,j],exp(logsigma[6,j]),xi[6,j],lambda[6,j],alpha[j])
    y_jun[1:len[7],j] ~ dmyModel_lpost(threshold[7,j],exp(logsigma[7,j]),xi[7,j],lambda[7,j],alpha[j])
    y_jul[1:len[8],j] ~ dmyModel_lpost(threshold[8,j],exp(logsigma[8,j]),xi[8,j],lambda[8,j],alpha[j])
    y_aug[1:len[9],j] ~ dmyModel_lpost(threshold[9,j],exp(logsigma[9,j]),xi[9,j],lambda[9,j],alpha[j])
    y_sep[1:len[10],j] ~ dmyModel_lpost(threshold[10,j],exp(logsigma[10,j]),xi[10,j],lambda[10,j],alpha[j])
    y_oct[1:len[11],j] ~ dmyModel_lpost(threshold[11,j],exp(logsigma[11,j]),xi[11,j],lambda[11,j],alpha[j])
    y_nov[1:len[12],j] ~ dmyModel_lpost(threshold[12,j],exp(logsigma[12,j]),xi[12,j],lambda[12,j],alpha[j])}
  
  for(j in 1:S){
    for(i in 1:M){
      constraint_data[i,j] ~ dconstraint( xi[i,j] > -sigma[i,j]/(maxim[i,j]-threshold[i,j]))
      xi[i,j] <- gamma_xi[i]+ eps_xi[j]
      sigma[i,j] <- exp(logsigma[i,j])
      logsigma[i,j] <- gamma_sigma[i]+ eps_sigma[j]
    }
  }
  for(i in 1:M){
    gamma_sigma[i] ~ dnorm(0,tau_sigma)
    gamma_xi[i] ~ dnorm(0,tau_xi)
  }
  for(j in 1:S){
    eps_sigma[j] ~ dnorm(a_sigma,phi_sigma)
    eps_xi[j] ~ dnorm(a_xi,phi_xi)
    alpha[j] ~ dunif(0,1)
  }
  a_sigma ~ dnorm(b,c)
  a_xi ~ dnorm(b,c)
 # phi_sigma ~ dgamma(f,f)
 # phi_xi ~ dgamma(f,f)
  phi_sigma ~ dunif(d,e)
  phi_xi ~ dunif(f,g)

  tau_sigma ~ dunif(d,e)
  tau_xi ~ dunif(f,g)
  
})



pumpConsts <- list(M = M, S = S, threshold = thresholds, lambda = lambda, maxim = maxim, len = r/5, b=0, c=10^-3,d=5*10^-3,e=20*10^-3,f=1*10^-2,g=15*10^-2)
pumpData <-list(y_dec = Ireland_monthly[[1]], y_jan = Ireland_monthly[[2]], y_feb = Ireland_monthly[[3]],
                y_mar = Ireland_monthly[[4]], y_apr = Ireland_monthly[[5]], y_may = Ireland_monthly[[6]],
                y_jun = Ireland_monthly[[7]], y_jul = Ireland_monthly[[8]], y_aug = Ireland_monthly[[9]],
                y_sep = Ireland_monthly[[10]], y_oct = Ireland_monthly[[11]], y_nov = Ireland_monthly[[12]])

Hmodel <- nimbleModel(code=code, name ="Hmodel", constants = pumpConsts, data=pumpData)

configureMCMC(Hmodel)

mcmc.out <- nimbleMCMC(model = Hmodel,
                       niter = 10000, nchains = 4, thin = 10, nburnin = 3000, 
                       monitors = c("a_sigma","a_xi","phi_sigma","phi_xi","tau_sigma",
                                    "tau_xi","gamma_sigma", "gamma_xi", "eps_sigma",
                                    "eps_xi","logsigma","sigma", "xi", "alpha"),
                       summary = TRUE, WAIC = TRUE, setSeed = TRUE, samplesAsCodaMCMC = TRUE)



coda_chain <- mcmc.out$samples
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

sigmacb<-mcmc.out$samples$chain4[,"sigma[1, 4]"]
xicb <- mcmc.out$samples$chain4[,"xi[1, 4]"]

estim = cbind( sigmacb, xicb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]))
ests2 <- array(0, c(2, 2))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), quantile(estim[, 2], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), quantile(estim[, 2], 0.975))
results <- list(posterior = estim, data = Ireland_monthly[[1]][,4], postmean = ests, 
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

