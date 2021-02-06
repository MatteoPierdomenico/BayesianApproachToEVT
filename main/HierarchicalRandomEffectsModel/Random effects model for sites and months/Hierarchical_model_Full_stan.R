## STAN version of the Hierarchical random effects model for daily maximum wind gusts
## varying by months and sites.


library(rstan)
library(coda)
library(tidyverse)
library(POT)
library(evir)
library(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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

## Building a list of matrices, containing monthly observations of daily maximum 
## wind gusts for each site


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


thresholds <- rbind(c(70,   80,   75,   80,   61),
                    c(78,   70,   80,   80,   60),
                    c(79,   72,   75,   76,   56),
                    c(63,   70,   70,   72,   61),
                    c(70,   63,   65,   68,   60),
                    c(62,   60,   60,   64,   58),
                    c(55,   58,   63,   58,   50),
                    c(58,   56,   55,   55,   40),
                    c(56,   57,   60,   56,   45),
                    c(63,   70,   65,   78,   60),
                    c(68,   70,   68,   69,   57),
                    c(75,   70,   70,   84,   64))

## The thresholds above, for every month and site, have been empirically estimated 
## by the following commented lines, which represent the empirical mean residual life plot (mrlplot) 
## and the Threshold Choice Plot (tcplot)

#i=7 # i varies from 1 to 12 for the months (from December to November)
#j=1 # j varies from 1 to 5 for the sites (respectively Clare, Cork, Dublin, Kerry and Westmeath)

## initial data selection
#initial_threshold <- 56
#events0 <- c(Ireland_monthly[[i]][which(Ireland_monthly[[i]][,j]>initial_threshold),j])
#events0 
#events0 <- data.frame(obs= events0, idx=which(Ireland_monthly[[i]][,j]>initial_threshold))
## threshold selection (try 12)
#par(mfrow = c(2, 2))
#mrlplot(events0[,"obs"]) # should be linear after threshold
#abline(v = initial_threshold, col = "green")
##diplot(events0) # should be close to 1
##abline(v = initial_threshold, col = "green")
#tcplot(events0[, "obs"], which = 1) # should be constant after threshold
#abline(v = initial_threshold, col = "green")
#tcplot(events0[, "obs"], which = 2) # should be constant after threshold
#abline(v = initial_threshold, col = "green")
#graphics.off()
#
#
#thresholds[i,j] <- initial_threshold
#lambda[i,j] <- length(Ireland_monthly[[i]][which(Ireland_monthly[[i]][,j]>thresholds[i,j]),j])/length(Ireland_monthly[[i]][,j])
#


## Computing lambda for every choice of threshold. Lambda is the normalized number of
## peaks over threshold for each month and each site.

lambda <- matrix(NA,M,S)

for(i in 1:M){
  for(j in 1:S){
    lambda[i,j] <- length(Ireland_monthly[[i]][which(Ireland_monthly[[i]][,j]>thresholds[i,j]),j])/r[i]
  }
}

## Computing the maximum daily wind gust for each month and each site


maxim <- matrix(NA,M,S)

for(i in 1:M){
  for(j in 1:S){
    maxim[i,j] <- max(Ireland_monthly[[i]][,j])
  }
}


#MCMC with STAN

data_win <- list( S = S, M = M, r = r/5, 
                  y_dec = Ireland_monthly[[1]], y_jan = Ireland_monthly[[2]], y_feb = Ireland_monthly[[3]],
                  y_mar = Ireland_monthly[[4]], y_apr = Ireland_monthly[[5]], y_may = Ireland_monthly[[6]],
                  y_jun = Ireland_monthly[[7]], y_jul = Ireland_monthly[[8]], y_aug = Ireland_monthly[[9]],
                  y_sep = Ireland_monthly[[10]], y_oct = Ireland_monthly[[11]], y_nov = Ireland_monthly[[12]],
                  threshold = thresholds, lambda = lambda,maxim=maxim)


fit <- stan(file = "Hierarchical_model_Full.stan", data = data_win, warmup = 400, iter = 600, 
            chains = 2, thin = 1,seed = 19, init_r= 0.01) 
help(stan)
is(fit)
print(fit, par = c('a_sigma','phi_sigma','a_xi','phi_xi','sigma', 'xi','alpha')) 




#PLOT dell 3 catene generate per ogni parametro
rstan::traceplot(fit, pars='a_sigma', inc_warmup = TRUE)
rstan::traceplot(fit, pars='phi_sigma', inc_warmup = TRUE)
rstan::traceplot(fit, pars='a_xi', inc_warmup = TRUE)
rstan::traceplot(fit, pars='phi_xi', inc_warmup = TRUE)
rstan::traceplot(fit, pars='sigma', inc_warmup = TRUE)
rstan::traceplot(fit, pars='xi', inc_warmup = TRUE)
rstan::traceplot(fit, pars='alpha', inc_warmup = TRUE)

par <- rstan::extract(fit, pars = c('a_sigma','phi_sigma', 'a_xi','phi_xi',"sigma", "xi","alpha", "lp__"), inc_warmup = FALSE)
is(par)
names(par)

sigmacb <- par$sigma
xicb <- par$xi
alphacb <- par$alpha

lp <- par$lp__

plot(sigmacb)
plot(xicb)
plot(alphacb)
plot(lp)

#PLOT della stima dei parametri

plot_post <- fit %>% 
  rstan::extract(c("sigma", "xi")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  theme(legend.position="bottom")

pairs(fit, pars = c( 'sigma', 'xi', 'alpha'), condition = 0.5)

#CODA analysis of the chains

coda_chain <- As.mcmc.list(fit, pars = c('a_sigma','phi_sigma','a_xi','phi_xi',"sigma", "xi", 'alpha'))
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

## Change the first four lines to analyze different sites in different month,
## in example, we are considering i = 2, which is january and j=4 which is Kerry

sigmacb=sigmacb[,8]
xicb=xicb[,8]
threshold=thresholds[2,4]
obs <- Ireland_monthly[[2]][,4]


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
  ylim(80,160) +
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
