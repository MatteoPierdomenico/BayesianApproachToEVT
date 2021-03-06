## STAN version of the Hierarchical random effects model for daily maximum wind gusts
## varying only by sites.

library(rstan)
library(coda)
library(tidyverse)
library(POT)
library(evir)
library(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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
#It is possible to see an higher lag 1 correlation
plot(pacf(spots[[i]]$obs))      
readline(prompt = "press [enter] to continue or ESC to stop")

graphics.off()
}



thresholds <- rep(NA,m)
lambda <- rep(NA,m)

thresholds <- c(70,70,68,68,63)  #sono state scelte empiricamente

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


for(i in 1:m){
  lambda[i] <- length(spots[[i]][which(spots[[i]]$obs>thresholds[i]),2])/length(spots[[i]]$obs)
}



#MCMC

y <- cbind(spots[[1]][,2],spots[[2]][,2],spots[[3]][,2],spots[[4]][,2],spots[[5]][,2])
N <- length(y[,1])
M <- length(spots)
b <- 0
c <- 10^-3
data_win <- list( N = N, M=M, y = y, threshold = thresholds, lambda=lambda,b=b,c=c)

fit <- stan(file = "Hierarchical_model.stan", data = data_win, warmup = 1000, iter = 1500, 
            chains = 2, thin = 1,seed = 199, init_r = .01) 
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

## Change the first four lines to analyze different sites,
## in the following, for example, we are chosing i=4 which is Kerry

sigmacb <- sigmacb[,4]
xicb <- xicb[,4]
threshold <- thresholds[4]
obs <- y[y[,4]>threshold, 4]

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

