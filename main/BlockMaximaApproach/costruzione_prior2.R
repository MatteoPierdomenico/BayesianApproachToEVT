library(evd)
library(ismev)
library(ppcc)
library(evir)
library(fitdistrplus)
library(rstan)
library(lubridate)
library(coda)
library(Rcpp)
library(loo)
library(betafunctions)
library(tidyverse)

# LOAD DATA

ireland <- read.csv("../../Dataset/Ireland_daily_from1990.csv")
unique(ireland$spot)
# "clare"     "cork"      "dublin"    "kerry"     "westmeath"

cities <- (as.vector(unique(ireland$spot))) 

# SELECT ONE COUNTY (CLARE in this case)

city_data <- ireland %>%
  filter(spot == cities[1]) %>%
  select(obs=hg)

M <- dim(city_data)[1]

#TAKE ONLY THE FIRST HALF TO EVALUATE THE PRIORS

testerP <- city_data[1:(M/2),]

#FIT THE GEV DISTRIBUTION (Block Maxima approach)

Values= evir::gev(testerP,14)      # 14 represents the dim. of the block (two weeks)
y = Values$data
N=length(y)
plot(y)
data_win <- list(N = N,  y = y, mu_mean=0, mu_var=100, sigma_mean=0, 
                 sigma_var=100, xi_mean=0, xi_var=100)                # I generate the values to give as input in stan
fit <- stan(file = "GEV2.stan", data = data_win, warmup = 1500, iter = 6000, 
            chains = 4, cores = 2, thin = 10, seed = 199)              #I use the log-posterior in the file GEV.stan        
is(fit)
print(fit, par = c('mu', 'sigma', 'xi')) 

#PLOT OF THE CHAINS FOR EACH PARAMETER

rstan::traceplot(fit, pars='mu', inc_warmup = TRUE)    #OK
rstan::traceplot(fit, pars='sigma', inc_warmup = TRUE) #OK
rstan::traceplot(fit, pars='xi', inc_warmup = TRUE)    #OK

# EXTRACT THE PARAMETERS

par <- rstan::extract(fit, pars = c("mu", "sigma", "xi", "lp__"), inc_warmup = FALSE)
is(par)
names(par)

mucb <- par$mu
sigmacb <- par$sigma
xicb <- par$xi
lp <- par$lp__
plot(mucb)
plot(sigmacb)
plot(xicb)
plot(lp)

#PLOT OF THE ESTIMATION OF THE PARAMETERS

plot_post <- fit %>% 
  rstan::extract(c("mu", "sigma", "xi")) %>% 
  as_tibble() %>%
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  theme_minimal() + 
  theme(legend.position="bottom")

pairs(fit, pars = c('mu', 'sigma', 'xi'), condition = 0.5)

#CODA ANALYSIS FOR THE CHAIN

coda_chain <- As.mcmc.list(fit, pars = c("mu","sigma", "xi"))
summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

plot(coda_chain)
acfplot(coda_chain)
cumuplot(coda_chain)

geweke.diag(coda_chain)
geweke.plot(coda_chain, frac1 = 0.1, frac2 = 0.5, nbins = 20)


#CREATING A LIST OF THE MAIN RESULTS OF THE CHAIN
###

estim = cbind(mucb, sigmacb, xicb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]), mean(estim[,3]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]), median(estim[,3]))
ests2 = array(0, c(2, 3))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), quantile(estim[,2], 0.025), quantile(estim[, 3], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), quantile(estim[,2], 0.975), quantile(estim[, 3], 0.975))

out <- list(posterior = estim, data = y, postmean = ests, 
                postmedian = ests1, postCI = ests2, block = 100)

names(out$postmean) <- c("mu", "sigma", "xi")
names(out$postmedian) <- c("mu", "sigma", "xi")
dimnames(out$postCI) <- list(c("lower bound", "upper bound"), 
                             c("mu", "sigma", "xi"))

class(out) <- "gevp"

# I ESTIMATE MEAN AND VARIANCE OF EACH PARAMETER

# MU #
mean_mu <- mean(mucb)  
var_mu <- var(mucb)    

# SIGMA #
mean_sigma <- mean(sigmacb) 
var_sigma <- var(sigmacb)   

# XI #
mean_xi <- mean(xicb)  
var_xi <- var(xicb)    

###
#WAIC criteria to choose the model
###

loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)

waic(loglik)

#####################################################################
# Once we have estimated the priors using the first half of our data,
# we want to compare the WAIC index of two models for the second half
# of the data, using different priors:
#
# 1. Non-informative priors as before;
# 2. Estimated informative priors.

testerP2 <- city_data[((M/2)+1):M,] # TAKE THE SECOND HALF OF THE DATA
Values= evir::gev(testerP2,14)      # 14 represents the dim. of the block (two weeks)
y = Values$data
N=length(y)
plot(y)

# 1. NON-INFORMATIVE PRIORS

data_win <- list(N = N,  y = y, mu_mean=0, mu_var=100, sigma_mean=0, 
                 sigma_var=100, xi_mean=0, xi_var=100)                # I generate the values to give as input in stan
fit1 <- stan(file = "GEV2.stan", data = data_win, warmup = 1500, iter = 6000, 
            chains = 4, cores = 2, thin = 10, seed = 199)              #I use the log-posterior in the file GEV.stan        
is(fit1)
print(fit1, par = c('mu', 'sigma', 'xi')) 

#PLOT OF THE CHAINS FOR EACH PARAMETER

rstan::traceplot(fit1, pars='mu', inc_warmup = TRUE)    #OK
rstan::traceplot(fit1, pars='sigma', inc_warmup = TRUE) #OK
rstan::traceplot(fit1, pars='xi', inc_warmup = TRUE)    #OK

# EXTRACT THE PARAMETERS

par1 <- rstan::extract(fit1, pars = c("mu", "sigma", "xi", "lp__"), inc_warmup = FALSE)
is(par1)
names(par1)

mucb1 <- par1$mu
sigmacb1 <- par1$sigma
xicb1 <- par1$xi
lp1 <- par1$lp__
plot(mucb1)
plot(sigmacb1)
plot(xicb1)
plot(lp1)

#CODA analysis of the chains

coda_chain1 <- As.mcmc.list(fit1, pars = c("mu","sigma", "xi"))
summary(coda_chain1)
gelman.diag(coda_chain1, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

plot(coda_chain1)
acfplot(coda_chain1)
cumuplot(coda_chain1)

geweke.diag(coda_chain1)
geweke.plot(coda_chain1, frac1 = 0.1, frac2 = 0.5, nbins = 20)

###
#WAIC criteria to choose the model
###

loglik1 <- extract_log_lik(fit1, parameter_name = "log_lik", merge_chains = TRUE)

waic(loglik1)


# 2. INFORMATIVE PRIORS

data_win <- list(N = N,  y = y, mu_mean=mean_mu, mu_var=var_mu, sigma_mean=mean_sigma, 
                 sigma_var=var_sigma, xi_mean=mean_xi, xi_var=var_xi)                # I generate the values to give as input in stan
fit2 <- stan(file = "GEV2.stan", data = data_win, warmup = 1500, iter = 6000, 
             chains = 4, cores = 2, thin = 10, seed = 199)              #I use the log-posterior in the file GEV.stan        
is(fit2)
print(fit2, par = c('mu', 'sigma', 'xi')) 

#PLOT OF THE CHAINS FOR EACH PARAMETER

rstan::traceplot(fit2, pars='mu', inc_warmup = TRUE)    #OK
rstan::traceplot(fit2, pars='sigma', inc_warmup = TRUE) #OK
rstan::traceplot(fit2, pars='xi', inc_warmup = TRUE)    #OK

# EXTRACT THE PARAMETERS

par2 <- rstan::extract(fit2, pars = c("mu", "sigma", "xi", "lp__"), inc_warmup = FALSE)
is(par2)
names(par2)

mucb2 <- par2$mu
sigmacb2 <- par2$sigma
xicb2 <- par2$xi
lp2 <- par2$lp__
plot(mucb2)
plot(sigmacb2)
plot(xicb2)
plot(lp2)

#CODA analysis of the chains

coda_chain2 <- As.mcmc.list(fit2, pars = c("mu","sigma", "xi"))
summary(coda_chain2)
gelman.diag(coda_chain2, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

plot(coda_chain2)
acfplot(coda_chain2)
cumuplot(coda_chain2)

geweke.diag(coda_chain2)
geweke.plot(coda_chain2, frac1 = 0.1, frac2 = 0.5, nbins = 20)

###
#WAIC criteria to choose the model
###

loglik2 <- extract_log_lik(fit2, parameter_name = "log_lik", merge_chains = TRUE)

waic(loglik2)


