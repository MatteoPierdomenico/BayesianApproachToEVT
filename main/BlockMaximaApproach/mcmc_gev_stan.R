library(lubridate)
library(evir)
library(rstan)
library(coda)
library(Rcpp)
library(ismev)
library(loo)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)

#Carico il dataset

Ireland <- read.csv('C:/Users/matte/OneDrive/Documenti/GitHub/Bayesian-Approach-to-EVT/Dataset/Ireland_daily_from1990.csv')
head(Ireland)
unique(Ireland$spot)

# fix a location
county <- Ireland %>%
  filter(spot == "clare") %>%
  select(time = year, obs = `hg`)

summary(county$obs) # note 3rd Qu. vs Max.
summary(county$time)

ggplot(county, aes(x=time, y=obs)) +
  geom_point()


#MCMC GEV with STAN: vento estremo

#BlockMaximaApproach
length_block <- 30     #Approximately, 30 is a month
N <- round(length(county$time)/length_block)
Extremes <- rep(NaN, N+1)
for(i in 1:N){
  Extremes[i] <- max(county$obs[(length_block*(i-1)):(length_block*(i))])
}
Extremes[N+1]<-max(county$obs[(30*N):length(county$obs)])


y <- Extremes
N <- length(y)

data_win <- list(N = N,  y = y)     #genero i valori da dare in input a stan
fit <- stan(file = "GEV.stan", data = data_win, warmup = 1500, iter = 6000, 
            chains = 4, cores = 2, thin = 10,seed = 215)              #nel file GEV.stan, uso la log posterior        
is(fit)
print(fit, par = c('mu', 'sigma', 'xi')) 




#PLOT dell 3 catene generate per ogni parametro
rstan::traceplot(fit, pars='mu', inc_warmup = TRUE)
rstan::traceplot(fit, pars='sigma', inc_warmup = TRUE)
rstan::traceplot(fit, pars='xi', inc_warmup = TRUE)


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

#PLOT della stima dei parametri

plot_post <- fit %>% 
  rstan::extract(c("mu", "sigma", "xi")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  theme(legend.position="bottom")

pairs(fit, pars = c('mu', 'sigma', 'xi'), condition = 0.5)

#CODA analysis of the chains

coda_chain <- As.mcmc.list(fit, pars = c("mu","sigma", "xi"))
summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
plot(coda_chain)
acfplot(coda_chain)
cumuplot(coda_chain)

geweke.diag(coda_chain)
geweke.plot(coda_chain, frac1 = 0.1, frac2 = 0.5, nbins = 20)


###
#CREATING A LIST OF THE MAIN RESULTS OF THE CHAIN
###


estim = cbind(mucb, sigmacb, xicb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]), mean(estim[, 
                                                         3]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]), median(estim[, 
                                                                3]))
ests2 <- array(0, c(2, 3))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), quantile(estim[, 
                                                            2], 0.025), quantile(estim[, 3], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), quantile(estim[, 
                                                            2], 0.975), quantile(estim[, 3], 0.975))
results <- list(posterior = estim, data = y, postmean = ests, 
            postmedian = ests1, postCI = ests2, block = 100)
names(results$postmean) <- c("mu", "sigma", "xi")
names(results$postmedian) <- c("mu", "sigma", "xi")
dimnames(results$postCI) <- list(c("lower bound", "upper bound"), 
                             c("mu", "sigma", "xi"))
class(results) <- "gevp"

###
#RETURN LEVEL COMPUTATION AND PLOT
###


t <- 2  #start ret level
k <- 100   #end ret level
x <- results
sampl <- qgev(1 - 1/t, x$posterior[, 3], x$posterior[, 
                                                    1], x$posterior[, 2])
res <- quantile(sampl, 0.5)
ta <- seq(1, k, 1)
n <- length(ta)
li <- array(0, c(n))
ls <- array(0, c(n))
pred <- array(0, c(n))
for (s in 1:n) {
  sampl <- qgev(1 - 1/s, x$posterior[, 3], x$posterior[, 
                                                      1], x$posterior[, 2])
  li[s] <- quantile(sampl, 0.025)
  ls[s] <- quantile(sampl, 0.975)
  pred[s] <- quantile(sampl, 0.5)
}
ggplot(data = data.frame( level= ta, intensity=pred), mapping = aes(x=level,y=intensity)) +
  geom_line( color = "red") +
  xlab("level") + 
  ylab("wind speed (kmh)") + 
  ggtitle("Return level") +
  xlim(0,100)+
  ylim(80,180) +
  geom_line(aes(y = li), color = "black", linetype = "dotted") +
  geom_line(aes(y = ls), color = "black", linetype = "dotted")


##
#log predictive density
##

dat <- x$data
linf <- max(min(dat) - 1, 0)
lsup <- 11 * max(dat)/10
dat1 <- seq(linf, lsup, (lsup - linf)/100)
n <- length(dat1)
int <- length(x$posterior[, 1])
res <- array(0, c(n))
for (i in 1:n) {
  for (j in 1:int) {
    if ((x$posterior[j, 3] > 0) && (dat1[i] > (x$posterior[j, 1] - x$posterior[j, 2]/x$posterior[j, 3]))) 
      res[i] <- res[i] + (1/int) * dgev(dat1[i], x$posterior[j, 3], x$posterior[j, 1], x$posterior[j, 2])
    if ((x$posterior[j, 3] < 0) && (dat1[i] < (x$posterior[j, 
                                                           1] - x$posterior[j, 2]/x$posterior[j, 3])))
      res[i] <- res[i] + (1/int) * dgev(dat1[i], 
                                       x$posterior[j, 3], x$posterior[j, 1], x$posterior[j, 2])
  }
}
dataf <- data.frame(log_pos = dat)
dataf2 <- data.frame(data = dat1, res = res)

ggplot(dataf, aes(x=log_pos)) + 
  geom_histogram(aes(y = ..density.. ), bins = 20, color = "black", fill = "lightblue") +
  xlab("data") + 
  ylab("density") +
  ggtitle("density") +
  geom_line(aes(x=data,y = res), dataf2,color = "red") 




###
#WAIC criteria to choose the model
###


loglik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)

waic(loglik)
