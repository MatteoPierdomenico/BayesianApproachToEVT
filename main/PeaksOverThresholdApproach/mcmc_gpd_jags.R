library(runjags)
library(coda)
library(tidyverse)
library(POT)
library(evir)
library(parallel)

load.runjagsmodule()
clust <- makeCluster(detectCores())

# load world data
spot_data_0 <- read_csv("../../Dataset/Ireland_daily_from1990.csv") %>%
  filter(spot == "kerry") %>%                    #Here is where to change the county to be analyzed
  select(time = year, obs = `hg`)

ggplot(spot_data_0, aes(x=time, y=obs)) +
  geom_point()

acf(spot_data_0$obs)

# pick a very low threshold and decluster
t1 <- 30
spot_data_1 <- as_tibble(clust(as.data.frame(spot_data_0), 
                             u = t1, 
                             tim.cond = 3/365, 
                             clust.max = T,
                             plot = T))
acf(spot_data_1$obs)

ggplot(spot_data_1, aes(x=time, y=obs)) +
  geom_point()

# jitter
delta <- 1.852         # delta between each distinct value
jitter <- runif(length(spot_data_1$obs), min=0, max=delta)
spot_data <- spot_data_1 %>%
  mutate(obs = obs + jitter)

ggplot(spot_data, aes(x=time, y=obs)) +
  geom_point()

# select a threshold
threshold <- 70
par(mfrow = c(2, 2))
mrlplot(spot_data$obs) # should be linear after threshold
abline(v = threshold, col="green")
tcplot(spot_data$obs, which = 1) # should be constant after threshold
abline(v = threshold, col="green")
tcplot(spot_data$obs, which = 2) # should be constant after threshold
abline(v = threshold, col="green")
graphics.off()

spot_data %>%
  mutate(Color = ifelse(obs>threshold, "black", "gray")) %>%
  ggplot(aes(x=time, y=obs, color=Color)) +
  geom_point() +
  scale_color_identity() +
  geom_hline(yintercept=threshold, color="green")




# frequentist analysis

mle <- fitgpd(spot_data$obs, thresh = threshold, est = "mle")
npy <- length(spot_data$obs)/(diff(range(spot_data$time)))
par(mfrow = c(2, 2))
plot(mle, npy = npy)
graphics.off()
mle
# Estimates:
#   scale    shape
# 14.6751  -0.1008
# 
# Standard Errors:
#   scale    shape  
# 1.91443  0.08554


## 
## BAYESIAN
## 

y <- spot_data$obs
max(y)
min(y)
N <- length(y)
Ones <- rep(1,N)

t2 <- 140

data_win <- list("N" = N, "y" = y, 
                 "ymax" = max(y), "Ones" = Ones,
                 "t1" = t1, "t2" = t2,
                 "t_guess" = threshold, "t_guess_var" = 900)

#inits
u <- function(chain) +
  return( switch(chain, "1"= t1+1, "2"= t2-1, "3"=threshold, "4" =threshold) )

posterior <- run.jags("GPD.jags", 
                      data = data_win,
                      n.chains = 4,
                      adapt = 1500,
                      burnin = 1500,
                      sample = 3000,
                      thin = 10,
                      cl = clust)
summary(posterior)
plot(posterior, vars="sigma")
plot(posterior, var="u")
plot(posterior, var="xi")
 


coda_chain <- as.mcmc.list(posterior, pars = c("sigma", "u", "xi"))
summary(coda_chain)
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
plot(coda_chain)
acfplot(coda_chain)
cumuplot(coda_chain)
#geweke.diag(coda_chain)
#geweke.plot(coda_chain, frac1 = 0.1, frac2 = 0.5, nbins = 20)






###
#CREATING A LIST OF THE MAIN RESULTS OF THE CHAIN
###
data <- runjags::extract(posterior, "model")

sigmacb <- coda_chain[[4]][,1]
ucb <- coda_chain[[4]][,2]
xicb <- coda_chain[[4]][,3]


estim = cbind(sigmacb, xicb, ucb)
ests <- c(mean(estim[, 1]), mean(estim[, 2]), mean(estim[, 3]))
ests1 <- c(median(estim[, 1]), median(estim[, 2]), median(estim[, 2]))
ests2 = array(0, c(2, 3))
ests2[1, ] <- c(quantile(estim[, 1], 0.025), 
                quantile(estim[, 2], 0.025),
                quantile(estim[, 3], 0.025))
ests2[2, ] <- c(quantile(estim[, 1], 0.975), 
                quantile(estim[, 2], 0.975),
                quantile(estim[, 3], 0.975))
results <- list(posterior = estim, data = y, postmean = ests, 
                postmedian = ests1, postCI = ests2, block = 100)
names(results$postmean) <- c("sigma", "xi", "u")
names(results$postmedian) <- c("sigma", "xi", "u")
dimnames(results$postCI) <- list(c("lower bound", "upper bound"), 
                                 c("sigma", "xi", "u"))
class(results) <- "gpdp"

###
#RETURN LEVEL COMPUTATION AND PLOT
###

t <- 1  # start return level
k <- 100 # finish return level
ta <- t:k
n <- length(ta)

li <- numeric(n)
ls <- numeric(n)
pred <- numeric(n)
li.mle <- numeric(n)
ls.mle <- numeric(n)
pred.mle <- numeric(n)
for (s in 1:n) {
  p <- rp2prob(ta[s], npy)[,"prob"]
  sampl <- evir::qgpd(p, 
                      mu=results$posterior[, 3], 
                      beta=results$posterior[, 1], xi=results$posterior[, 2])
  li[s] <- quantile(sampl, 0.025)
  ls[s] <- quantile(sampl, 0.975)
  pred[s] <- quantile(sampl, 0.5)
  
  conf.mle <- gpd.firl(mle, p, conf = 0.95)
  li.mle[s] <- conf.mle[["conf.inf"]]
  ls.mle[s] <- conf.mle[["conf.sup"]]
  pred.mle[s] <- evir::qgpd(p, 
                            mu=threshold, 
                            beta=mle$fitted.values[1], xi=mle$fitted.values[2])
}

data.frame(ta,
           pred, li, ls,
           pred.mle, li.mle, ls.mle) %>%
  ggplot(aes(x=ta,y=pred)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin=li, ymax=ls), fill="red", alpha=0.2) +
  #geom_line(aes(x=ta, y=pred.mle), color="black", alpha=0.5) +
  #geom_ribbon(aes(ymin=li.mle, ymax=ls.mle), fill="black", alpha=0.1) +
  #scale_x_continuous(trans="log", limits=c(1,k)) +
  ylim(80,210) +
  xlab("return period (years)") + 
  ylab("return level (wind speed (kmh))") + 
  ggtitle("Return level")



##
#log predictive density
##

x <- results
dat = x$data
linf = max(min(dat) - 1, 0)
lsup = 11 * max(dat)/10
dat1 = seq(linf, lsup, (lsup - linf)/100)
n = length(dat1)
int = length(x$posterior[, 1])
res = array(0, c(n))
for (i in 1:n) {
  for (j in 1:int) {
    if (dat1[i] < x$posterior[j, 3]) {
      res[i] = res[i] + (1/int) * 0.5 * dunif(dat1[i], min(dat1), x$posterior[j, 3])
    } else {   
      res[i] = res[i] + (1/int) * 0.5 * POT::dgpd(dat1[i], 
                                            shape=x$posterior[j, 2], 
                                            loc=x$posterior[j, 3], 
                                            scale=x$posterior[j, 1])
    }
    #if ((x$posterior[j, 3] > 0) && (dat1[i] > (x$posterior[j, 1] - x$posterior[j, 2]/x$posterior[j, 3]))) 
    #  res[i] = res[i] + (1/int) * dgpd(dat1[i], x$posterior[j, 3], x$posterior[j, 1], x$posterior[j, 2])
    #if ((x$posterior[j, 3] < 0) && (dat1[i] < (x$posterior[j, 
    #                                                       1] - x$posterior[j, 2]/x$posterior[j, 3])))
    #  res[i] = res[i] + (1/int) * dgpd(dat1[i], 
    #                                   x$posterior[j, 3], x$posterior[j, 1], x$posterior[j, 2])
  }
}
#hist(dat, freq = F, ylim = c(min(res), max(res)), main = NULL, 
#     xlab = "data", ylab = "density")
#lines(dat1, res)


dataf <- data.frame(pos = dat)
dataf2 <- data.frame(data = dat1, res = res)

ggplot(dataf, aes(x=pos)) + 
  stat_bin(aes(y = ..density.. ), bins = 20, 
           color = "black", fill = "lightblue",
           boundary=1) +
  #geom_density() +
  xlab("data") + 
  ylab("density") +
  ggtitle("density") +
  #xlim(c(30, max(dat))) +
  geom_line(aes(x=data,y = res), dataf2,color = "red") 






