## Bayesian Course Project: Bayesian modeling for a changepoint problem
## Yaohui Zeng <yaohui-zeng@uiowa.edu>, May 13, 2015

## Taken from the OpenBUGS example:
##      Stagnant: a changepoint problem and an illustration of how NOT
##      to do MCMC!

## Model 1: discrete prior for the index of the changepoint, k.

## function: compute the data log likelihood given k and other parameters
log.like <- function(k, alpha, beta1, beta2, s2) {
  log.like1 <- 0
  log.like2 <- 0
  for(i in 1:k) {
    log.like1 <- log.like1 + dnorm(Y[i], 
                                   mean = alpha + beta1 * (x[i] - x[k]),
                                   sd = s2, log = TRUE)
  }
  if (k < N) {
    for(i in (k+1):N) {
      log.like2 <- log.like2 + dnorm(Y[i], 
                                     mean = alpha + beta2 * (x[i] - x[k]),
                                     sd = s2, log = TRUE)
    }
  }
  log.like1 + log.like2
}

## function: normalize probabilities
normalize.prob <- function(loglike, epsilon = 1e-16) {
  prob <- rep(0, length(loglike))
  
  if (isTRUE(all.equal(exp(loglike), rep(0, length(loglike))))) {
    ## the log likelihood values are tooooo small, rescale them!!
    max.loglike <- max(loglike)
    loglike.norm <- loglike - max.loglike
    prec <- log(epsilon) - log(length(loglike))
    idx <- which(loglike.norm >= prec)
    
    like <- exp(loglike.norm[idx])
    prob[idx] <- like / sum(like)
  } else {
    like <- exp(loglike)
    prob <- like / sum(like)
  }  
  prob
}

## function: Gibbs sampling
gibbs.sampling <- function(para.init, seed, L = 1000) {
  set.seed(seed)
  
  ## Hyperparameters
  mu.alpha <- 0; s2.alpha <- 1e6
  mu.beta <- 0; s2.beta <- 1e6
  a <- 0.001; b <- 0.001
  
  ## Model Parameters and Initial Values
  alpha <- para.init$alpha
  beta1 <- para.init$beta1
  beta2 <- para.init$beta2
  s2    <- para.init$s2
  k     <- para.init$k
  
  ## save sampled values
  parms <- matrix(NA, ncol = 6, nrow = L)
  parms[1, ] <- c(alpha, beta1, beta2, s2, k, x[k])
  
  for (i in 2:L) {
    # alpha
    b.a <- 1 / (N / s2 + 1 / s2.alpha)
    if (k < N) {
      delta <- sum(Y) - beta1 * sum(x[1:k] - x[k]) - 
        beta2 * sum(x[(k+1):N] - x[k])
    } else {
      delta <- sum(Y) - beta1 * sum(x[1:k] - x[k])
    }
    a.a <- b.a * (delta / s2 + mu.alpha / s2.alpha)
    alpha <- rnorm(1, mean = a.a, sd = b.a ^ 0.5)
    
    # beta1
    b.b1 <- 1 / (sum((x[1:k] - x[k]) ^ 2) / s2 + 1 / s2.beta)
    a.b1 <- b.b1 * (sum((Y[1:k] - alpha) * (x[1:k] - x[k])) / s2 + 
                      mu.beta / s2.beta)
    beta1 <- rnorm(1, mean = a.b1, sd = b.b1 ^ 0.5)
    
    # beta2
    if (k < N) {
      b.b2 <- 1 / (sum((x[(k+1):N] - x[k]) ^ 2) / s2 + 1 / s2.beta)
      a.b2 <- b.b2 * (sum((Y[(k+1):N] - alpha) * (x[(k+1):N] - x[k])) / s2 
                      + mu.beta / s2.beta)
    } else {
      b.b2 <- 1 / (1 / s2.beta) # = s2.beta
      a.b2 <- b.b2 * (mu.beta / s2.beta) # = mu.beta
    }
    beta2 <- rnorm(1, mean = a.b2, sd = b.b2 ^ 0.5)
    
    # s2
    a.3 <- a + N / 2
    if (k < N) {
      b.3 <- sum((Y[1:k] - alpha - beta1 * (x[1:k] - x[k])) ^ 2) / 2 + 
        sum((Y[(k+1):N] - alpha - beta2 * (x[(k+1):N] - x[k])) ^ 2) / 2 + b
    } else {
      b.3 <- sum((Y[1:k] - alpha - beta1 * (x[1:k] - x[k])) ^ 2) / 2 + b
    }
    s2 <- 1 / rgamma(1, a.3, b.3)
    
    # k
    like <- unlist(lapply(1:N, FUN = log.like, alpha, beta1, beta2, s2))
    prob <- normalize.prob(like)
    k <- sample(1:N, 1, replace = TRUE, prob = prob)
    #         if (k == N) k = N - 1
    ## store sampled values
    parms[i, ] <- c(alpha, beta1, beta2, s2, k, x[k])
    
  }
  
  parms <- data.frame(parms)
  names(parms) <- c("alpha", "beta1", "beta2", "s2", "k", "xk")
  parms
}


library(coda)
library(xtable)
data <- list(Y = c(1.12, 1.12, 0.99, 1.03, 0.92, 0.90, 0.81, 0.83, 0.65, 0.67, 
                   0.60, 0.59, 0.51, 0.44, 0.43, 0.43, 0.33, 0.30, 0.25, 0.24, 
                   0.13, -0.01, -0.13, -0.14, -0.30, -0.33, -0.46, -0.43, -0.65),
             x = c(-1.39, -1.39, -1.08, -1.08, -0.94, -0.80, -0.63, -0.63, -0.25, 
                   -0.25, -0.12, -0.12, 0.01, 0.11, 0.11, 0.11, 0.25, 0.25, 0.34, 
                   0.34, 0.44, 0.59, 0.70, 0.70, 0.85, 0.85, 0.99, 0.99, 1.19),
             N = 29)
## Data
Y <<- data$Y
x <<- data$x
N <<- data$N

# plot(data$x, data$Y, xlab = 'x', ylab = 'y', main = 'Stagnant band height data')
para.ch1 <- list(
  alpha = 0.2,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  k     = 16
)

para.ch2 <- list(
  alpha = 0.6,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  k     = 8
)

## initial comparison
ch1 <- gibbs.sampling(para.init = para.ch1, seed = 123456, L = 10000)
ch1.mcmc <- as.mcmc(ch1)
plot(ch1.mcmc)
summary(ch1.mcmc)

acf(ch1$alpha)
acf(ch1$k)
cor(ch1[-6])

ch2 <- gibbs.sampling(para.init = para.ch2, seed = 123456, L = 10000)
ch2.mcmc <- as.mcmc(ch2)
plot(ch2.mcmc)
summary(ch2.mcmc)
cor(ch2[-6])

stagnant.mcmc <- mcmc.list(ch1.mcmc, ch2.mcmc)
plot(stagnant.mcmc, ask=TRUE)
summary(stagnant.mcmc)


# compare to OpenBUGS
png('compare2BUGS_alpha_123456.png', width = 800, height = 300)
plot(1:10000, ch1$alpha, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'alpha',
     ylim = c(0.4, 1.4),
     main = 'Comparison of alpha')
lines(1:10000, ch2$alpha, col = 'blue')
legend('topleft', legend = c('Chain 1', 'Chain 2'), lty = c(1, 1),
       col = c('red', 'blue'))
dev.off()

png('compare2BUGS_k_123456.png', width = 800, height = 300)
plot(1:10000, ch1$k, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'k',
     ylim = c(5, 25),
     main = 'Comparison of k')
lines(1:10000, ch2$k, col = 'blue')
legend('topleft', legend = c('Chain 1', 'Chain 2'), lty = c(1, 1),
       col = c('red', 'blue'))
dev.off()

## two more chains
para.ch3 <- list(
  alpha = 0.1,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  k     = 20
)

para.ch4 <- list(
  alpha = 0.8,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  k     = 4
)

ch3 <- gibbs.sampling(para.init = para.ch3, seed = 12345, L = 10000)
ch3.mcmc <- as.mcmc(ch3)
plot(ch3.mcmc)
summary(ch3.mcmc)

ch4 <- gibbs.sampling(para.init = para.ch4, seed = 12345, L = 10000)
ch4.mcmc <- as.mcmc(ch4)
plot(ch4.mcmc)
summary(ch4.mcmc)

# compare to OpenBUGS
png('compare2BUGS_alpha_2_more.png', width = 800, height = 300)
plot(1:10000, ch1$alpha, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'alpha',
     ylim = c(0.4, 1.4),
     main = 'Comparison of alpha')
lines(1:10000, ch2$alpha, col = 'blue')
legend('topleft', legend = c('Chain 1', 'Chain 2'), lty = c(1, 1),
       col = c('red', 'blue'))
dev.off()

png('compare2BUGS_k_2_more.png', width = 800, height = 300)
plot(1:10000, ch1$k, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'k',
     ylim = c(5, 25),
     main = 'Comparison of k')
lines(1:10000, ch2$k, col = 'blue')
legend('topleft', legend = c('Chain 1', 'Chain 2'), lty = c(1, 1),
       col = c('red', 'blue'))
dev.off()
