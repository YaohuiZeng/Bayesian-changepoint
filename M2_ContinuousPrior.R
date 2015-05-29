## Bayesian Course Project: Bayesian modeling for a changepoint problem
## Yaohui Zeng <yaohui-zeng@uiowa.edu>, May 13, 2015

## Taken from the OpenBUGS example:
##      Stagnant: a changepoint problem and an illustration of how NOT
##      to do MCMC!

## Model 2: continuous prior for the changepoint value, x_k.

## Slice sampling (written by Dr. Brian Smith <brian-j-smith@uiowa.edu>)
slice <- function(x0, f, width, log=FALSE, ...) {
  n <- length(x0)
  y <- ifelse(log, f(x0, ...) - rgamma(1, 1, 1), runif(1, 0, f(x0, ...)))
  l <- x0 - runif(n, 0, width)
  r <- l + width
  repeat {
    x1 <- runif(n, l, r)
    if(f(x1, ...) > y) return(x1)
    else {
      idx <- (x1 < x0)
      l[idx] <- x1[idx]
      r[!idx] <- x1[!idx]
    }
  }
}

## Slice Sampling Density Kernel for xk
gxk <- function(xk, a, b, pars) {
  if(xk < a || xk > b) return(-Inf)
  s2 <- pars$s2
  beta1 <- pars$beta1
  beta2 <- pars$beta2
  alpha <- pars$alpha
  k <- sum(x <= xk)
  
  delta <- sum((Y[1:k] - alpha - beta1 * x[1:k] + beta1 * xk)^2) + 
    sum((Y[(k+1):N] - alpha - beta2 * x[(k+1):N] + beta2 * xk)^2)
  
  -1 / (2*s2) * delta
}

## Rejection sampling
freject<- function(xk, pars){
  s2 <- pars$s2
  beta1 <- pars$beta1
  beta2 <- pars$beta2
  alpha <- pars$alpha
  
  k <- sum(x <= xk)
  delta <- sum((Y[1:k] - alpha - beta1 * x[1:k] + beta1 * xk)^2) + 
    sum((Y[(k+1):N] - alpha - beta2 * x[(k+1):N] + beta2 * xk)^2)
  
  -1 / (2*s2) * delta
}

## Rejection sampling support
support <- function(xk, pars){
  (xk > -1.3) & (xk < 1.1)
}


gibbs.sampling.m2 <- function(para.init, seed, L = 1000,
                              method = 1, verbose = FALSE) {
  set.seed(seed)
  
  ## Hyperparameters
  mu.alpha <- 0; s2.alpha <- 1e6
  mu.beta <- 0; s2.beta <- 1e6
  a.s2 <- 0.001; b.s2 <- 0.001
  a <- -1.3; b <- 1.1 # parameters of Unif(a, b) for changepoint xk
  
  ## Model Parameters and Initial Values
  alpha <- para.init$alpha
  beta1 <- para.init$beta1
  beta2 <- para.init$beta2
  s2    <- para.init$s2
  xk    <- para.init$xk
  
  ## Divide data into two groups based on xk
  k <- sum(x <= xk)
  
  ## save sampled values
  parms <- matrix(NA, ncol = 6, nrow = L)
  parms[1, ] <- c(alpha, beta1, beta2, s2, xk, k)
  
  for (i in 2:L) {
    # alpha
    b.a <- 1 / (N / s2 + 1 / s2.alpha)
    if (k < N) {
      delta <- sum(Y) - beta1 * sum(x[1:k] - xk) - 
        beta2 * sum(x[(k+1):N] - xk)
    } else {
      delta <- sum(Y) - beta1 * sum(x[1:k] - xk)
    }
    a.a <- b.a * (delta / s2 + mu.alpha / s2.alpha)
    alpha <- rnorm(1, mean = a.a, sd = b.a ^ 0.5)
    
    # beta1
    b.b1 <- 1 / (sum((x[1:k] - xk) ^ 2) / s2 + 1 / s2.beta)
    a.b1 <- b.b1 * (sum((Y[1:k] - alpha) * (x[1:k] - xk)) / s2 + 
                      mu.beta / s2.beta)
    beta1 <- rnorm(1, mean = a.b1, sd = b.b1 ^ 0.5)
    
    # beta2
    if (k < N) {
      b.b2 <- 1 / (sum((x[(k+1):N] - xk) ^ 2) / s2 + 1 / s2.beta)
      a.b2 <- b.b2 * (sum((Y[(k+1):N] - alpha) * (x[(k+1):N] - xk)) / s2 
                      + mu.beta / s2.beta)
    } else {
      b.b2 <- 1 / (1 / s2.beta) # = s2.beta
      a.b2 <- b.b2 * (mu.beta / s2.beta) # = mu.beta
    }
    beta2 <- rnorm(1, mean = a.b2, sd = b.b2 ^ 0.5)
    
    # s2
    a.3 <- a.s2 + N / 2
    if (k < N) {
      b.3 <- sum((Y[1:k] - alpha - beta1 * (x[1:k] - xk)) ^ 2) / 2 + 
        sum((Y[(k+1):N] - alpha - beta2 * (x[(k+1):N] - xk)) ^ 2) / 2 + b.s2
    } else {
      b.3 <- sum((Y[1:k] - alpha - beta1 * (x[1:k] - xk)) ^ 2) / 2 + b.s2
    }
    s2 <- 1 / rgamma(1, a.3, b.3)
    
    if (method == 1) { ## slice sampling
      pars <- list(
        s2 = s2,
        beta1 = beta1,
        beta2 = beta2,
        alpha = alpha
      )
      xk <- slice(xk, f = gxk, width = 10, log = TRUE, a = a, b = b, pars = pars)
    } else if (method == 2) { ## rejection sampling
      parss <- list(
        alpha = alpha,
        beta1 = beta1,
        beta2 = beta2,
        s2    = s2
      )
      xk <- arms(0, freject, support, pars=parss, 1)    
    } else {
      stop('Method must be equal 1 or 2. 1 = slice sampling; 2 = rejection sampling!')
    }
    
    # Update index k
    k <- sum(x <= xk)
    
    if(verbose) {
      print(i)
      print(c(alpha, beta1, beta2, s2, xk, k))
    }
    
    ## store sampled values
    parms[i, ] <- c(alpha, beta1, beta2, s2, xk, k)
    
  }
  parms <- data.frame(parms)
  names(parms) <- c("alpha", "beta1", "beta2", "s2", "xk", "k")
  parms
}

library(coda)
library(xtable)
library(ars)
library(HI)
data <- list(Y = c(1.12, 1.12, 0.99, 1.03, 0.92, 0.90, 0.81, 0.83, 0.65, 0.67, 
                   0.60, 0.59, 0.51, 0.44, 0.43, 0.43, 0.33, 0.30, 0.25, 0.24, 
                   0.13, -0.01, -0.13, -0.14, -0.30, -0.33, -0.46, -0.43, -0.65),
             x = c(-1.39, -1.39, -1.08, -1.08, -0.94, -0.80, -0.63, -0.63, 
                   -0.25, -0.25, -0.12, -0.12, 0.01, 0.11, 0.11, 0.11, 0.25, 
                   0.25, 0.34, 0.34, 0.44, 0.59, 0.70, 0.70, 0.85, 0.85, 0.99, 
                   0.99, 1.19),
             N = 29)
## Data
Y <<- data$Y
x <<- data$x
N <<- data$N

# plot(data$x, data$Y, xlab = 'x', ylab = 'y', main = 'Stagnant band height data')

# chain 1
para.ch1 <- list(
  alpha = 0.47,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  xk    = 0.5
)

# chain 2
para.ch2 <- list(
  alpha = 0.47,
  beta1 = -0.45,
  beta2 = -1.0,
  s2    = 1 / 5,
  xk    = -0.5
)

## initial comparison
ch1 <- gibbs.sampling.m2(para.init = para.ch1, seed = 12345, L = 10000, 
                         method = 1, verbose = T)
ch1.mcmc <- as.mcmc(ch1)
plot(ch1.mcmc)
summary(ch1.mcmc)

ch2 <- gibbs.sampling.m2(para.init = para.ch2, seed = 12345, L = 10000,
                         method = 1, verbose = T)
ch2.mcmc <- as.mcmc(ch2)
plot(ch2.mcmc)
summary(ch2.mcmc)
cor(ch2[-6])

stagnant.mcmc <- mcmc.list(ch1.burn.mcmc, ch2.burn.mcmc)
plot(stagnant.mcmc, ask=TRUE)
summary(stagnant.mcmc)


ch1 <- gibbs.sampling.m2(para.init = para.ch1, seed = 12345, L = 10000, 
                         method = 3, verbose = F)
ch2 <- gibbs.sampling.m2(para.init = para.ch2, seed = 12345, L = 10000,
                         method = 3, verbose = F)

# compare to OpenBUGS
png('compare2BUGS_alpha_m2_fixed_rejection sampling.png', width = 800, height = 300)
plot(1:10000, ch1$alpha, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'alpha', lty = 3,
     ylim = c(0, 1),
     main = 'Comparison of alpha')
lines(1:10000, ch2$alpha, col = 'blue')
legend('bottomright', legend = c('Chain 1', 'Chain 2'), lty = c(3, 1),
       col = c('red', 'blue'))
dev.off()

png('compare2BUGS_xk_m2_fixed_rejection sampling.png', width = 800, height = 300)
plot(1:10000, ch1$xk, type = 'l', col = 'red', 
     xlab = 'iteration', ylab = 'xk', lty = 3,
     ylim = c(-0.5, 1),
     main = 'Comparison of xk')
lines(1:10000, ch2$xk, col = 'blue')
legend('topright', legend = c('Chain 1', 'Chain 2'), lty = c(3, 1),
       col = c('red', 'blue'))
dev.off()

## more iterations, doesn't help
ch1 <- gibbs.sampling.m2(para.init = para.ch1, seed = 12345, L = 100000, 
                         method = 1, verbose = T)
ch1.mcmc <- as.mcmc(ch1)
plot(ch1.mcmc)
summary(ch1.mcmc)

ch1.burn <- ch1[-(1:10000), ]
ch1.burn.mcmc <- as.mcmc(ch1.burn)
plot(ch1.burn.mcmc)
summary(ch1.burn.mcmc)
cor(ch1[-6])

ch2 <- gibbs.sampling.m2(para.init = para.ch2, seed = 12345, L = 1000000,
                         method = 1, verbose = T)
ch2.mcmc <- as.mcmc(ch2)
plot(ch2.mcmc)
summary(ch2.mcmc)
cor(ch2[-6])

ch2.burn <- ch2[-(1:10000), ]
ch2.burn.mcmc <- as.mcmc(ch2.burn)
plot(ch2.burn.mcmc)
summary(ch2.burn.mcmc)

stagnant.mcmc <- mcmc.list(ch1.burn.mcmc, ch2.burn.mcmc)
plot(stagnant.mcmc, ask=TRUE)
summary(stagnant.mcmc)
