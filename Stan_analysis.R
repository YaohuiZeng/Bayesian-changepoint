## RStan implementation

library(coda)
library(rstan)
library(lattice)

N <- 29
x <- c(-1.39, -1.39, -1.08, -1.08, -0.94, -0.8, -0.63, -0.63, -0.25, 
       -0.25, -0.12, -0.12, 0.01, 0.11, 0.11, 0.11, 0.25, 0.25, 0.34, 
       0.34, 0.44, 0.59, 0.7, 0.7, 0.85, 0.85, 0.99, 0.99, 1.19)
Y <- c(1.12, 1.12, 0.99, 1.03, 0.92, 0.9, 0.81, 0.83, 0.65, 0.67, 
       0.6, 0.59, 0.51, 0.44, 0.43, 0.43, 0.33, 0.3, 0.25, 0.24, 0.13, 
       -0.01, -0.13, -0.14, -0.3, -0.33, -0.46, -0.43, -0.65)

# function to convert stan to coda
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

init.par <- list(
  list(alpha = 0.47, beta = c(-0.45, -1.0), sigma = 1/5, x_change = 0.5),
  list(alpha = 0.47, beta = c(-0.45, -1.0), sigma = 1/5, x_change = -0.5)
)

stan_model <- "
    data { 
        int<lower=0> N; 
        real x[N]; 
        real Y[N]; 
    } 

    parameters { 
        real<lower=0> sigma; 
        real<lower=0> alpha; 
        real beta[2]; 
        real<lower=min(x), upper=max(x)> x_change; 
    } 

    model { 
        real mu[N]; 
        alpha ~ normal(0, 1.0E6); 
        beta ~ normal(0, 1.0E6); 
        sigma ~ inv_gamma(0.001, 0.001); 
        x_change ~ uniform(-1.3, 1.1);
  
    for (n in 1:N) 
        mu[n] <- alpha + 
            if_else(x[n] < x_change, beta[1], beta[2]) * (x[n] - x_change); 
    Y ~ normal(mu,sigma); 
    } 
"
    
# 1: no init values; 10 warmup
fit1 <- stan(model_code = stan_model, data = list(N=N, x=x, Y=Y),
            chains=2, iter=10000, init= 0, 
            warmup = 10, diagnostic_file = 'diag_1.txt'
)


fit_mcmc1 <- stan2coda(fit1)
summary(fit_mcmc1)
plot(fit_mcmc1)


# 2: no init values, warmup = 5000;
fit2 <- stan(model_code = stan_model, data = list(N=N, x=x, Y=Y),
             chains=2, iter=15000, init= 0, warmup = 5000,
             diagnostic_file = 'diag_2.txt'
)

fit_mcmc2 <- stan2coda(fit2)
summary(fit_mcmc2)
plot(fit_mcmc2)

# 3: FINAL MODEL: same init values as OpenBugs, warmup = 5000;
seed <- 1234
fit3 <- stan(model_code = stan_model, data = list(N=N, x=x, Y=Y),
             chains=2, iter=15000, init= init.par, warmup = 5000, 
             diagnostic_file = 'diag_3.txt', seed = seed
)

fit_mcmc3 <- stan2coda(fit3)
fit3.summary <- summary(fit_mcmc3)
plot(fit_mcmc3)

ess <- effectiveSize(fit_mcmc3)
hpd <- HPDinterval(fit_mcmc3)
ccorr <- crosscorr(fit_mcmc3)
autocorr <- autocorr(fit_mcmc3)
mcse <- sqrt(spectrum0(fit_mcmc3)$spec / 20000)[-1]

## plots
plot(fit_mcmc3)
traceplot(fit3)
densityplot(fit_mcmc3)
xyplot(fit_mcmc3)
cumuplot(fit_mcmc3)
qqmath(fit_mcmc3)
acfplot(fit_mcmc3)
autocorr.plot(fit_mcmc3)
crosscorr.plot(fit_mcmc3)

## Diagnostics
gel.diag <- gelman.diag(fit_mcmc3)
gelman.plot(fit_mcmc3)
gewe.diag <- geweke.diag(fit_mcmc3)
geweke.plot(fit_mcmc3)
raft.diag <- raftery.diag(fit_mcmc3)
raft.diag2 <- raftery.diag(fit_mcmc3, q = 0.025, r = 0.02, s = 0.9)

save(list = ls(), file = 'stan_results_final.RData')
