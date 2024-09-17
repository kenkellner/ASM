
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -------------------------------------------------------
# Chapter 14  --  Poisson generalized linear mixed model,
#                 or Poisson GLMM
# -------------------------------------------------------

# Last changes: 11 June 2024


# 14.1 Introduction
# -----------------
# (no code)


# 14.2 Data generation 
# --------------------

set.seed(14)
nPops <- 16
nYears <- 30
n <- nPops * nYears                  # n = 480
pop <- gl(n = nPops, k = nYears)

orig.year <- rep(1:nYears, nPops)
year <- (orig.year-1)/29             # Squeeze between 0 and 1

Xmat <- model.matrix(~ pop * year - 1 - year)
print(Xmat[1:91,], 2)                # Print top 91 rows of 480

# Choose values for hyperparams and draw Normal random numbers
mu.alpha <- 3                        # Mean of intercepts
sigma.alpha <- 1                     # SD of intercepts
mu.beta <- -2                        # Mean of slopes
sigma.beta <- 0.6                    # SD of slopes
alpha <- rnorm(n = nPops, mean = mu.alpha, sd = sigma.alpha)
beta <- rnorm(n = nPops, mean = mu.beta, sd = sigma.beta)
all.effects <- c(alpha, beta)        # All together

# Save true parameter values
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta,
  sigma.alpha = sigma.alpha, sigma.beta = sigma.beta)

lin.pred <- Xmat[,] %*% all.effects  # Value of lin.predictor (eta)
C <- rpois(n = n, lambda = exp(lin.pred)) # Exponentiate and add Poisson noise
hist(C, col = "grey")                # Inspect what we’ve created (not shown)

# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB); library(lattice);
library(glmmTMB); library(DHARMa); library(abind); library(pracma)

xyplot(C ~ orig.year | pop, ylab = "Red-backed shrike counts", xlab = "Year",
  pch = 16, cex = 1.2, col = rgb(0, 0, 0, 0.4))


# 14.3 Likelihood analysis with canned functions in R
# ---------------------------------------------------

gtmb.data <- data.frame(C = C, year = year, pop = pop) # bundle data
out14.3 <- glmmTMB(C ~ year + (year || pop), data = gtmb.data, family = poisson)
summary(out14.3)                     # Inspect results
sds <- attr(VarCorr(out14.3)$cond$pop, 'stddev')       # Save results
gtmb_est <- c(fixef(out14.3)$cond, sds)

# Print estimates of random effects
ranef(out14.3) # ran-ef. estimates from zero-mean dist. (not shown)
coef(out14.3) # 'full' ran-ef estimates (not shown)

# Assess goodness-of-fit of the model (not shown)
par(mfrow = c(1, 2)) # Distribution of estimates of alpha and beta
hist(coef(out14.3)$cond$pop[,1], main = 'alpha', breaks = 12)
hist(coef(out14.3)$cond$pop[,2], main = 'beta', breaks = 12)
library(DHARMa)
simOut <- simulateResiduals(out14.3, n = 1000, plot = TRUE)


# 14.4 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# Write JAGS model file
cat(file = "model14.4.txt", "
model {
# Priors
for (i in 1:nPops){
  # Models for the sets of random effects alpha and beta
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # Intercepts
  beta[i] ~ dnorm(mu.beta, tau.beta)    # Slopes
}
mu.alpha ~ dnorm(0, 0.0001)             # Hyperparam. for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 0.1, 1)T(0,)        # half-Cauchy prior
mu.beta ~ dnorm(0, 0.0001)              # Hyperparam. for random slopes
tau.beta <- pow(sigma.beta, -2)
sigma.beta ~ dt(0, 0.1, 1)T(0,)         # half-Cauchy prior

# 'Likelihood'
for (i in 1:n){
  C[i] ~ dpois(lambda[i])
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]] * year[i])
  # log(lambda[i]) <- alpha[pop[i]] + beta[pop[i]]* year[i] # same
}
# Posterior predictive simulations (used in DHARMa GoF assessment)
for (i in 1:n) {
  C.new[i] ~ dpois(lambda[i])
}
}
")

# Function to generate starting values
inits <- function(){
list(mu.alpha = rnorm(1), mu.beta = rnorm(1),
  sigma.alpha = runif(1), sigma.beta = runif(1))
}

# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta", "alpha", "beta", "C.new")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out14.4 <- jags(dataList, inits, params, "model14.4.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out14.4) # not shown
print(out14.4, 3) # partially shown

# Distribution of random effects estimates (alpha, beta)
par(mfrow = c(1,2))                  # Very low power with small samples. . .
hist(out14.4$mean$alpha, main = 'alpha', breaks = 12)
hist(out14.4$mean$beta, main = 'beta', breaks = 12)


# Do quantile residual assessments (not shown)
C.new <- out14.4$sims.list$C.new
sim <- createDHARMa(simulatedResponse = t(C.new), observedResponse = C,
fittedPredictedResponse = apply(C.new, 2, median), integerResponse = T)
plot(sim)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- unlist(out14.4$mean[1:4])
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est)
print(comp, 4)


# 14.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# Write NIMBLE model file
model14.5 <- nimbleCode( {
  # Priors
  for (i in 1:nPops){    # Priors for random effects
    alpha[i] ~ dnorm(mu.alpha, sd = sigma.alpha) # Intercepts
    beta[i] ~ dnorm(mu.beta, sd = sigma.beta)    # Slopes
  }
  mu.alpha ~ dnorm(0, sd = 100)                  # Same as in JAGS
  sigma.alpha ~ T(dt(0, 0.1, 1), 0,)

  mu.beta ~ dnorm(0, sd = 100)
  sigma.beta ~ T(dt(0, 0.1, 1), 0, )

  # 'Likelihood'
  for (i in 1:n){
    C[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]] * year[i])
  }
} )


# Inits
inits <- function(){
  list(mu.alpha = rnorm(1), mu.beta = rnorm(1),
       sigma.alpha = runif(1), sigma.beta = runif(1))
}

# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta", "alpha", "beta")

# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 50 sec), check convergence, summarize posteriors and save estimates
system.time( out14.5 <- 
    nimbleMCMC(code = model14.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out14.5)  # not shown
(nsum <- nimble_summary(out14.5, params))    # not shown
nimble_est <- nsum[1:4,1]                    # Save estimates


# 14.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, pop = as.numeric(pop), year = year, nPops = nPops, n = n) )

# Write Stan model
cat(file = "model14_6.stan", "
data{
  int n;
  int nPops;
  array[n] int C;
  vector[n] year;
  array[n] int pop;
}

parameters{
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[nPops] alpha;
  vector[nPops] beta;
}

model{
  vector[n] lambda;
  mu_alpha ~ normal(0, 100);
  sigma_alpha ~ cauchy(0, sqrt(10));
  mu_beta ~ normal(0, 100);
  sigma_beta ~ cauchy(0, sqrt(10));
  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
    beta[i] ~ normal(mu_beta, sigma_beta);
  }
  for (i in 1:n){
    lambda[i] = exp(alpha[pop[i]] + beta[pop[i]] * year[i]);
    C[i] ~ poisson(lambda[i]);
  }
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 67/12 sec), assess convergence, show results and save estimates
system.time(
out14.6 <- stan(file = "model14_6.stan", data = dataList,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out14.6) # not shown
print(out14.6, dig = 3) # not shown
stan_est <- summary(out14.6)$summary[1:4,1] # Save estimates


# 14.7 Do-it-yourself MLEs
# ------------------------

# Definition of custom functions f and g (note f is needed in g,
#   and g is needed in definition of NLL function below)
f <- function(alpha_j, beta_j, mu.alpha, sig.alpha, 
  mu.beta, sig.beta, C_j, year_j){
    lambda <- exp(alpha_j + beta_j * year_j)
    prod(dpois(C_j, lambda)) * dnorm(alpha_j, mu.alpha, sig.alpha) *
      dnorm(beta_j, mu.beta, sig.beta)
}

g <- function(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, year_j){
  integral2(fun = f, xmin = -10, xmax = 10, ymin = -10, ymax = 10,
    vectorized = FALSE,
    mu.alpha = mu.alpha, sig.alpha = sig.alpha, mu.beta = mu.beta,
    sig.beta = sig.beta, C_j = C_j, year_j = year_j)$Q
}

# Define NLL for Poisson GLMM (using custom functions f and g)
NLL <- function(pars, data) {
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sig.alpha <- exp(pars[3])
  sig.beta <- exp(pars[4])
  nll <- 0 # Initialize nll at 0 prior to summation in loop below
  for (j in 1:data$nPops){
    # Subset data to just pop j
    C_j <- data$C[data$pop == j]
    year_j <- data$year[data$pop == j]
    lik <- g(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, year_j)
    nll <- nll - log(lik)
  }
  return(nll)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates (ART 122 sec)
inits <- c(mu.alpha = 0, mu.beta = 0, log.sigma.alpha = 0, log.sigma.beta = 0)
system.time(
  out14.7 <- optim(inits, NLL, data = dataList, hessian = TRUE,
    method = 'BFGS', control = list(trace = 1, REPORT = 1)) )
get_MLE(out14.7, 4)
diy_est <- c(out14.7$par[1:2], exp(out14.7$par[3:4])) # Save estimates


# 14.8 Likelihood analysis with TMB
# ---------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop - 1 #indices start at 0 in TMB
str(tmbData)

# Write TMB model file
cat(file = "model14_8.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(C); //response
  DATA_VECTOR(year); //covariate
  DATA_IVECTOR(pop); //population index (of integers; IVECTOR)
  DATA_INTEGER(n); //Number of observations
  DATA_INTEGER(nPops); //Number of groups
 
 //Describe parameters
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER(log_sigma_alpha);
  PARAMETER(log_sigma_beta);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type sigma_alpha = exp(log_sigma_alpha);
  Type sigma_beta = exp(log_sigma_beta);

  Type LL = 0.0; //Initialize log-likelihood at 0

  //Random effects
  for (int i = 0; i<nPops; i++){
    LL += dnorm(alpha(i), mu_alpha, sigma_alpha, true);
    LL += dnorm(beta(i), mu_beta, sigma_beta, true);
  }

  for (int i = 0; i<n; i++){
    Type lambda = exp(alpha(pop(i)) + beta(pop(i)) * year(i));
    LL += dpois(C(i), lambda, true);
  }
  return -LL;
}
")

# Compile and load TMB function
compile("model14_8.cpp")
dyn.load(dynlib("model14_8"))

# Provide dimensions and starting values for parameters
params <- list(mu_alpha = 0, mu_beta = 0, log_sigma_alpha = 0,
  log_sigma_beta = 0, alpha = rep(0, tmbData$nPops),
  beta = rep(0, tmbData$nPops))

# Create TMB object
out14.8 <- MakeADFun(data = tmbData,
  parameters = params,
  random = c("alpha", "beta"),
  DLL = "model14_8", silent = TRUE)

# Optimize TMB object and print and save results
opt <- optim(out14.8$par, fn = out14.8$fn, gr = out14.8$gr, method = "BFGS")
(tsum <- tmb_summary(out14.8)) # not shown
tmb_est <- c(tsum[1:2,1], exp(tsum[3:4,1])) # save estimates


# 14.9 Comparison of the parameter estimates
# ------------------------------------------

# Compare estimates with truth
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 14.10 Summary
# -------------
# (no code)
