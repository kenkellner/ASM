
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# For more information, see www.kenkellner.org/ASM
# ---------------------------------------------------------------

# ----------------------------------------------------
# Chapter 12  --  Overdispersion, zero inflation, and
#                 offsets in a Poisson GLM
# ----------------------------------------------------

# Last changes: 11 June 2024


# 12.1 Introduction
# -----------------
# (no code)


# -------------------
# 12.2 Overdispersion
# -------------------
# (no code)


#   12.2.1 Data generation
# ------------------------

set.seed(122)
nSites <- 50
n <- 2 * nSites
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))
eps <- rnorm(2*nSites, mean = 0, sd = 0.5)             # Normal random effect
lambda.OD <- exp(0.69 +(0.92*(as.numeric(x)-1) + eps) )
lambda.Poisson <- exp(0.69 +(0.92*(as.numeric(x)-1)) ) # For comparison

# Save true parameter values
truth <- c(Intercept = 0.69, arable = 0.92, eps_sd = 0.5)

C_OD <- rpois(n = n, lambda = lambda.OD)               # Counts with OD
C_Poisson <- rpois(n = n, lambda = lambda.Poisson)     # Counts without OD

par(mfrow = c(1,2)) # Fig. 12.1
boxplot(C_OD ~ x, col = "grey", xlab = "Land-use", main = "With overdispersion",
ylab = "Hare count", las = 1, ylim = c(0, max(C_OD)), frame = FALSE)
boxplot(C_Poisson ~ x, col = "grey", xlab = "Land-use", main = "Without overdispersion",
ylab = "Hare count", las = 1, ylim = c(0, max(C_OD)), frame = FALSE )

# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


#   12.2.2 Likelihood analysis with canned functions in R
# -------------------------------------------------------

glm.fit.no.OD <- glm(C_OD ~ x, family = poisson)
glm.fit.with.OD <- glm(C_OD ~ x, family = quasipoisson)
summary(glm.fit.no.OD)
glm_est <- c(coef(glm.fit.no.OD), NA) # Save point estimates

# Canned function for Poisson lognormal
library(PLNmodels)
PLN.fit <- PLN(matrix(C_OD, ncol = 1) ~ x)
PLN_est <- c(coef(PLN.fit), sqrt(sigma(PLN.fit)))


#   12.2.3 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Write JAGS model file
cat(file = "model12.2.3.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
tau <- pow(sigma, -2)
sigma ~ dt(0, 0.1, 1)T(0, ) # Half-Cauchy

# Likelihood
for (i in 1:n) {
  C_OD[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha + beta *x[i] + eps[i]
  eps[i] ~ dnorm(0, tau)
}
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1), sigma = rlnorm(1, 0.1))}

# Parameters to estimate
params <- c("alpha", "beta", "sigma")

# MCMC settings
na <- 2000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART <1 min), check convergence, summarize and save results
out12.2.3 <- jags(dataList, inits, params, "model12.2.3.txt",
n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out12.2.3)                # not shown
print(out12.2.3, 3)
jags_est <- as.numeric(out12.2.3$mean[1:3]) # save estimates

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, glm = glm_est, PLN = PLN_est, JAGS = jags_est)
print(comp, 4)


#   12.2.4 Bayesian analysis with NIMBLE
# --------------------------------------

library(nimble)

# Bundle and summarize data (same as before)
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Write NIMBLE model file
model12.2.4 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
sigma ~ T(dt(0, 0.1, 1), 0, )

# Likelihood
for (i in 1:n) {
  C_OD[i] ~ dpois(lambda[i]) 
  log(lambda[i]) <- alpha + beta *x[i] + eps[i]
  eps[i] ~ dnorm(0, sd = sigma)
}
} )

# Inits
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1), sigma = rlnorm(1, 0.1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma")

# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 25 sec), check convergence, summarize posteriors and save estimates
system.time( out12.2.4 <- 
    nimbleMCMC(code = model12.2.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out12.2.4) # not shown
(nsum <- nimble_summary(out12.2.4, params))   # not shown
nimble_est <- nsum[1:3,1]                     # save estimates


#   12.2.5 Bayesian analysis with Stan
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Write Stan model
cat(file = "model12_2_5.stan", "
data{
  int n;
  array[n] int C_OD;
  array[n] int x;
}

parameters{
  real alpha;
  real beta;
  real<lower=0> sigma;
  vector[n] eps;
}

model{
  vector[n] lambda;
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  for (i in 1:n){
    eps[i] ~ normal(0, sigma);
    lambda[i] = exp(alpha + beta * x[i] + eps[i]);
    C_OD[i] ~ poisson(lambda[i]);
  }
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 60/6 sec), assess convergence, print results table and save estimates
system.time(
  out12.2.5 <- stan(file = "model12_2_5.stan", data = dataList,
                    warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out12.2.5)          # not shown
print(out12.2.5, dig = 3)            # not shown
stan_est <- summary(out12.2.5)$summary[1:3,1]


#   12.2.6 Do-it-yourself MLEs
# ----------------------------

# Bundle data
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Definition of NLL for overdispersed Poisson GLM with 1 covariate
# Note this is technically a Poisson GLMM, as in Chapter 14
NLL <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  sigma <- exp(param[3])
  L <- numeric(n)
  for(i in 1:data$n){
    L[i] <- integrate(function(eps){
    lambda <- exp(alpha + beta * data$x[i] + eps)
    dpois(data$C_OD[i], lambda)* dnorm(eps, 0, sigma)},
      lower = -Inf, upper = Inf)$value
  }
  NLL <- -sum(log(L))
  NLL
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('alpha' = log(mean(C_OD)), 'beta(arable)' = 0, 'log.sigma' = 0)
out12.2.6 <- optim(inits, NLL, data = dataList, hessian = TRUE)
get_MLE(out12.2.6, 5)
diy_est <- c(out12.2.6$par[1:2], exp(out12.2.6$par[3]))


#   12.2.7 Likelihood analysis with TMB
# -------------------------------------

# Bundle and summarize data
str(dataList <- list(C_OD = C_OD, x = as.numeric(x)-1, n = length(x)))

# Write TMB model file
cat(file = "model12_2_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{


  //Describe input data
  DATA_VECTOR(C_OD);                 //response
  DATA_VECTOR(x);                    //covariate
  DATA_INTEGER(n);                   //Number of observations

  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(eps);

  Type sigma = exp(log_sigma);
  Type LL = 0.0;                     //Initialize log-likelihood at 0
  for (int i= 0; i<n; i++){
    LL += dnorm(eps(i), Type(0.0), sigma, true);
    Type lambda = exp(alpha + beta * x(i) + eps(i));

    //Calculate log-likelihood of observation and add to total
    LL += dpois(C_OD(i), lambda, true);
  }
  return -LL;                        //Return negative log likelihood
}
")

# Compile and load TMB function
compile("model12_2_7.cpp")
dyn.load(dynlib("model12_2_7"))

# Provide dimensions and starting values for parameters (incl. the random effects 'eps')
params <- list(alpha = 0, beta = 0, log_sigma = 0, eps = rep(0, dataList$n))

# Create TMB object
out12.2.7 <- MakeADFun(data = dataList,
                       parameters = params,
                       random = "eps",
                       DLL = "model12_2_7", silent = TRUE)

# Optimize TMB object, print and save results
starts <- rep(0, 3)
opt <- optim(starts, fn = out12.2.7$fn, gr = out12.2.7$gr, method = "BFGS")
(tsum <- tmb_summary(out12.2.7)) # not shown
tmb_est <- c(tsum[1:2,1], exp(tsum[3,1]))


#   12.2.8 Comparison of the parameter estimates
# ----------------------------------------------

# Compare truth with estimates from all engines and methods
comp <- cbind(truth = truth, glm = glm_est, PLN = PLN_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# -------------------
# 12.3 Zero inflation
# -------------------
# (no code)


#   12.3.1 Data generation
# ------------------------

set.seed(123)
psi <- 0.2                           # NOTE: here psi is probability of zero inflation
nSites <- 50
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))

w <- rbinom(n = 2 * nSites, size = 1, prob = 1 - psi)

lambda <- exp(0.69 +(0.92*(as.numeric(x)-1)) )

C <- rpois(n = 2*nSites, lambda = w *lambda)
data.frame('habitat' = x, 'suitability' = w, 'count' = C) # Look at data (not shown)

# Save true parameter values
truth <- c(Intercept = 0.69, arable = 0.92, psi = psi)


#   12.3.2 Likelihood analysis with canned functions in R
# -------------------------------------------------------

library(pscl)
out12.3.2 <- zeroinfl(C ~ x | 1, dist = "poisson")
summary(out12.3.2)
pscl_est <- c(coef(out12.3.2)[1:2], plogis(coef(out12.3.2)[3]))


#   12.3.3 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, n = length(x)) )

# Write JAGS model file
cat(file = "model12.3.3.txt", "
model {
# Priors
psi ~ dunif(0,1)
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
# Likelihood
for (i in 1:n) {
  w[i] ~ dbern(1-psi) # Habitat suitability
  C[i] ~ dpois(eff.lambda[i])
  eff.lambda[i] <- w[i] * lambda[i]
  log(lambda[i]) <- alpha + beta * x[i] # expected abundance at suitable sites
 }
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1), w = rep(1, 2 * nSites))}

# Parameters to estimate
params <- c("alpha", "beta", "psi")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save results
out12.3.3 <- jags(dataList, inits, params, "model12.3.3.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out12.3.3) # not shown
print(out12.3.3, 3)
jags_est <- as.numeric(out12.3.3$mean[1:3])

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, pscl = pscl_est, JAGS = jags_est)
print(comp, 4)


#   12.3.4 Bayesian analysis with NIMBLE
# --------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, n = length(x)) )

# Write NIMBLE model file
model12.3.4 <- nimbleCode( {
# Priors
psi ~ dunif(0, 1)
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:n) {
  w[i] ~ dbern(1 - psi)
  C[i] ~ dpois(eff.lambda[i]) 
  eff.lambda[i] <- w[i] * lambda[i]
  log(lambda[i]) <- alpha + beta * x[i]
}
} )

# Inits
inits <- function(){ list(alpha  =rlnorm(1), beta = rlnorm(1), w = rep(1, 2 * nSites))}

# Parameters to estimate
params <- c("alpha", "beta", "psi")

# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 25 sec), check convergence, summarize posteriors and save results
system.time( out12.3.4  <- 
    nimbleMCMC(code = model12.3.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow = c(3, 3))  ;  coda::traceplot(out12.3.4) # not shown
(nsum <- nimble_summary(out12.3.4, params))         # not shown
nimble_est <- nsum[1:3,1]


#   12.3.5 Bayesian analysis with Stan
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, n = length(x)) )

# Write Stan model
cat(file = "model12_3_5.stan", "
data{
  int n;
  array[n] int C;
  array[n] int x;
}

parameters{
  real alpha;
  real beta;
  real<lower = 0,upper = 1> psi;
}

model{
  vector[n] lambda;
  vector[n] lik;

  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);

  for (i in 1:n){
    lambda[i] = exp(alpha + beta * x[i]);

    if(C[i] == 0){ // Either w state is possible
      lik[i] = psi + (1-psi) * exp(poisson_lpmf(0 | lambda[i]));
    } else { // Only w = 1 is possible
      lik[i] = (1-psi) * exp(poisson_lpmf(C[i] | lambda[i]));
    }
    target += log(lik[i]); //target is the total log likelihood
  }
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 56/4 sec), assess convergence, print results table and save results
system.time(
  out12.3.5 <- stan(file = "model12_3_5.stan", data = dataList,
                    warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out12.3.5) # not shown
print(out12.3.5, dig = 3) # not shown
stan_est <- summary(out12.3.5)$summary[1:3,1]


#   12.3.6 Do-it-yourself MLEs
# ----------------------------

# Definition of NLL for ZIP
NLL <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  psi <- plogis(param[3])
  lambda <- exp(alpha + beta * data$x)
  lik <- numeric(data$n)

  for (i in 1:data$n){
    if(data$C[i] == 0){
      lik[i] <- psi + (1-psi) * dpois(0, lambda[i])
    } else {
      lik[i] <- (1-psi) * dpois(data$C[i], lambda[i])
    }
  }

  NLL <- -sum(log(lik))
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('Intercept' = 0, 'arable' = 0, 'logit.psi' = 0)
out12.3.6 <- optim(inits, NLL, data = dataList, hessian = TRUE)
get_MLE(out12.3.6, 4)
diy_est <- c(out12.3.6$par[1:2], plogis(out12.3.6$par[3]))


#   12.3.7 Likelihood analysis with TMB
# -------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, n = length(x)) )

# Write TMB model file
cat(file = "model12_3_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(C);                    //response
  DATA_VECTOR(x);                    //covariate
  DATA_INTEGER(n);                   //Number of observations

  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(logit_psi);

  Type psi = invlogit(logit_psi);

  Type LL = 0.0;                     //Initialize log-likelihood at 0
  for (int i=0; i< n; i++){
    Type lik;
    Type lambda = exp(alpha + beta * x(i));
    if(C(i) == 0){
      lik = psi + (1-psi) * dpois(Type(0), lambda, false);
    } else {
      lik = (1-psi) * dpois(C[i], lambda, false);
    }
    LL += log(lik);
  }
  return -LL;                        //Return negative log-likelihood
}
")

# Compile and load TMB function
compile("model12_3_7.cpp")
dyn.load(dynlib("model12_3_7"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0, logit_psi = 0)

# Create TMB object
out12.3.7 <- MakeADFun(data = dataList,
                       parameters = params,
                       DLL = "model12_3_7", silent = TRUE)

# Optimize TMB object and print results
starts <- rep(0, 3)
opt <- optim(starts, fn = out12.3.7$fn, gr = out12.3.7$gr, method = "BFGS")
(tsum <- tmb_summary(out12.3.7)) # not shown
tmb_est <- c(tsum[1:2,1], plogis(tsum[3,1]))


#   12.3.8 Comparison of the parameter estimates
# ----------------------------------------------

# Compare estimates with truth
comp <- cbind(truth = truth, pscl = pscl_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# ------------
# 12.4 Offsets
# ------------
# (no code)


#   12.4.1 Data generation
# ------------------------

set.seed(124)
nSites <- 50
A <- runif(n = 2 * nSites, 2, 5)     # Areas range in size from 2 to 5 km2
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))
linear.predictor <- log(A) + 0.69 +(0.92*(as.numeric(x)-1))
lambda <- exp(linear.predictor)
lambda <- A * exp(0.69 +(0.92*(as.numeric(x)-1))) # exactly the same !
C <- rpois(n = 2 * nSites, lambda = lambda)       # Add Poisson noise

# Save true parameter values
truth <- c(Intercept = 0.69, arable = 0.92)


#   12.4.2 Likelihood analysis with canned functions in R
# -------------------------------------------------------

glm.fit.no.offset <- glm(C ~ x, family = poisson)
glm.fit.with.offset <- glm(C ~ x, family = poisson, offset = log(A))
summary(glm.fit.no.offset)           # not shown
summary(glm.fit.with.offset)         # not shown

# Compare to truth
glmNoOff <- coef(glm.fit.no.offset)
glmOff <- coef(glm.fit.with.offset)
comp <- cbind(truth, glmNoOff = glmNoOff, glmOff = glmOff)
print(comp, 4)


#   12.4.3 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, A = A, n = length(x)))

# Write JAGS model file
cat(file = "model12.4.3.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:n) {
  C[i] ~ dpois(lambda[i])
  log(lambda[i]) <- 1 * log(A[i]) + alpha + beta *x[i] # Note offset
  # lambda[i] <- A * exp(alpha + beta *x[i])           # exactly the same
}
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out12.4.3 <- jags(dataList, inits, params, "model12.4.3.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out12.4.3) # not shown
print(out12.4.3, 3)
jags_est <- as.numeric(out12.4.3$mean[1:2])

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth, glmOff = glmOff, JAGS = jags_est)
print(comp, 4)


#   12.4.4 Bayesian analysis with NIMBLE
# --------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, A = A, n = length(x)))

# Write NIMBLE model file
model12.4.4 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:n) {
  C[i] ~ dpois(lambda[i]) 
  log(lambda[i]) <- 1 * log(A[i]) + alpha + beta *x[i]  # Note offset
  # lambda[i] <- A * exp(alpha + beta *x[i])   # exactly the same
}

} )

# Inits
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1))}

# Parameters monitored: same as before
params <- c("alpha", "beta")

# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 20 sec), check convergence, summarize posteriors and save estimates
system.time( out12.4.4 <- 
    nimbleMCMC(code = model12.4.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out12.4.4)  # not shown
(nsum <- nimble_summary(out12.4.4, params))    # not shown
nimble_est <- nsum[1:2,1]


#   12.4.5 Bayesian analysis with Stan
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, x = as.numeric(x)-1, A = A, n = length(x)))

# Write Stan model
cat(file = "model12_4_5.stan", "
data{
  int n;
  array[n] int C;
  array[n] int x;
  vector[n] A;
}

parameters{
  real alpha;
  real beta;
}

model{
  vector[n] lambda;

  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);

  for (i in 1:n){
    lambda[i] = exp(alpha + beta * x[i] + log(A[i]));
    C[i] ~ poisson(lambda[i]);
  }
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 44/2 sec), assess convergence, print results table and save estimates
system.time(
out12.4.5 <- stan(file = "model12_4_5.stan", data = dataList,
                  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out12.4.5)          # not shown
print(out12.4.5, dig = 3)            # not shown
stan_est <- summary(out12.4.5)$summary[1:2,1]


#   12.4.6 Do-it-yourself MLEs
# ----------------------------

# Bundle data into list
str(dataList <- list(C = C, x = as.numeric(x) - 1, A = A, n = length(x)))

# Definition of NLL for Poisson GLM with offset
NLL <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  mu <- exp(log(data$A) + alpha + beta * data$x)
  # mu <- data$A * exp(alpha + beta * data$x) # Exactly the same
  L <- dpois(data$C, mu)                      # Likelihood contr. for 1 observation
  LL <- log(L)                                # Loglikelihood contr. for 1 observation
  NLL <- -sum(LL)                             # NLL for all observations
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('Intercept' = 0, 'arable' = 0)
out12.4.6 <- optim(inits, NLL, data = dataList, hessian = TRUE, method = "BFGS")
get_MLE(out12.4.6, 4)
diy_est <- out12.4.6$par


#   12.4.7 Likelihood analysis with TMB 
# -------------------------------------

# Bundle data into list
str(dataList <- list(C = C, x = as.numeric(x) - 1, A = A, n = length(x)))

# Write TMB model file
cat(file = "model12_4_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(C);                    //response
  DATA_VECTOR(x);                    //covariate
  DATA_VECTOR(A);                    //offset
  DATA_INTEGER(n);                   //Number of observations

  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);

  Type LL = 0.0;                     //Initialize log-likelihood at 0

  for (int i=0; i< n; i++){
	Type mu = exp(alpha + beta * x(i) + log(A(i)));
    LL += dpois(C[i], mu, true);     //arg true means log-lik is returned
  }
  return -LL;                        //Return negative log likelihood}
}
")

# Compile and load TMB function
compile("model12_4_7.cpp")
dyn.load(dynlib("model12_4_7"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0)

# Create TMB object
out12.4.7 <- MakeADFun(data = dataList,
                       parameters = params,
                       DLL = "model12_4_7", silent = TRUE)

# Optimize TMB object, print and save results
starts <- rep(0, 2)
opt <- optim(starts, fn = out12.4.7$fn, gr = out12.4.7$gr, method = "BFGS")
(tsum <- tmb_summary(out12.4.7)) # not shown
tmb_est <- tsum[1:2,1]


#   12.4.8 Comparison of the parameter estimates
# ----------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth, glmOff = glmOff, JAGS = jags_est, NIMBLE = nimble_est,
  Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 12.5 Summary and outlook
# ------------------------
# (no code)
