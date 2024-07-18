
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# --------------------------------------------------------
# Chapter 17  --  Binomial generalized linear mixed model
# --------------------------------------------------------

# Last changes: 11 June 2024


# 17.1 Introduction
# -----------------
# (no code)


# 17.2 Data generation
# --------------------

set.seed(17)
nPops <- 16
nYears <- 10
n <- nPops * nYears # n = 160
pop <- gl(n = nPops, k = nYears)

precip <- runif(n, -1, 1)

N <- round(runif(n, 20, 50) )

Xmat <- model.matrix(~ pop * precip - 1 - precip)
print(Xmat[1:91,], dig = 2) # Print top 91 rows (not shown)

mu.alpha <- 0 # Select hyperparams
mu.beta <- -2
sigma.alpha <- 1
sigma.beta <- 1
alpha <- rnorm(n = nPops, mean = mu.alpha, sd = sigma.alpha)
beta <- rnorm(n = nPops, mean = mu.beta, sd = sigma.beta)
all.pars <- c(alpha, beta) # All parameters together

# Save vector of true parameter values
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta,
sigma.alpha = sigma.alpha, sigma.beta = sigma.beta)

lin.pred <- Xmat %*% all.pars # Value of lin.predictor
exp.p <- exp(lin.pred)/(1 + exp(lin.pred)) # Expected proportion

library(lattice)
xyplot(exp.p ~ precip | pop, ylab = "Expected woodchat shrike breeding success ",
  xlab = "Spring precipitation index", main = "Expected breeding success", pch = 16,
  cex = 1.2, col = rgb(0, 0, 0, 0.4))

C <- rbinom(n = n, size = N, prob = exp.p) # Add binomial variation
xyplot(C/N ~ precip | pop, ylab = "Realized woodchat shrike breeding success ",
  xlab = "Spring precipitation index", main = "Realized breeding success",
  pch = 16, cex = 1.2, col = rgb(0, 0, 0, 0.4))

# Required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 17.3 Likelihood analysis with canned functions in R
# ---------------------------------------------------

# Fit model, present and save estimates
library(glmmTMB)
gtmb.data <- data.frame(C = C, N = N, precip = precip, pop = pop) # bundle data
out17.3 <- glmmTMB(cbind(C, N-C) ~ precip + (precip || pop), data = gtmb.data, family = binomial)
summary(out17.3) # Inspect results
sds <- attr(VarCorr(out17.3)$cond$pop, 'stddev') # Save results
gtmb_est <- c(fixef(out17.3)$cond, sds)

# Model goodness of fit
par(mfrow = c(1, 2)) # Distribution of estimates of alpha and beta
hist(ranef(out17.3)$cond$pop[,1], main = 'alpha', breaks = 12)
hist(ranef(out17.3)$cond$pop[,2], main = 'beta', breaks = 12)

library(DHARMa)
simOut <- simulateResiduals(out17.3, n = 1000, plot = TRUE)


# 17.4 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, N = N, pop = as.numeric(pop), precip = precip, nPops = nPops, n = n) )

# Write JAGS model file
cat(file = "model17.4.txt", "
model {
# Priors
# Models for the random effects alpha and beta
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # Intercepts
  beta[i] ~ dnorm(mu.beta, tau.beta) # Slopes
}
# (Hyper-)Priors for hyperparameters
mu.alpha ~ dnorm(0, 0.0001) # Hyperparameter for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 0.1, 1)I(0,)

mu.beta ~ dnorm(0, 0.0001) # Hyperparameter for random slopes
tau.beta <- pow(sigma.beta, -2)
sigma.beta ~ dt(0, 0.1, 1)I(0, )

# Likelihood
for (i in 1:n) {
  C[i] ~ dbin(p[i], N[i])
  logit(p[i]) <- alpha[pop[i]] + beta[pop[i]]* precip[i]
}

# Posterior predictive simulations (used in DHARMa GoF assessment)
for (i in 1:n) {
  C.new[i] ~ dbin(p[i], N[i])
}
}
")

# Function to generate starting values
inits <- function(){ list(mu.alpha = rnorm(1, 0, 1), mu.beta = rnorm(1, 0, 1))}

# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
"alpha", "beta", "C.new")

# MCMC settings
na <- 1000 ; ni <- 12000 ; nb <- 2000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out17.4 <- jags(dataList, inits, params, "model17.4.txt", n.iter = ni, n.burnin = nb,
n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out17.4) # not shown
print(out17.4, 3)

# Distribution of random effects estimates (alpha, beta)
par(mfrow = c(1,2)) # Low power with small samples. . .
hist(out17.4$mean$alpha, main = 'alpha', breaks = 12)
hist(out17.4$mean$beta, main = 'beta', breaks = 12)

# Do quantile residual assessments (not shown)
#library(DHARMa)
#library(abind)
C.new <- out17.4$sims.list$C.new
sim <- createDHARMa(simulatedResponse = t(C.new), observedResponse = C,
fittedPredictedResponse = apply(C.new, 2, median), integerResponse = T)
plot(sim)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- out17.4$summary[1:4,1]
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est)
print(comp, 4)


# 17.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(C = C, N = N, pop = as.numeric(pop), precip = precip, nPops = nPops, n = n) )

# Write NIMBLE model file
model17.5 <- nimbleCode( {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, sd = sigma.alpha)   # Intercepts
  beta[i] ~ dnorm(mu.beta, sd = sigma.beta)      # Slopes
}
mu.alpha ~ dnorm(0, sd = 1000) # Hyperparam. for random intercepts

sigma.alpha ~ T(dt(0, 0.1, 1), 0, 100)

mu.beta ~ dnorm(0, sd = 1000)     # Hyperparameter for random slopes

sigma.beta ~ T(dt(0, 0.1, 1), 0, 100)

# 'Likelihood'
for (i in 1:n) {
  C[i] ~ dbin(p[i], N[i])
  logit(p[i]) <- alpha[pop[i]] + beta[pop[i]]* precip[i]
}
} )

# Inits
inits <- function(){
  list(mu.alpha = rnorm(1, 0, 1), mu.beta = rnorm(1, 0, 1),
       sigma.alpha = rlnorm(1), sigma.beta = rlnorm(1))
}

# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
            "alpha", "beta")

# MCMC settings
ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 36 sec), check convergence, summarize posteriors and save estimates
system.time( out17.5 <- 
    nimbleMCMC(code = model17.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out17.5) # not shown
(nsum <- nimble_summary(out17.5, params))   # not shown
nimble_est <- nsum[1:4, 1]                  # Save estimates


# 17.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, N = N, pop = as.numeric(pop), precip = precip, nPops = nPops, n = n) )

# Write Stan model
cat(file = "model17_6.stan", "
data{
  int n; //Number of samples
  int nPops; //Number of populations
  array[n] int N; //Number of trials in each sample
  array[n] int C; //Successes in each sample
  vector[n] precip; //covariate
  array[n] int pop; //Population index
}

parameters{
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_beta;
  vector[nPops] alpha;
  vector[nPops] beta;
}

model{
  vector[n] p; //Estimated success probability
  mu_alpha ~ normal(0, 100);
  mu_beta ~ normal(0, 100);
  sigma_alpha ~ cauchy(0, sqrt(10));
  sigma_beta ~ cauchy(0, sqrt(10));
  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
    beta[i] ~ normal(mu_beta, sigma_beta);
  }
  
  for(i in 1:n){
    p[i] = inv_logit(alpha[pop[i]] + beta[pop[i]] * precip[i]);
    C[i] ~ binomial(N[i], p[i]);
  }
}
")

# Parameters to estimate
params <- c("mu_alpha", "mu_beta", "sigma_alpha", "sigma_beta",
            "alpha", "beta")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 50/5 sec), assess convergence, print and save results
system.time(
out17.6 <- stan(file = "model17_6.stan", data = dataList, pars = params,
                warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out17.6) # not shown
print(out17.6, dig = 3) # not shown
stan_est <- summary(out17.6)$summary[1:4,1] # Save estimates


# 17.7 Do-it-yourself MLEs
# ------------------------

f <- function(alpha_j, beta_j, mu.alpha, sig.alpha,
  mu.beta, sig.beta, C_j, N_j, precip_j){
  p <- plogis(alpha_j + beta_j * precip_j)
  prod(dbinom(C_j, N_j, p)) * dnorm(alpha_j, mu.alpha, sig.alpha) *
  dnorm(beta_j, mu.beta, sig.beta)
}

library(pracma)
g <- function(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, N_j, precip_j){
tryCatch({
  integral2(fun = f, xmin = -10, xmax = 10, ymin = -10, ymax = 10,
    vectorized = FALSE,
    mu.alpha = mu.alpha, sig.alpha = sig.alpha, mu.beta = mu.beta,
    sig.beta = sig.beta, C_j = C_j, N_j = N_j, precip_j = precip_j)$Q
  }, error = function(e) return(Inf))
}

# Define NLL for Binomial GLMM
NLL <- function(pars, data) {
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sig.alpha <- exp(pars[3])
  sig.beta <- exp(pars[4])
  
  nll <- 0 # Initialize at zero
  
  for (j in 1:data$nPops){
    #Subset data to just pop j
    C_j <- data$C[data$pop == j]
    precip_j <- data$precip[data$pop == j]
    N_j <- data$N[data$pop == j]
    lik <- g(mu.alpha, sig.alpha, mu.beta, sig.beta, C_j, N_j, precip_j)
    nll <- nll - log(lik)
  }
  return(nll)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates (ART 100 sec)
inits <- c(mu.alpha = 0, mu.beta = 0, log.sigma.alpha = 0, log.sigma.beta = 0)
system.time(
  out17.7 <- optim(inits, NLL, data = dataList, hessian = TRUE,
    method = 'BFGS', control = list(trace = 1, REPORT = 1)) )
get_MLE(out17.7, 4)
diy_est <- c(out17.7$par[1:2], exp(out17.7$par[3:4])) # Save estimates


# 17.8 Likelihood analysis with TMB
# ---------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop - 1 #indices start at 0 in TMB
str(tmbData)

# Write TMB model file
cat(file = "model17_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(C);
  DATA_VECTOR(N);
  DATA_IVECTOR(pop);
  DATA_VECTOR(precip);
  DATA_INTEGER(nPops);
  DATA_INTEGER(n);
  
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
  
  for (int j = 0; j < nPops; j++){
    LL += dnorm(alpha(j), mu_alpha, sigma_alpha, true);
    LL += dnorm(beta(j), mu_beta, sigma_beta, true);
  }

  vector <Type> p(n);
  for (int i = 0; i < n; i++){
    p(i) = invlogit(alpha(pop(i)) + beta(pop(i)) * precip(i));
    LL += dbinom(C(i), N(i), p(i), true);
  }

  return -LL;
}
")

# Compile and load TMB function
compile("model17_8.cpp")
dyn.load(dynlib("model17_8"))

# Provide dimensions and starting values for parameters
params <- list(mu_alpha = 0, mu_beta = 0, log_sigma_alpha = 0,
              log_sigma_beta = 0, alpha = rep(0, tmbData$nPops),
              beta = rep(0, tmbData$nPops))

# Create TMB object
out17.8 <- MakeADFun(data = tmbData, parameters = params,
              random = c("alpha", "beta"),
              DLL = "model17_8", silent = TRUE)

# Optimize TMB object and print and save results
opt <- optim(out17.8$par, fn = out17.8$fn, gr = out17.8$gr, method = "BFGS")
(tsum <- tmb_summary(out17.8))
tmb_est <- c(tsum[1:2,1], exp(tsum[3:4,1])) # Save estimates


# 17.9 Comparison of the parameter estimates
# ------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, gtmb = gtmb_est, JAGS = jags_est,
              NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 17.10 Summary
# -------------
# (no code)
