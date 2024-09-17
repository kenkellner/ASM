
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ----------------------------------------------------------------------
# Chapter 7  --  Models with a single categoricalcovariate with more 
#                than two levels
# ----------------------------------------------------------------------

# Last changes: 11 June 2024


# 7.1 Introduction: fixed and random effects
# ------------------------------------------
# (no code)


# 7.2 Fixed-effects models
# ------------------------
# (no code)


#   7.2.1 Data generation
# -----------------------

# Simulate a data set
set.seed(72) # Initialize RNGs
nPops <- 5 # Number of populations
nSample <- 10 # Number of snakes in each
pop.means <- c(50, 40, 45, 55, 60) # Population mean SVL
sigma <- 5 # Residual sd
n <- nPops * nSample # Total number of data points
eps <- rnorm(n, 0, sigma) # Residuals
pop <- factor(rep(1:5, rep(nSample, nPops))) # Indicator for population
means <- rep(pop.means, rep(nSample, nPops))
X <- as.matrix(model.matrix(~ pop-1)) # Create design matrix
X # Inspect design matrix
y <- as.numeric(X %*% as.matrix(pop.means) + eps)
# %*% denotes matrix multiplication

# Save true values for later comparisons
truth <- c(pop.means, sigma)
names(truth) <- c(paste0("pop", 1:5), "sigma")

# Make plot (Fig. 7–2)
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ pop, col = "grey", xlab = "Population", ylab = "SVL", main = "", las = 1, frame = FALSE)

# Load required libraries
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)


#   7.2.2 Likelihood analysis with canned functions in R
# ------------------------------------------------------

# Default treatment contrast, or effects, parameterization (not shown)
summary(out <- lm(y ~ pop))
# Means parameterization of factor levels
summary(out72.2 <- lm(y ~ pop - 1))

# Save least-squares estimates
lm_est <- c(coef(out72.2), sigma = sigma(out72.2))

# Check goodness-of-fit using traditional and quantile residuals
plot(out72.2) # Traditional residual check for comparison
simOut <- simulateResiduals(out72.2, n = 1000, plot = TRUE)


#   7.2.3 Bayesian analysis with JAGS
# -----------------------------------

# Bundle and summarize data
# Note that JAGS requires us to convert pop from a factor to numeric
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# Write JAGS model file
cat(file = "model72.3.txt", "
model {
# Priors
for (i in 1:nPops){ # Define alpha as a vector
  alpha[i] ~ dnorm(0, 0.0001)
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)

# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mean[i], tau)
  mean[i] <- alpha[pop[i]]
}

# Derived quantities (effects estimates)
effe2 <- alpha[2] - alpha[1]
effe3 <- alpha[3] - alpha[1]
effe4 <- alpha[4] - alpha[1]
effe5 <- alpha[5] - alpha[1]

# Custom hypothesis test/Define your own contrasts
test1 <- (effe2 + effe3) - (effe4 + effe5) # Equals 0 when 2 + 3 = 4 + 5
test2 <- effe5 - 2 * effe4 # Equals 0 when effe5 = 2*effe4
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(nPops, mean = mean(y)), sigma = rlnorm(1) )}

# Parameters to estimate
params <- c("alpha", "sigma", "effe2", "effe3", "effe4", "effe5", "test1", "test2")
# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 sec), check convergence and summarize posteriors
out72.3 <- jags(dataList, inits, params, "model72.3.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out72.3) # not shown
print(out72.3, 3)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- out72.3$summary[1:6, 1]
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


#   7.2.4 Bayesian analysis with NIMBLE
#--------------------------------------

library(nimble)

# Bundle and summarize data (same as before)
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# Write NIMBLE model file
model72.4 <- nimbleCode( {
# Priors
for (i in 1:nPops){                 # Define alpha as a vector
  alpha[i] ~ dnorm(0, 0.001)
}
sigma ~ dunif(0, 10)

# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mean[i], sd = sigma) # Note use of sigma
  mean[i] <- alpha[pop[i]]
}

# Derived quantities
effe2 <- alpha[2] - alpha[1]
effe3 <- alpha[3] - alpha[1]
effe4 <- alpha[4] - alpha[1]
effe5 <- alpha[5] - alpha[1]

# Custom hypothesis test / Define your own contrasts
test1 <- (effe2+effe3) - (effe4+effe5) # Equals 0 when 2+3 = 4+5
test2 <- effe5 - 2 * effe4             # Equals 0 when effe5 = 2*effe4
}
)

# Inits
inits <- function(){ list(alpha = rnorm(nPops, mean = mean(y)), sigma = rlnorm(1) )}

# Parameters monitored: same as before
params <- c("alpha", "sigma", "effe2", "effe3", "effe4", "effe5", "test1", "test2")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 24 sec), check convergence, summarize posteriors and save estimates
system.time(
  out72.4 <- nimbleMCMC(code = model72.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out72.4)     # not shown
(nsum <- nimble_summary(out72.4, params))       # not shown
nimble_est <- nsum[1:6,1]

# 7.2.5 Bayesian analysis with Stan

# Bundle and summarize data (same as for NIMBLE)
str(dataList <- list(y = y, pop = as.numeric(pop), n = n, nPops = nPops))

# Write text file with model description in BUGS language
cat(file = "model72_5.stan",
"
data {
  int <lower = 1> n; // Declare all data
  int <lower = 1> nPops;
  vector[n] y;
  array[n] int <lower = 1> pop;
}

parameters { // Define format for all parameters
  vector [nPops] alpha; // Mean value for each pop
  real <lower = 0> sigma; // Standard deviation must be positive
}

model {
  // Priors
  alpha ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  // Likelihood
  for(i in 1:n) {
    y[i] ~ normal(alpha[pop[i]], sigma);
  }
}

generated quantities {
  // Derived quantities
  real effe2 = alpha[2] - alpha[1];
  real effe3 = alpha[3] - alpha[1];
  real effe4 = alpha[4] - alpha[1];
  real effe5 = alpha[5] - alpha[1];
  // Custom hypothesis test / Define your own contrasts
  // test1 equals 0 when 2 + 3 = 4 + 5
  real test1 = (effe2 + effe3) - (effe4 + effe5);
  real test2 = effe5 - 2 * effe4;
}
" )

# HMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# Call STAN (ART 46/2 sec)
system.time(
out72.5 <- stan(file = "model72_5.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out72.5) # not shown
print(out72.5, dig = 2) # not shown
stan_est <- summary(out72.5)$summary[1:6,1]


#   7.2.6 Do-it-yourself maximum likelihood estimates
# ---------------------------------------------------

# Definition of NLL for a one-factor linear model with Gaussian errors
NLL <- function(param, y, X) {
  alpha <- param[1:5] # Population means
  sigma <- param[6] # Residual SD
  mu <- X %*% alpha # Get pop mean for each datum
  LL <- dnorm(y, mu, sigma, log = TRUE) # log-likelihood each datum
  NLL <- -sum(LL) # NLL for all data points
  return(NLL)
}

# Get desired design matrix (means parameterization)
X <- model.matrix(~ pop - 1)

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('mu1' = 50, 'mu2' = 50, 'mu3' = 50, 'mu4' = 50, 'mu5' = 50, 'sigma' = 10)
out72.6 <- optim(inits, NLL, y = y, X = X, hessian = TRUE)
get_MLE(out72.6, 4)
diy_est <- out72.6$par


#   7.2.7 Likelihood analysis with TMB
# ------------------------------------

# Bundle and summarize data (similar as before, except for pop)
str(tmbData <- list(y = y, pop = as.numeric(pop) - 1, n = n))

# Write TMB model file
cat(file = "model72_7.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(y); //response
  DATA_IVECTOR(pop); //Population index: note IVECTOR, not VECTOR
  DATA_INTEGER(n); //Number of obs

  //Describe parameters
  PARAMETER_VECTOR(alpha); //Population means
  PARAMETER(log_sigma); //log(residual standard deviation)
  Type sigma = exp(log_sigma); //Standard deviation
  Type LL = 0; //Initialize total log likelihood at 0

  // Likelihood of each datapoint
  for (int i = 0; i < n; i++){ //Note index starts at 0 instead of 1!
    LL += dnorm(y(i), alpha(pop(i)), sigma, true);
  }

  // Derived effects (note index starts at 0!)
  Type effe2 = alpha(1) - alpha(0);
  Type effe3 = alpha(2) - alpha(0);
  Type effe4 = alpha(3) - alpha(0);
  Type effe5 = alpha(4) - alpha(0);
  
  // Custom hypothesis test
  Type test1 = (effe2 + effe3) - (effe4 + effe5);
  Type test2 = effe5 - 2 * effe4;
  
  // Tell TMB to save derived parameter values
  ADREPORT(effe2);
  ADREPORT(effe3);
  ADREPORT(effe4);
  ADREPORT(effe5);
  ADREPORT(test1);
  ADREPORT(test2);
  return -LL; //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model72_7.cpp") # gibberish....
dyn.load(dynlib("model72_7"))

# Provide dimensions and starting values for parameters
params <- list(alpha = rep(0, nPops), log_sigma = 0)

# Create TMB object
out72.7 <- MakeADFun(data = tmbData, parameters = params,
                     DLL = "model72_7", silent = TRUE)

# Optimize TMB object, print results and save estimates
opt <- optim(out72.7$par, fn = out72.7$fn, gr = out72.7$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out72.7))
tmb_est <- c(opt$par[1:5], exp(opt$par[6]))


#   7.2.8 Comparison of the parameter estimates
# ---------------------------------------------

# Compare results with truth
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est,
NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 7.3 Random-effects model
# ------------------------
# (no code)


#   7.3.1 Data generation
# -----------------------

# Simulate a data set
set.seed(73)
nPops <- 10 # Number of populations: choose 10 rather than 5
nSample <- 12 # Number of snakes in each
n <- nPops * nSample # Total number of data points

pop.grand.mean <- 50 # Grand mean SVL
pop.sd <- 3 # sd of population effects about mean
pop.means <- rnorm(n = nPops, mean = pop.grand.mean, sd = pop.sd)
sigma <- 5 # Residual sd
eps <- rnorm(n, 0, sigma) # Draw residuals

pop <- factor(rep(1:nPops, rep(nSample, nPops)))
Xmat <- as.matrix(model.matrix(~ pop - 1))
y <- as.numeric(Xmat %*% as.matrix(pop.means) + eps)

# Save true values for later comparisons
truth <- c(pop.grand.mean = pop.grand.mean, pop.sd = pop.sd, residual.sd = sigma)

# Make a plot (Fig. 7–3)
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ pop, col = "grey", xlab = "Population", ylab = "SVL", main = "", las = 1, frame = FALSE)
abline(h = pop.grand.mean)


#   7.3.2 Likelihood analysis with canned functions in R
# ------------------------------------------------------

library(glmmTMB) # Load glmmTMB

# Fit model and inspect results
summary(out73.2ML <- glmmTMB(y ~ 1 + 1 | pop, REML = FALSE)) # ML
summary(out73.2REML <- glmmTMB(y ~ 1 + 1 | pop, REML = TRUE)) # REML
ranef(out73.2REML) # Look: THESE are the estimated random effects !

summary(out73.2ML <- glmmTMB(y ~ 1 + 1 | pop, REML = FALSE)) # ML
ranef(out73.2REML) # Look: THESE are the estimated random effects !

# Save estimates from glmmTMB and make interim comparison
gtmb_est <- c(fixef(out73.2ML)$cond, sqrt(VarCorr(out73.2ML)$cond$pop[1]),
sigma(out73.2ML))
gtmb_REMLest <- c(fixef(out73.2REML)$cond, sqrt(VarCorr(out73.2REML)$cond$pop[1]),
sigma(out73.2REML))
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest)
print(comp, 4) # not shown, but note REML estimate of pop.sd greater

# Check goodness-of-fit using quantile residuals
simOut <- simulateResiduals(out73.2ML, n = 1000, plot = TRUE)


#   7.3.3 Bayesian analysis with JAGS ... and a quick fixed-random comparison
# ---------------------------------------------------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# Write JAGS model file
cat(file = "model73.3.txt", "
model {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # Prior for population means
  pop.mean[i] ~ dnorm(pop.grand.mean, pop.tau) # Random effects !
  # Pop. effects as derived quantities
  effe[i] <- pop.mean[i] - pop.grand.mean # Ranef as in lme4
}
pop.grand.mean ~ dnorm(0, 0.0001) # Hyperprior for grand mean svl
pop.sd ~ dt(0, 0.01, 1)T(0,) # Hyperprior for sd of population effects
residual.sd ~ dunif(0, 10) # Prior for residual sd

pop.tau <- pow(pop.sd, -2)
residual.tau <- pow(residual.sd, -2)

# 'Likelihood'
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], residual.tau)
  mu[i] <- pop.mean[pop[i]]
}
}
")

# Function to generate starting values
inits <- function(){ list(pop.grand.mean = runif(1, 1, 100),
  pop.sd = runif(1), residual.sd = runif(1)) }

# Parameters to estimate
params <- c("pop.grand.mean", "pop.sd", "residual.sd", "pop.mean", "effe")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART 1 sec), check convergence, summarize posteriors and save estimates
out73.3 <- jags(dataList, inits, params, "model73.3.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI:: traceplot(out73.3) # not shown
print(out73.3, 3)
jags_est <- out73.3$summary[1:3,1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest, JAGS = jags_est)
print(comp, 4)

# Check shape of posterior distribution for pop.sd (Fig. 7.4)
hist(out73.3$sims.list$pop.sd, col = 'grey', breaks = 30, main = '', freq = F)

# Compute the posterior mode of pop.sd as a better point estimate
library(MCMCglmm)
posterior.mode(mcmc(out73.3$sims.list$pop.sd))

# Write JAGS model file
cat(file = "model73.3fix.txt", "
model {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # Prior for population means: no estimated hyperparameters now !
  pop.mean[i] ~ dnorm(0, 0.0001)
}
residual.sd ~ dunif(0, 10) # Prior for residual sd
residual.tau <- pow(residual.sd, -2)

# 'Likelihood'
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], residual.tau)
  mu[i] <- pop.mean[pop[i]]
}
}
")

# Function to generate starting values
inits <- function(){ list(pop.mean = runif(nPops, 30, 70), residual.sd = runif(1)) }

# Parameters to estimate
params <- c("pop.mean", "residual.sd")

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out73.3fix <- jags(dataList, inits, params, "model73.3fix.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out73.3fix) # not shown
print(out73.3fix, 3) # also not shown

# Compare fixed- and the random-effects population means (Fig. 7.5)
off <- 0.15
plot((1:nPops)-off, pop.means, pch = 16, cex = 2, col = 'red', frame = FALSE,
  xlab = 'Population', ylab = 'Mean SVL in population', las = 1, ylim = c(38, 62))
abline(h = out73.3$mean$pop.grand.mean, lty = 2, lwd = 3, col = 'blue')
points(1:nPops, out73.3$mean$pop.mean, pch = 16, cex = 2, col = 'blue')
segments(1:nPops, out73.3$q2.5$pop.mean, 1:nPops, out73.3$q97.5$pop.mean, lwd = 2, col = 'blue')
points((1:nPops) + off, out73.3fix$mean$pop.mean, pch = 16, cex = 2, col = 'black')
segments((1:nPops) + off, out73.3fix$q2.5$pop.mean, (1:nPops) + off,
  out73.3fix$q97.5$pop.mean, lwd = 2, col = 'black')
legend('topleft', pch = 16, cex = 1.5, col = c("red", "blue", "black"), legend = c
("Truth", "Random effects", "Fixed effects"), bty = 'n')


#   7.3.4 Bayesian analysis with NIMBLE
# -------------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# Write NIMBLE model file
model73.4 <- nimbleCode( {
# Priors (and also some derived quantities)
for (i in 1:nPops){
  # Prior for population means
  pop.mean[i] ~ dnorm(pop.grand.mean, sd = pop.sd)
  # Pop. effects as derived quantities
  effe[i] <- pop.mean[i] - pop.grand.mean
}
pop.grand.mean ~ dnorm(0, sd = 1000) # Hyperprior for grand mean svl
pop.sd ~ T(dt(0, 0.01, 1), 0,) # Hyperprior for sd of population effects
residual.sd ~ dunif(0, 10) # Prior for residual sd

# 'Likelihood'
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], sd = residual.sd)
  mu[i] <- pop.mean[pop[i]]
}
}
)

# Inits
inits <- function(){ list(pop.grand.mean = runif(1, 0, 100), 
  pop.sd = runif(1), residual.sd = runif(1)) }

# Parameters monitored: same as before
params <- c("pop.grand.mean","pop.sd","residual.sd","pop.mean","effe")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 32 sec), check convergence, summarize posteriors and save estimates
system.time(
  out73.4 <- nimbleMCMC(code = model73.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(3,3)); coda::traceplot(out73.4)     # not shown
(nsum <- nimble_summary(out73.4, params))       # not shown
nimble_est <- nsum[1:3,1]


#   7.3.5 Bayesian analysis with Stan
# -----------------------------------

# Bundle and summarize data (same as before)
str(dataList <- list(y = y, pop = as.numeric(pop), nPops = nPops, n = n))

# Write text file with model description in BUGS language
cat(file = "model73_5.stan",
"data { // Describe input data
  int n;
  int nPops;
  vector[n] y;
  array[n] int pop; // A vector of integers of length n
}

parameters { // Define parameters
  real <lower = 0> pop_grand_mean;
  real <lower = 0> pop_sd;
  real <lower = 0> residual_sd;
  vector <lower = 0> [nPops] pop_mean;
}

model {
  // Vector to hold expected value for each datapoint
  vector[n] mu;

  // Priors
  pop_grand_mean ~ normal(0, 100);
  pop_sd ~ cauchy(0, 10);
  residual_sd ~ cauchy(0, 10);

  for (i in 1:nPops) {
    pop_mean[i] ~ normal(pop_grand_mean, pop_sd);
  }

  // 'Likelihood'
  for(i in 1:n) {
    mu[i] = pop_mean[pop[i]];
    y[i] ~ normal(mu[i], residual_sd);
  }
}

generated quantities {
  vector[nPops] effe;
  for (i in 1:nPops){
    effe[i] = pop_mean[i] - pop_grand_mean;
  }
}
" )

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 45/2 sec), assess convergence, print and save estimates
system.time(
out73.5 <- stan(file = "model73_5.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out73.5) # not shown
print(out73.5, dig = 2) # not shown
stan_est <- summary(out73.5)$summary[1:3,1]


#   7.3.6 Do-it-yourself maximum likelihood estimates
# ---------------------------------------------------

# Definition of NLL for a model with one random-effects factor and Gaussian errors
NLL <- function(pars) {
  pop.grand.mean <- pars[1]
  pop.sd <- exp(pars[2])
  sigma <- exp(pars[3])
  LL <- 0
  for (i in 1:nPops){
    #Subset data to just pop i
    ysub <- y[pop == i]
    L <- integrate(function(r){
      tot <- 1 #Starting value for product over J
      #Iterate over each sample j (1–12) in pop i
      for (j in 1:nSample){
        tot = tot * dnorm(ysub[j], pop.grand.mean + r, sigma)
      }
      tot <- tot * dnorm(r, 0, pop.sd)
      tot
    }, lower = -Inf, upper = Inf, subdivisions = 20)$value
    LL <- LL + log(L)
  }
  return(-LL)
}

# Minimize NLL to find MLEs, also get SEs and CIs and save estimates
inits <- c('pop.grand.mean' = mean(y), 'log(pop.sd)' = 0, 'log(residual.sd)' = log(sd(y)))
out73.6 <- optim(inits, NLL, hessian = TRUE, method = 'BFGS')
get_MLE(out73.6, 4)
exp(out73.6$par[2:3]) # Backtransform the SD's (ignore 'log' in name)
diy_est <- c(out73.6$par[1], exp(out73.6$par[2:3]))


#   7.3.7 Likelihood analysis with TMB
# ------------------------------------

# Bundle and summarize data; adjust pop so it starts at 0 instead of 1
tmbData <- list(y = y, pop = as.numeric(pop)-1, nPops = nPops, n = n)

# Write TMB model file
cat(file = "model73_7.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(y); //response
  DATA_IVECTOR(pop); //Population index: note IVECTOR, not VECTOR
  DATA_INTEGER(nPops); //Number of populations
  DATA_INTEGER(n); //Number of obs
  
  //Describe parameters
  PARAMETER(pop_grand_mean);
  PARAMETER(log_pop_sd)
  PARAMETER(log_residual_sd);
  PARAMETER_VECTOR(pop_mean);
  
  Type pop_sd = exp(log_pop_sd);
  Type residual_sd = exp(log_residual_sd);
  
  Type LL = 0; //Initialize total log likelihood at 0
  
  // Random intercepts
  vector <Type> effe(nPops); //Pop effects
  for (int i = 0; i < nPops; i++){
    LL += dnorm(pop_mean(i), pop_grand_mean, pop_sd, true);
    effe(i) = pop_mean(i) - pop_grand_mean;
  }

  // Likelihood of obs
  for (int i = 0; i < n; i++){ //Note index starts at 0 instead of 1!
    LL += dnorm(y(i), pop_mean(pop(i)), residual_sd, true);
  }
  ADREPORT(effe);
  return -LL; //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model73_7.cpp")
dyn.load(dynlib("model73_7"))

# Provide dimensions and starting values for parameters
params <- list(pop_grand_mean = 0, log_pop_sd = 0, log_residual_sd = 0,
               pop_mean = rep(0, tmbData$nPops))

# Create TMB object
out73.7 <- MakeADFun(data = tmbData,
  parameters = params, random = "pop_mean",
  DLL = "model73_7", silent = TRUE)

# Optimize TMB object and print results (including random effects)
opt <- optim(out73.7$par, fn = out73.7$fn, gr = out73.7$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out73.7)) # not shown

# Save estimates from TMB
tmb_est <- c(opt$par[1], exp(opt$par[2:3]))


#   7.3.8 Comparison of the parameter estimates
# ---------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, ML = gtmb_est, REML = gtmb_REMLest, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 7.4 Summary and outlook
# -----------------------
# (no code)
