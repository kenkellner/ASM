
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -----------------------------------------------------
# Chapter 8  --  Comparisons along two classifications
#                in a model with two factors
# -----------------------------------------------------

# Last changes: 11 June 2024


# 8.1 Introduction: main and interaction effects
# ----------------------------------------------
# (no code)


# 8.2 Data generation
# -------------------

set.seed(8)

# Choose sample size
nPops <- 5
nHab <- 3
nSample <- 4
n <- nPops * nHab * nSample

# Create factor levels
pop <- gl(n = nPops, k = nHab * nSample, length = n)
hab <- gl(n = nHab, k = nSample, length = n)

# Choose effects
baseline <- 40                       # Intercept
pop.eff <- c(-10, -5, 5, 10)         # Population effects
hab.eff <- c(5, 10)                  # Hab effects
interaction.eff <- c(-2, 3, 0, 4, 4, 0, 3, -2) # Interaction effects
all.eff <- c(baseline, pop.eff, hab.eff, interaction.eff)
sigma <- 3                           # Residual standard deviation

eps <- rnorm(n, 0, sigma)            # Residuals
Xmat <- as.matrix(model.matrix(~ pop * hab) ) # Create design matrix
Xmat     # Have a look at design matrix, make sure you understand it!
str(Xmat)                            # not shown

wing <- as.vector(Xmat %*% all.eff +  eps)
boxplot(wing ~ hab*pop, col = "grey", xlab = "Habitat-by-Population combination",
  ylab = "Wing length", main = "Simulated data set", las = 1,
  ylim = c(20, 70), frame = FALSE)   # Plot of generated data
abline(h = 40)

library(lattice)     # load the lattice library
xyplot(wing ~ hab | pop, ylab = "Wing length", xlab = "Habitat",
  main = "Population-specific relationship between wing and habitat type",
  pch = 16, cex = 2, col = rgb(0,0,0,0.4))
xyplot(wing ~ pop | hab, ylab = "Wing length", xlab = "Population", 
  main = "Habitat-specific relationship between wing and population",    
  pch = 16, cex = 2, col = rgb(0,0,0,0.6)) # Fig. 8.3

# Save true values for later comparisons
truth <- c(all.eff, sigma)
( true.group.means <- c(baseline, baseline + pop.eff, baseline + hab.eff[1],
  baseline + pop.eff + hab.eff[1] + interaction.eff[1:4], baseline + hab.eff[2],
  baseline + pop.eff + hab.eff[2] + interaction.eff[5:8]) )

#Required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB); library(DHARMa)


# 8.3 Likelihood analysis with canned functions in R
# --------------------------------------------------

# Fit main-effects model and save estimates
summary( out83.main <- lm(wing ~ hab + pop) )
lm_est_me <- c(coef(out83.main), sigma = sigma(out83.main))

# Interaction-effects model (means parameterization)
summary(out83.intX <- lm(wing ~ pop:hab-1) ) # not shown

# Interaction-effects model (treatment contrast parameterization)
summary(out83.int <- lm(wing ~ pop*hab) )

# GOF
# Check goodness-of-fit using quantile residuals
# Main-effects model
plot(out83.main) # Traditional residual check for comparison
simOut <- simulateResiduals(out83.main, n = 1000, plot = TRUE) # DHARMa
# Interaction-effects model
plot(out83.int) # Traditional residual check for comparison
simOut <- simulateResiduals(out83.int, n = 1000, plot = TRUE) # DHARMa


# 8.4 An aside: using simulation to assess bias and precision 
#     of an estimator ... and to understand what a standard error is 
# ------------------------------------------------------------------

# Compare estimates with truth
lm_est <- c(coef(out83.int), sigma = sigma(out83.int))
comp <- cbind(truth = truth, lm = lm_est)
print(comp, 4)

simrep <- 10^4 # Desired number of iterations
estimates <- array(dim = c(simrep, length(all.eff))) # Data structure to hold results
for(i in 1:simrep) {                 # Run simulation simrep times
  if(i %% 100 == 0) cat(paste('\n iter', i))  # Counter
  eps <- rnorm(n, 0, sigma)          # Residuals
  y <- as.vector(Xmat %*% all.eff + eps) # Assemble data
  fit.model <- lm(y ~ pop*hab)       # Break down data
  estimates[i,] <- fit.model$coef    # Keep values of coefficients
}

# Compare true values and mean of estimates
data.frame(params = names(fit.model$coefficients), truth = all.eff,
  mean.of.estimates = round(apply(estimates, 2, mean), 3))

# Compare theoretical and bootstrapped standard errors
SE.lm <- summary(out83.int)$coef[,2]
  data.frame(SE.lm = round(SE.lm, 3), bootstrapped.SE = round(apply(estimates, 2, sd), 3))

# Plot bootstrapped sampling distributions of four estimators (Fig. 8.4)
par(mfrow = c(2,2), mar = c(5,6,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(estimates[,1], xlab = '', freq = F, main = 'Intercept', col = 'grey', breaks = 40)
hist(estimates[,2], xlab = '', freq = F, main = 'pop2', col = 'grey', breaks = 40)
hist(estimates[,6], xlab = '', freq = F, main = 'hab2', col = 'grey', breaks = 40)
hist(estimates[,8], xlab = '', freq = F, main = ' pop2:hab2', col = 'grey', breaks = 40)


# 8.5 Bayesian analysis with JAGS
# -------------------------------
# (no code)


#   8.5.1 Main-effects model
# --------------------------

# Bundle and summarize data
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing)))

# Write JAGS model file
cat(file = "model85.main.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)             # Intercept
beta.pop[1] <- 0                     # set to zero effect of 1st level
for(j in 2:5){                       # Loop over remaining factor levels
  beta.pop[j] ~ dnorm(0, 0.0001)
}
beta.hab[1] <- 0                     # ditto
for(j in 2:3){                       # Loop over remaining factor levels
  beta.hab[j] ~ dnorm(0, 0.0001)
}
sigma ~ dunif(0, 100)
tau <- pow(sigma, -2)

# Likelihood
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau)
  mean[i] <- alpha + beta.pop[pop[i]] + beta.hab[hab[i]]
}
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(1), sigma = rlnorm(1) )}

# Parameters to estimate
params <- c("alpha", "beta.pop", "beta.hab", "sigma")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out85me <- jags(dataList, inits, params, "model85.main.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out85me) # not shown
print(out85me, 3)
jags_est_me <- out85me$summary[c(1, 8:9, 3:6, 10), 1]

# Compare likelihood and Bayesian estimates
comp <- cbind(lm = lm_est_me, JAGS = jags_est_me)
print(comp, 4)


#   8.5.2 Interaction-effects model
# ---------------------------------

# Bundle and summarize data
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = length(unique(hab)), nPops = length(unique(pop))) )

# Write JAGS model file
cat(file = "model85.int.txt", "
model {
# Priors
for (i in 1:nPops){
  for(j in 1:nHab) {
    group.mean[i,j] ~ dnorm(0, 0.0001)
  }
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# Likelihood
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau)
  mean[i] <- group.mean[pop[i], hab[i]]
}
}
")

# Function to generate starting values
inits <- function(){list(sigma = rlnorm(1) )}

# Parameters to estimate
params <- c("group.mean", "sigma")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out85ie <- jags(dataList, inits, params, "model85.int.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out85ie) # not shown
print(out85ie, 3)
jags_est_ie <- out85ie$summary[1:16,1]

# Align true vals of all params
truthMP <- c(true.group.means, sigma) # Truth for means parameterization

# Compare likelihood with Bayesian estimates and with truth
lm_est_ie <- c(coef(out83.intX), sigma = sigma(out83.intX))
comp <- cbind(truth = truthMP, lm = lm_est_ie, JAGS = jags_est_ie)
print(comp, 4)


#   8.5.3 Forming predictions 
# ---------------------------
par(mfrow = c(1,1), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
ord <- c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15) # select order
plot(ord, out85ie$mean$group.mean, xlab = "Hab-by-Population", las = 1,
  ylab = "Predicted wing length", cex = 2, ylim = c(20, 70), frame = FALSE,
  pch = 16, col = rgb(0,0,0,0.5))
segments(ord, out85ie$q2.5$group.mean, ord, out85ie$q97.5$group.mean,
  col = rgb(0,0,0,0.5), lwd = 2)
abline(h = 40)


# 8.6 Bayesian analysis with NIMBLE
# ---------------------------------
# (no code)


#   8.6.1 Main-effects model
# --------------------------

library(nimble)

# Bundle and summarize data (same as for JAGS)
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing)))

# Write NIMBLE model file
model86.main <- nimbleCode( {
# Priors
alpha ~ dnorm(0, 0.0001)             # Intercept
beta.pop[1] <- 0                     # set to zero effect of 1st level
for(j in 2:5){                       # Loop over remaining factor levels
  beta.pop[j] ~ dnorm(0, 0.0001)
}
beta.hab[1] <- 0                     # ditto
for(j in 2:3){                       # Loop over remaining factor levels
  beta.hab[j] ~ dnorm(0, 0.0001)
}
sigma ~ dunif(0, 100)
tau <- pow(sigma, -2)

# Likelihood
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau) 
  mean[i] <- alpha + beta.pop[pop[i]] + beta.hab[hab[i]]
}
} )

# Inits
inits <- function(){ list(alpha = rnorm(1), sigma = rlnorm(1) )}

# Parameters monitored: same as before
params <- c("alpha", "beta.pop", "beta.hab", "sigma")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 30 sec), check convergence, summarize posteriors and save estimates
system.time( out86me <- 
    nimbleMCMC(code = model86.main,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out86me)   # not shown
(nsum <- nimble_summary(out86me, params))     # not shown
nimble_est_me <- nsum[c(1, 8:9, 3:6, 10), 1]


#   8.6.2 Interaction-effects model 
# ---------------------------------

# Bundle and summarize data (same as for JAGS)
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = length(unique(hab)), nPops = length(unique(pop))) )

# Write NIMBLE model file
model86.int <- nimbleCode( {
# Priors
for (i in 1:nPops){
  for(j in 1:nHab) {
    group.mean[i,j] ~ dnorm(0, 0.0001)
  }
}
tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# Likelihood
for (i in 1:n) {
  wing[i] ~ dnorm(mean[i], tau) 
  mean[i] <- group.mean[pop[i], hab[i]]
}
} )

# Inits
inits <- function(){list(sigma = rlnorm(1) )}

# Parameters monitored: same as before
params <- c("group.mean", "sigma")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 20 sec), check convergence, summarize posteriors and save results
system.time( out86ie <- 
    nimbleMCMC(code = model86.int,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out86ie)   # not shown
(nsum <- nimble_summary(out86ie, params))     # not shown
nimble_est_ie <- nsum[1:16, 1]


# 8.7 Bayesian analysis with Stan 
# -------------------------------
# (no code)


#   8.7.1 Main-effects model
# --------------------------
# Bundle and summarize data (same as for JAGS)
dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop), n = length(wing))


# Write text file with model description in Stan
cat(file = "model87_main.stan",      # This line is R code
"
data {
  int n;
  vector[n] wing;
  int pop[n];
  int hab[n];
}

parameters {
  real alpha;
  vector[5] beta_pop_raw;
  vector[3] beta_hab_raw;
  real<lower=0> sigma;
}

transformed parameters {
  vector[5] beta_pop = beta_pop_raw;
  vector[3] beta_hab = beta_hab_raw;
  beta_pop[1] = 0;
  beta_hab[1] = 0;
}

model {
  vector [n] mu;
  beta_pop_raw ~ normal(0, 100);
  beta_hab_raw ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  for (i in 1:n){
    mu[i] = alpha + beta_pop[pop[i]] + beta_hab[hab[i]];
    wing[i] ~ normal(mu[i], sigma);
  }
}                                    // This is the last line of Stan code
" )

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Parameters to save
params <- c("alpha", "beta_pop", "beta_hab", "sigma")

# Call STAN (ART 50/3 sec)
system.time(
  out87me <- stan(file = "model87_main.stan", data = dataList,
    pars = params, chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out87me)            # not shown
print(out87me, dig = 2)              # not shown
stan_est_me <- summary(out87me)$summary[c(1, 8, 9, 3, 4, 5, 6, 10), 1]


#   8.7.2 Interaction-effects model
# ---------------------------------

# Bundle and summarize data
str(dataList <- list(wing = wing, hab = as.numeric(hab), pop = as.numeric(pop),
  n = length(wing), nHab = nHab, nPops = nPops))

# Write text file with model description in Stan
cat(file = "model87_int.stan",       # This line is R code
"
data {
  int n;
  int nPops;
  int nHab;
  vector[n] wing; int hab[n];
  int pop[n];
}

parameters {
  matrix[nPops, nHab] group_mean;
  real<lower=0> sigma;
}

model {
  to_vector(group_mean) ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  for (i in 1:n) {
    wing[i] ~ normal(group_mean[pop[i], hab[i]], sigma);
  }
}                                    // This is the last line of Stan code
" )

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 48/3 sec)
system.time(
  out87ie <- stan(file = "model87_int.stan", data = dataList, chains = nc, iter = ni,
  warmup = nb, thin = nt) )
rstan::traceplot(out87ie)            # not shown
print(out87ie, dig = 2)              # not shown
# Save estimates
# Stan sorts summary by row, then by column, opposite of JAGS
stan_est_ie <- summary(out87ie)$summary[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15,16),1]


# 8.8 Do-it-yourself MLEs
# -----------------------
# (no code)


#   8.8.1 Main-effects model
# --------------------------

# Design matrix for main-effects model
head( Xmat <- model.matrix(~ hab + pop))

# NLL for main-effects normal linear model with two factors
NLL <- function(param, y, Xmat) {
  beta <- param[1:7]                 # Intercept and slope coefficients
  sigma <- exp(param[8])             # Residual SD
  mu <- Xmat %*% beta                # Multiply each row of Xmat by beta
  LL <- dnorm(y, mu, sigma, log = TRUE) # Log-likelihood for each obs
  NLL <- -sum(LL)                    # NLL for all observations (whole data set)
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c(mean(wing), rep(0, 6), log(sd(wing)))
names(inits) <- c(names(coef(out83.main)), 'log-sigma')
out88me <- optim(inits, NLL, y = wing, Xmat = Xmat, hessian = TRUE, method = 'BFGS')
get_MLE(out88me, 4)
diy_est_me <- c(out88me$par[-8], exp(out88me$par[8]))


#   8.8.2 Interaction-effects model
# ---------------------------------

# Design matrix for interaction-effects model
head( Xmat <- model.matrix(~ pop:hab-1)) # not shown

# NLL for interaction-effects normal linear model with two factors
NLL <- function(param, y, Xmat) {
  beta <- param[1:15]                # Interaction parameters
  sigma <- exp(param[16])            # Residual SD
  mu <- Xmat %*% beta
  LL <- dnorm(y, mu, sigma, log = TRUE) # Log-lik for each observation
  NLL <- -sum(LL)                    # NLL for all observations (whole data set)
  return(NLL)
}
# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c(rep(mean(wing), 15), log(sd(wing)))
names(inits) <- c(colnames(Xmat), 'log-sigma')
out88.int <- optim(inits, NLL, y = wing, Xmat = Xmat, hessian = TRUE, method = 'BFGS')
get_MLE(out88.int, 4)
diy_est_ie <- c(out88.int$par[-16], exp(out88.int$par[16]))


# 8.9 Likelihood analysis with TMB
# --------------------------------
# (no code)


#   8.9.1 Main-effects model 
# --------------------------
# Bundle and summarize data
Xmat <- model.matrix(~ hab + pop)
str(tmbData <- list(wing = wing, Xmat = Xmat, n = length(wing)))

# Write TMB model file
cat(file = "model89_main.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(wing);                 //response
  DATA_MATRIX(Xmat);                 //Model matrix
  DATA_INTEGER(n);                   //Number of obs

  //Describe parameters
  PARAMETER_VECTOR(beta);            //Intercept and slopes
  PARAMETER(log_sigma);              //log(residual standard deviation)

  Type sigma = exp(log_sigma);
  Type LL = 0;                       //Initialize total log likelihood at 0
  
  matrix<Type> mu = Xmat * beta;     //Matrix multiply Xmat and beta
  
  for (int i=0; i< n; i++){          //Note index starts at 0 instead of 1
    LL += dnorm(wing(i), mu(i), sigma, true);
  }
  return -LL;                        //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model89_main.cpp")
dyn.load(dynlib("model89_main"))

# Provide dimensions and starting values for parameters
params <- list(beta = rep(0, ncol(Xmat)), log_sigma = 0)

# Create TMB object
out89me <- MakeADFun(data = tmbData, parameters = params,
DLL = "model89_main", silent = TRUE)

# Optimize TMB object, print results and save estimates
opt1 <- optim(out89me$par, fn = out89me$fn, gr = out89me$gr, method="BFGS", hessian = TRUE)
(tsum <- tmb_summary(out89me))
tmb_est <- c(opt1$par[1:7], exp(opt1$par[8]))


#   8.9.2 Interaction-effects model
# ---------------------------------

# Bundle and summarize data
str(tmbData <- list(wing = wing,
pop = as.numeric(pop)-1,             # TMB indices start at 0
hab = as.numeric(hab)-1,             # TMB indices start at 0
n= length(wing)))

# Write TMB model file
cat(file = "model89_int.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(wing);                 //response
  DATA_IVECTOR(pop);                 //Population of each obs
  DATA_IVECTOR(hab);                 //Habitat group of each obs
  DATA_INTEGER(n);                   //Number of obs
  
  //Describe parameters
  PARAMETER_MATRIX(group_mean);
  PARAMETER(log_sigma);              //log(residual standard deviation)

  Type sigma = exp(log_sigma);       //Type = match type of function output
  Type LL = 0.0;                     //Initialize total log likelihood at 0

  for (int i= 0; i<n; i++){          //Note index starts at 0 instead of 1!
    LL += dnorm(wing(i), group_mean(pop(i), hab(i)), sigma, true);
  }
  return -LL;                        //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model89_int.cpp")
dyn.load(dynlib("model89_int"))

# Provide dimensions and starting values for parameters
params <- list(group_mean = matrix(0, nPops, nHab), log_sigma=0)

# Create TMB object
out89ie <- MakeADFun(data = tmbData, parameters = params,
  DLL = "model89_int", silent = TRUE)

# Optimize TMB object, print results and save estimates
starts <- rep(0, nPops*nHab+ 1)
opt2 <- optim(starts, fn = out89ie$fn, gr = out89ie$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out89ie)) # not shown
tmb_est_ie <- c(opt2$par[1:15], exp(opt2$par[16]))


# 8.10 Comparison of the parameter estimates 
# ------------------------------------------

# Compare interaction-effects model estimates with truth
comp <- cbind(truth = truthMP, lm = lm_est_ie, JAGS = jags_est_ie, NIMBLE = nimble_est_ie,
  Stan = stan_est_ie, DIY = diy_est_ie, TMB = tmb_est_ie)
print(comp, 4)


# 8.11 Summary and outlook 
# ------------------------
# (no code)
