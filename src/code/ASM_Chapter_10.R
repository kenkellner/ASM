
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# --------------------------------------------
# Chapter 10  --  Linear mixed-effects model
# --------------------------------------------

# Last changes: 11 June 2024


# 10.1 Introduction
# -----------------
# (no code)


# 10.2 Data generation
# --------------------

set.seed(10)
nPops <- 56 # Number of populations
nSample <- 10 # Number of vipers in each pop
n <- nPops * nSample # Total number of data points
pop <- gl(n = nPops, k = nSample) # Indicator for population

# Body length (cm)
orig.length <- runif(n, 45, 70)
mn <- mean(orig.length)
sd <- sd(orig.length)
cat("Mean and sd used to normalise original length:", mn, sd, "\n\n")
lengthN <- (orig.length - mn)/sd # N for 'normalized'
hist(lengthN, col = "grey") # Not shown

Xmat <- model.matrix(~ pop * lengthN - 1 - lengthN)
print(Xmat[1:21,], dig = 2) # Print top 21 rows (not shown)

mu.alpha <- 260
sigma.alpha <- 20
mu.beta <- 60
sigma.beta <- 30

alpha <- rnorm(n = nPops, mean = mu.alpha, sd = sigma.alpha)
beta <- rnorm(n = nPops, mean = mu.beta, sd = sigma.beta)
all.pars <- c(alpha, beta) # Put them all together

sigma <- 30 # Residual standard deviation
lin.pred <- Xmat[,] %*% all.pars # Value of lin.predictor
eps <- rnorm(n = n, mean = 0, sd = sigma) # residuals
mass <- lin.pred + eps # response = lin.pred + residual
mass <- drop(mass) # Turn 1xn matrix into vector of length n
hist(mass, col = "grey") # Inspect what we’ve created (not shown)

# Save true values for comparisons later
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta,
sigma.alpha = sigma.alpha, sigma.beta = sigma.beta, residual.sd = sigma)

library(lattice)
xyplot(mass ~ lengthN | pop, xlab = 'Length', ylab = 'Mass', main = 'Realized mass-length relationships',
  pch = 16, cex = 1.2, col = rgb(0, 0, 0, 0.4))

library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 10.3 Analysis under a random-intercepts model
# ---------------------------------------------
# (no code)


#   10.3.1 Likelihood analysis using canned functions in R
# --------------------------------------------------------

library(glmmTMB)
gtmbData <- data.frame(mass = mass, lengthN = lengthN, pop = pop)
out10.3.ML <- glmmTMB(mass ~ lengthN + (1 | pop), data = gtmbData, REML = FALSE,
                      control = glmmTMBControl(optimizer=optim, 
                                               optArgs=list(method="L-BFGS-B")))  # ML
out10.3.REML <- glmmTMB(mass ~ lengthN + (1 | pop), data = gtmbData, REML = TRUE) # REML
summary(out10.3.ML)
summary(out10.3.REML)

# Save estimates from ML and from REML fit of random-intercepts model
gtmbML_est <- c(fixef(out10.3.ML)$cond,
                attr(VarCorr(out10.3.ML)$cond$pop, "stddev"),
                sigma(out10.3.ML))
gtmbREML_est <- c(fixef(out10.3.REML)$cond,
                  attr(VarCorr(out10.3.REML)$cond$pop, "stddev"),
                  sigma(out10.3.REML))

# Assess goodness-of-fit of the model (not shown)
hist(coef(out10.3.ML)$cond$pop[,1], main = 'alpha', breaks = 12)
library(DHARMa)
simOut <- simulateResiduals(out10.3.ML, n = 1000, plot = TRUE)


#   10.3.2 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(mass = mass, pop = as.numeric(pop), lengthN = lengthN, nPops = nPops, n = n) )

# Write JAGS model file
cat(file = "model10.3.2.txt", "
model {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # Random intercepts
}

mu.alpha ~ dnorm(0, 0.000001) # Mean for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 1, 0.0001)I(0,) # SD for random intercepts

mu.beta ~ dnorm(0, 0.000001) # Common slope
tau.residual <- pow(sigma.residual, -2) # Residual precision
sigma.residual ~ dt(0, 1, 0.0001)I(0,) # Residual standard deviation

# Likelihood
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual) # The observed random variables
  mu[i] <- alpha[pop[i]] + mu.beta * lengthN[i] # Expectation
}

# Simulation of new data not conditional on random effects for GoF assessment using DHARM
for (i in 1:nPops){
  alpha.new[i] ~ dnorm(mu.alpha, tau.alpha)
}
for (i in 1:n) {
  mass.new.cond[i] ~ dnorm(mu[i], tau.residual) # Conditional
  mass.new.un[i] ~ dnorm(mu.new[i], tau.residual) # Unconditional
  mu.new[i] <- alpha.new[pop[i]] + mu.beta * lengthN[i]
}
}
")

# Function to generate starting values
inits <- function(){list(mu.alpha = rnorm(1, 0, 1),
  mu.beta = rnorm(1, 0, 1), sigma.alpha = rlnorm(1),
  sigma.residual = rlnorm(1)) }

# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.residual", "alpha",
  "mass.new.cond", "mass.new.un")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out10.3.2 <- jags(dataList, inits, params, "model10.3.2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
par(mfrow = c(2, 2)); jagsUI::traceplot(out10.3.2) # not shown
print(out10.3.2, 2)
jags_est <- out10.3.2$summary[c(1:4),1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth[-4], gtmbML = gtmbML_est, JAGS = jags_est)
print(comp, 4)

# Distribution of random effects estimates (alpha)
hist(out10.3.2$mean$alpha, main = 'alpha')

# Quantile residual assessments (not shown)
mass.new.un <- out10.3.2$sims.list$mass.new.un # Unconditional
mass.new.cond <- out10.3.2$sims.list$mass.new.cond # Conditional

sim.un <- createDHARMa(simulatedResponse = t(mass.new.un), observedResponse = mass,
  fittedPredictedResponse = apply(mass.new.un, 2, median), integerResponse = FALSE)
plot(sim.un) # Unconditional new data
sim.cond <- createDHARMa(simulatedResponse = t(mass.new.cond), observedResponse = mass,
  fittedPredictedResponse = apply(mass.new.cond, 2, median), integerResponse = FALSE)
plot(sim.cond) # New data conditional on estimated random effects


#   10.3.3 Bayesian analysis with NIMBLE
# --------------------------------------

library(nimble)

# Summarize data again
str(dataList)

# Write NIMBLE model file
model10.3.3 <- nimbleCode( {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha)   # Random intercepts
}

mu.alpha ~ dnorm(0, 0.000001)  # Mean for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ T(dt(0, 1, 0.0001), 0, )    # SD for random intercepts

mu.beta ~ dnorm(0, 0.000001)   # Common slope
tau.residual <- pow(sigma.residual, -2)   # Residual precision
sigma.residual ~ T(dt(0, 1, 0.0001), 0, ) # Residual standard deviation

# Likelihood
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual) # The observed random variables
  mu[i] <- alpha[pop[i]] + mu.beta * lengthN[i]  # Expectation
}
} )

# Inits
inits <- function(){list(mu.alpha = rnorm(1, 0, 1), 
  mu.beta = rnorm(1, 0, 1), sigma.alpha = rlnorm(1), 
  sigma.residual = rlnorm(1)) }

# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.residual", "alpha")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 32 sec), check convergence, summarize posteriors and save point estimates
system.time( out10.3.3 <-
    nimbleMCMC(code = model10.3.3,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out10.3.3)
(nsum <- nimble_summary(out10.3.3, params))   # not shown
nimble_est <- nsum[1:4,1]


#   10.3.4 Bayesian analysis with Stan
# ------------------------------------

# Summarize data set again
str(dataList) # not shown

# Write Stan model
cat(file = "model10_3_4.stan", "
data {
  int n; //Number of observations
  int nPops; //Number of populations
  vector[n] mass; //Response variable
  vector[n] lengthN; //Covariate
  array[n] int pop; //Population assignment of each obs
}

parameters {
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_residual;
  vector[nPops] alpha;
}

model {
  vector[n] mu; //Expected value

  mu_alpha ~ normal(0, 1000);
  mu_beta ~ normal(0, 1000);
  sigma_alpha ~ cauchy(0, 100);
  sigma_residual ~ cauchy(0, 100);

  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
  }

  //'Likelihood'
  for (i in 1:n){
    mu[i] = alpha[pop[i]] + mu_beta * lengthN[i];
    mass[i] ~ normal(mu[i], sigma_residual);
  }
}
")

# Function to generate starting values
inits <- function(){list(mu_alpha = rnorm(1, 0, 1),
  mu_beta = rnorm(1, 0, 1), sigma_alpha = rlnorm(1),
  sigma_residual = rlnorm(1)) }

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 50 sec/22 sec) and save results
system.time(
out10.3.4 <- stan(file = "model10_3_4.stan", data = dataList, init = inits,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out10.3.4) # not shown
print(out10.3.4, dig = 3) # not shown
stan_est <- summary(out10.3.4)$summary[1:4,1]


#   10.3.5 Do-it-yourself MLEs
# ----------------------------

f <- function(alpha.j, mu.alpha, sigma.alpha, mu.beta, sigma.residual,
  mass.j, lengthN.j){
  mu <- alpha.j + mu.beta * lengthN.j
  prod(dnorm(mass.j, mu, sigma.residual)) *
  dnorm(alpha.j, mu.alpha, sigma.alpha)
}

f <- function(alpha.j, mu.alpha, sigma.alpha, mu.beta, sigma.residual,
              mass.j, lengthN.j){
  out <- numeric(length(alpha.j))
  for (i in 1:length(alpha.j)){
    mu <- alpha.j[i] + mu.beta * lengthN.j
    out[i] <- prod(dnorm(mass.j, mu, sigma.residual)) *
    dnorm(alpha.j[i], mu.alpha, sigma.alpha)
  }
  out
}

g <- function(mu.alpha, mu.beta, sigma.alpha, sigma.residual,
              mass.j, lengthN.j){
  integrate(f, lower = 100, upper = 400, mu.alpha = mu.alpha,
            sigma.alpha = sigma.alpha,
            mu.beta = mu.beta, sigma.residual = sigma.residual,
            mass.j = mass.j, lengthN.j = lengthN.j)$value
}

# Definition of NLL for random-intercepts model with Gaussian errors
NLL <- function(pars, data){
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sigma.alpha <- exp(pars[3])
  sigma.residual <- exp(pars[4])
  LL <- 0 # initialize nll at 0
  
  for (j in 1:data$nPops){
    mass.j <- data$mass[data$pop == j]
    lengthN.j <- data$lengthN[data$pop == j]
    L <- g(mu.alpha, mu.beta, sigma.alpha, sigma.residual,
           mass.j, lengthN.j) # likelihood for pop j
    LL <- LL + log(L)
  }
  
  return(-LL)
}

# Minimize that NLL to find MLEs and get SEs and CIs
inits <- c(mu.alpha = mean(dataList$mass), mu.beta = 0,
           log.sigma.alpha = 1,
           log.sigma.residual = log(sd(dataList$mass)))
out10.3.5 <- optim(inits, NLL, hessian = TRUE, data = dataList)
get_MLE(out10.3.5, 4)

# Save MLEs
diy_est <- out10.3.5$par
diy_est[3:4] <- exp(diy_est[3:4]) # Get SD params on ilog scale


#   10.3.6 Likelihood analysis with TMB
# -------------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop-1 # convert pop index to 0-based
str(tmbData)

# Write TMB model file
cat(file = "model10_3_6.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(mass); //response
  DATA_VECTOR(lengthN); //covariate
  DATA_IVECTOR(pop); //population index
  DATA_INTEGER(nPops); //Number of populations
  DATA_INTEGER(n); //Number of observations
  
  //Describe parameters
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER(log_sigma_alpha);
  PARAMETER(log_sigma_residual);
  PARAMETER_VECTOR(alpha);

  Type sigma_alpha = exp(log_sigma_alpha);
  Type sigma_residual = exp(log_sigma_residual);
  
  Type LL = 0.0; //Initialize log-likelihood at 0

  //Distribution of random intercepts – add to log-likelihood
  for (int i = 0; i < nPops; i++){ //Note index starts at 0 instead of 1!
    LL += dnorm(alpha(i), mu_alpha, sigma_alpha, true);
  }
  
  for (int i = 0; i < n; i++){
    Type mu = alpha(pop(i)) + mu_beta * lengthN(i);
    //Calculate log-likelihood of observation and add to total
    LL += dnorm(mass(i), mu, sigma_residual, true);
  }
  
  return -LL; //Return negative log likelihood
}
")

# Compile and load TMB function
compile("model10_3_6.cpp")
dyn.load(dynlib("model10_3_6"))

# Provide dimensions and starting values for parameters
params <- list(mu_alpha = 0, mu_beta = 0,
               log_sigma_alpha = 0, log_sigma_residual = 0,
               alpha = rep(0, tmbData$nPops))

# Create TMB object and tell TMB that alpha is random
out10.3.6 <- MakeADFun(data = tmbData, parameters = params,
              random = "alpha", #Identify which params(s) are random
              DLL = "model10_3_6", silent = TRUE)

# Optimize TMB object and print results
opt <- optim(out10.3.6$par, fn = out10.3.6$fn, gr = out10.3.6$gr,
             method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out10.3.6)) # not shown
# Save TMB estimates
tmb_est <- c(opt$par[1:2], exp(opt$par[3:4]))


#   10.3.7 Comparison of the parameter estimates
# ----------------------------------------------

# Compare all results with truth
comp <- cbind(truth = truth[-4], gtmbML = gtmbML_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 10.4 Analysis under a random-coefficients model without correlation between
# intercept and slope
# ---------------------------------------------------------------------------
# (no code)


#   10.4.1 Likelihood analysis using canned functions in R
# --------------------------------------------------------

# the || notation means no correlation between random effects
out10.4.ML <- glmmTMB(mass ~ lengthN + ( 1 + lengthN || pop),
                      data = gtmbData, REML = FALSE)
out10.4.REML <- glmmTMB(mass ~ lengthN + ( 1 + lengthN || pop),
                        data = gtmbData, REML = TRUE,
                        control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(out10.4.ML)
summary(out10.4.REML)

# Save estimates
gtmbML_est <- c(fixef(out10.4.ML)$cond,
                attr(VarCorr(out10.4.ML)$cond$pop, "stddev"),
                sigma(out10.4.ML))
gtmbREML_est <- c(fixef(out10.4.REML)$cond,
                  attr(VarCorr(out10.4.REML)$cond$pop, "stddev"),
                  sigma(out10.4.REML))

# Assess goodness-of-fit of the model (not shown)
par(mfrow = c(1, 2)) # Distribution of estimates of alpha and beta
hist(coef(out10.4.ML)$cond$pop[,1], main = 'alpha', breaks = 12)
hist(coef(out10.4.ML)$cond$pop[,2], main = 'beta', breaks = 12)
simOut <- simulateResiduals(out10.4.ML, n = 1000, plot = TRUE)


#   10.4.2 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(mass = mass, pop = as.numeric(pop), lengthN = lengthN,
  nPops = nPops, n = n) )

# Write JAGS model file
cat(file = "model10.4.2.txt", "
model {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha) # Random intercepts
  beta[i] ~ dnorm(mu.beta, tau.beta) # Random slopes
}

mu.alpha ~ dnorm(0, 0.000001) # Mean for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ dt(0, 1, 0.0001)I(0,) # SD for random intercepts

mu.beta ~ dnorm(0, 0.000001) # Mean for random slopes
tau.beta <- pow(sigma.beta, -2)
sigma.beta ~ dt(0, 1, 0.0001)I(0,) # SD for slopes

tau.residual <- pow(sigma.residual, -2) # Residual precision
sigma.residual ~ dt(0, 1, 0.0001)I(0,) # Residual standard deviation

# Likelihood
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthN[i]
}}
")

# Function to generate starting values
inits <- function(){ list(mu.alpha = rnorm(1, 0, 1),
  sigma.alpha = rlnorm(1), mu.beta = rnorm(1, 0, 1),
  sigma.beta = rlnorm(1), sigma.residual = rlnorm(1)) }

# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
            "sigma.residual", "alpha", "beta")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out10.4.2 <- jags(dataList, inits, params, "model10.4.2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out10.4.2) # not shown
print(out10.4.2, 2)
jags_est <- out10.4.2$summary[c(1:5),1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, gtmbML = gtmbML_est, JAGS = jags_est)
print(comp, 4)

coef(out10.4.ML)


#   10.4.3 Bayesian analysis with NIMBLE
# --------------------------------------

# Summarize data again
str(dataList)                     # not shown

# Write NIMBLE model file
model10.4.3 <- nimbleCode( {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(mu.alpha, tau.alpha)  # Random intercepts
  beta[i] ~ dnorm(mu.beta, tau.beta)     # Random slopes
}

mu.alpha ~ dnorm(0, 0.000001)  # Mean for random intercepts
tau.alpha <- pow(sigma.alpha, -2)
sigma.alpha ~ T(dt(0, 1, 0.0001), 0, ) # SD for random intercepts

mu.beta ~ dnorm(0, 0.000001)      # Mean for random slopes
tau.beta <- pow(sigma.beta, -2)
sigma.beta ~ T(dt(0, 1, 0.0001), 0, ) # SD for slopes

tau.residual <- pow(sigma.residual, -2)     # Residual precision
sigma.residual ~ T(dt(0, 1, 0.0001), 0, )      # Residual standard deviation

# 'Likelihood'
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau.residual)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthN[i]
}
} )

# Inits
inits <- function(){ list(mu.alpha = rnorm(1, 0, 1), 
  sigma.alpha = rlnorm(1), mu.beta = rnorm(1, 0, 1), 
  sigma.beta = rlnorm(1), sigma.residual = rlnorm(1)) }

# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
            "sigma.residual", "alpha", "beta")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 46 sec), check convergence, summarize posteriors and save estimates
system.time( out10.4.3 <- 
    nimbleMCMC(code = model10.4.3,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out10.4.3)
(nsum <- nimble_summary(out10.4.3, params))
nimble_est <- nsum[1:5,1]


#   10.4.4 Bayesian analysis with Stan
# ------------------------------------

# Summarize data set again
str(dataList) # not shown

# Write Stan model
cat(file = "model10_4_4.stan", "
data {
  int n; //Number of observations
  int nPops; //Number of populations
  vector[n] mass; //Response variable
  vector[n] lengthN; //Covariate
  array[n] int pop; //Population assignment of each obs
}

parameters {
  real mu_alpha;
  real mu_beta;
  real <lower = 0> sigma_alpha;
  real <lower = 0> sigma_beta;
  real <lower = 0> sigma_residual;
  vector[nPops] alpha;
  vector[nPops] beta;
}

model {
  vector[n] mu; //Expected value

  //Priors
  mu_alpha ~ normal(0, 1000);
  mu_beta ~ normal(0, 1000);
  sigma_alpha ~ cauchy(0, 100);
  sigma_beta ~ cauchy(0, 100);
  
  for (i in 1:nPops){
    alpha[i] ~ normal(mu_alpha, sigma_alpha);
    beta[i] ~ normal(mu_beta, sigma_beta);
  }

  //Likelihood
  for (i in 1:n){
    mu[i] = alpha[pop[i]] + beta[pop[i]] * lengthN[i];
    mass[i] ~ normal(mu[i], sigma_residual);
  }
}
")

# Initial values
inits <- function(){ list(mu_alpha = rnorm(1, 0, 1),
  sigma_alpha = rlnorm(1), mu_beta = rnorm(1, 0, 1),
  sigma_beta = rlnorm(1), sigma_residual = rlnorm(1)) }

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 50 sec/26 sec)
system.time(
out10.4.4 <- stan(file = "model10_4_4.stan", data = dataList, init = inits,
  warmup = nb, iter = ni, chains = nc, thin = nt)
)
rstan::traceplot(out10.4.4) # not shown
print(out10.4.4, dig = 3) # not shown
stan_est <- summary(out10.4.4)$summary[1:5,1]


#   10.4.5 Do-it-yourself MLEs
# ----------------------------

f <- function(alpha.j, beta.j, mu.alpha, sigma.alpha, mu.beta,
  sigma.beta, sigma.residual, mass.j, lengthN.j){
  mu <- alpha.j + beta.j * lengthN.j
  prod(dnorm(mass.j, mu, sigma.residual)) *
    dnorm(alpha.j, mu.alpha, sigma.alpha) *
    dnorm(beta.j, mu.beta, sigma.beta)
}

f <- function(alpha.j, beta.j, mu.alpha, sigma.alpha, mu.beta,
              sigma.beta, sigma.residual, mass.j, lengthN.j){
  matdim <- dim(alpha.j)
  alpha.j <- as.vector(alpha.j)
  beta.j <- as.vector(beta.j)

  par.mat <- cbind(alpha.j, beta.j)
  Xmat <- cbind(1, lengthN.j)
  mu <- Xmat %*% t(par.mat)
  mass.j.mat <- matrix(mass.j, nrow = nrow(mu), ncol = ncol(mu))
  L <- dnorm(mass.j.mat, mu, sigma.residual)
  L <- apply(L, 2, prod)
  L <- L * dnorm(alpha.j, mu.alpha, sigma.alpha) *
    dnorm(beta.j, mu.beta, sigma.beta)

  matrix(L, nrow = matdim[1], ncol = matdim[2])
}

library(pracma)
g <- function(mu.alpha, mu.beta, sigma.alpha, sigma.beta, sigma.residual, mass.j, lengthN.j){
  tryCatch({
  integral2(fun = f, xmin = 150, xmax = 400, ymin = 0, ymax = 150,
    mu.alpha = mu.alpha, sigma.alpha = sigma.alpha, mu.beta = mu.beta,
    sigma.beta = sigma.beta, sigma.residual = sigma.residual,
    mass.j = mass.j, lengthN.j = lengthN.j)$Q
  }, error = function(e) return(Inf))
}

# Definition of NLL for random-slopes model with Gaussian errors
NLL <- function(pars, data){
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sigma.alpha <- exp(pars[3])
  sigma.beta <- exp(pars[4])
  sigma.residual <- exp(pars[5])
  
  nll <- 0 # initialize nll at 0
  
  for (j in 1:data$nPops){
    mass.j <- data$mass[data$pop == j]
    lengthN.j <- data$lengthN[data$pop == j]
    L <- g(mu.alpha, mu.beta, sigma.alpha, sigma.beta, sigma.residual,
      mass.j, lengthN.j) # likelihood for pop j
  nll <- nll - log(L)
  }

  return(nll)
}

# Definition of NLL for random-slopes model with Gaussian errors
NLL <- function(pars, data){
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  sigma.alpha <- exp(pars[3])
  sigma.beta <- exp(pars[4])
  sigma.residual <- exp(pars[5])
  
  nll <- 0 # initialize nll at 0
  for (j in 1:data$nPops){
    mass.j <- data$mass[data$pop == j]
    lengthN.j <- data$lengthN[data$pop == j]
    L <- g(mu.alpha, mu.beta, sigma.alpha, sigma.beta, sigma.residual,
           mass.j, lengthN.j) # likelihood for pop j
    nll <- nll - log(L)
  }
  
  return(nll)
}

inits <- c(mu.alpha = mean(dataList$mass), mu.beta = 0,
           log.sigma.alpha = 1, log.sigma.beta = 1,
           log.sigma.residual = log(sd(dataList$mass)))

# Optimization (ART 206 sec) and calculation of ASEs and CIs
system.time(
out10.4.5 <- optim(inits, NLL, method = "BFGS", hessian = TRUE, data = dataList,
  control = list(trace = 1, REPORT = 5)) )
get_MLE(out10.4.5, 4)

# Save the MLEs
diy_est <- out10.4.5$par
diy_est[3:5] <- exp(diy_est[3:5]) # Get SD params on ilog scale


#   10.4.6 Likelihood analysis with TMB
# -------------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop-1 # convert pop index to 0-based

# Write TMB model file
cat(file = "model10_4_6.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(mass); //response
  DATA_VECTOR(lengthN); //covariate
  DATA_IVECTOR(pop); //population index
  DATA_INTEGER(nPops); //Number of populations
  DATA_INTEGER(n); //Number of observations

  //Describe parameters
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER(log_sigma_alpha);
  PARAMETER(log_sigma_beta);
  PARAMETER(log_sigma_residual);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);
  
  Type sigma_alpha = exp(log_sigma_alpha);
  Type sigma_beta = exp(log_sigma_beta);
  Type sigma_residual = exp(log_sigma_residual);
  
  Type LL = 0.0; //Initialize log-likelihood at 0

  //Distribution of random intercepts and random slopes
  for (int i = 0; i < nPops; i++){ //Note index starts at 0 instead of 1!
    LL += dnorm(alpha(i), mu_alpha, sigma_alpha, true);
    LL += dnorm(beta(i), mu_beta, sigma_beta, true);
  }

  for (int i = 0; i < n; i++){
    Type mu = alpha(pop(i)) + beta(pop(i)) * lengthN(i);
    //Calculate log-likelihood of observation and add to total
    LL += dnorm(mass(i), mu, sigma_residual, true);
  }
  
  return -LL; //Return negative log likelihood
}
")

# Compile and load TMB function
compile("model10_4_6.cpp")
dyn.load(dynlib("model10_4_6"))

# Provide dimensions and starting values for parameters
params <- list(mu_alpha = 0, mu_beta = 0, log_sigma_alpha = 0,
               log_sigma_beta = 0, log_sigma_residual = 0,
               alpha = rep(0, tmbData$nPops),
               beta = rep(0, tmbData$nPops))

# Create TMB object
out10.4.6 <- MakeADFun(data = tmbData, parameters = params,
  random = c("alpha", "beta"), #Identify which params(s) are random
  DLL = "model10_4_6", silent = TRUE)

# Optimize TMB object and print and save results
opt <- optim(out10.4.6$par, fn = out10.4.6$fn, gr = out10.4.6$gr,
             method = "L-BFGS-B", hessian = TRUE)
(tsum <- tmb_summary(out10.4.6)) # not shown
tmb_est <- c(opt$par[1:2], exp(opt$par[3:5]))


#   10.4.7 Comparison of the parameter estimates
# ----------------------------------------------

# Compare all results with truth
comp <- cbind(truth = truth, gtmbML = gtmbML_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 10.5 The random-coefficients model with correlation between intercept and slope
# -------------------------------------------------------------------------------
# (no code)


#   10.5.1 Introduction
# ---------------------
# (no code)


#   10.5.2 Data generation
# ------------------------

set.seed(10)
nPops <- 56
nSample <- 10
n <- nPops * nSample
pop <- gl(n = nPops, k = nSample)

orig.length <- runif(n, 45, 70) # Body length (cm)
mn <- mean(orig.length)
sd <- sd(orig.length)
cat("Mean and sd used to normalise original length:", mn, sd, "\n\n")
lengthN <- (orig.length - mn)/sd
hist(lengthN, col = "grey")

Xmat <- model.matrix(~ pop * lengthN - 1 - lengthN)
print(Xmat[1:21,], dig = 2) # Print top 21 rows (not shown)

library(MASS) # Load MASS

mu.alpha <- 260 # Values for five hyperparameters
sigma.alpha <- 20
mu.beta <- 60
sigma.beta <- 30
cov.alpha.beta <- -50 # Covariance
sigma.residual <- 30

mu.vector <- c(mu.alpha, mu.beta)
VC.matrix <- matrix(c(sigma.alpha^2, cov.alpha.beta, cov.alpha.beta, sigma.beta^2),2,2)

ranef.matrix <- mvrnorm(n = nPops, mu = mu.vector, Sigma = VC.matrix)
colnames(ranef.matrix) <- c("alpha", "beta")
head(ranef.matrix) # Look at what we’ve created:
# pairs of intercepts and slopes .... THESE are the random effects !
apply(ranef.matrix, 2, mean) # Compare with mu.vector above
var(ranef.matrix) # Compare with VC.matrix above

lin.pred <- Xmat[,] %*% as.vector(ranef.matrix) # Value of lin.predictor
eps <- rnorm(n = n, mean = 0, sd = sigma.residual) # residuals
mass <- lin.pred + eps # response = lin.pred + residual
mass <- drop(mass) # Turn matrix into a vector
hist(mass, col = "grey") # Inspect what we’ve created

# Save true values for comparisons later
truth <- c(mu.alpha = mu.alpha, mu.beta = mu.beta, sigma.alpha = sigma.alpha,
  sigma.beta = sigma.beta, cov.alpha.beta = cov.alpha.beta, sigma.residual = sigma.residual)

library(lattice) # not shown
xyplot(mass ~ lengthN | pop, pch = 16, cex = 1, col = rgb(0,0,0,0.5),
  main = 'Mass-length relationships of asp vipers by population')

# 10.5.3 Likelihood anlaysis using canned functions in R

gtmbData <- data.frame(mass = mass, pop = pop, lengthN = lengthN)
out10.5.ML <- glmmTMB(mass ~ lengthN + (lengthN | pop),
                      data = gtmbData, REML = FALSE) # with ML
out10.5.REML <- glmmTMB(mass ~ lengthN + (lengthN | pop),
                        data = gtmbData, REML = TRUE) # with REML
summary(out10.5.ML)
summary(out10.5.REML)

# Save the estimates from ML and REML
vcML <- summary(out10.5.ML)$varcor$cond$pop
gtmbML_est <- c(fixef(out10.5.ML)$cond, attr(vcML, "stddev"), vcML[2], sigma(out10.5.ML))

vcREML <- summary(out10.5.REML)$varcor$cond$pop
gtmbREML_est <- c(fixef(out10.5.REML)$cond, attr(vcREML, "stddev"), vcREML[2], sigma(out10.5.REML))

# Show covariance
gtmbML_est[5]


#   10.5.4 Bayesian analysis with JAGS
# ------------------------------------

# Bundle and summarize data
str(dataList <- list(mass = mass, pop = as.numeric(pop), lengthN = lengthN, nPops = nPops,
  n = n, V = diag(2), Wdf = 3) )

# Write JAGS model file
cat(file = "model10.5.4.txt", "
model {
# Priors
mu.alpha ~ dnorm(0, 0.00001)
mu.beta ~ dnorm(0,0.00001)

# Put them into a mean vector
Mu[1] <- mu.alpha
Mu[2] <- mu.beta

# Inverse variance-covariance matrix
Tau[1:2, 1:2] ~ dwish(V[1:2, 1:2], Wdf)

# Get variance-covariance matrix
Sigma <- inverse(Tau)
sigma.alpha <- sqrt(Sigma[1,1])
sigma.beta <- sqrt(Sigma[2,2])
cov.alpha.beta <- Sigma[1,2]

sigma.residual ~ dt(0, 1, 0.0001)I(0, )
tau.residual <- pow(sigma.residual, -2)

for (i in 1:nPops){
  B[i,1:2] ~ dmnorm(Mu[], Tau[,])
  alpha[i] <- B[i,1]
  beta[i] <- B[i,2]
}

# Likelihood
for (i in 1:n){
  mass[i] ~ dnorm(mu[i], tau.residual)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthN[i] # Expected value
}}
")

# Function to generate starting values
inits <- function(){ list(mu.alpha = rnorm(1, 0, 1),
  mu.beta = rnorm(1, 0, 1), sigma.residual = rlnorm(1)) }

# Parameters to estimate
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
  "cov.alpha.beta", "sigma.residual", "alpha", "beta")

# MCMC settings
na <- 5000 ; ni <- 10000 ; nb <- 5000 ; nc <- 4 ; nt <- 1

# Call JAGS from R (ART 1 min), check convergence, summarize posteriors and save estimates
out10.5.4 <- jags(dataList, inits, params, "model10.5.4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out10.5.4) # not shown
print(out10.5.4, 2)
jags_est <- out10.5.4$summary[c(1:6),1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, gtmbML = gtmbML_est, JAGS = jags_est)
print(comp, 4)


#   10.5.5 Bayesian analysis with NIMBLE
# --------------------------------------

# Summarize data (same as for JAGS)
str(dataList) # not shown

# Write NIMBLE model file
model10.5.5 <- nimbleCode( {
# Priors
mu.alpha ~ dnorm(0, 0.00001)
mu.beta ~ dnorm(0,0.00001)

# Put them into a mean vector
Mu[1] <- mu.alpha
Mu[2] <- mu.beta

# Inverse variance-covariance matrix
Tau[1:2, 1:2] ~ dwish(V[1:2, 1:2], Wdf)

# Get variance-covariance matrix
Sigma[1:2,1:2] <- inverse(Tau[1:2,1:2])

sigma.alpha <- sqrt(Sigma[1,1])
sigma.beta <- sqrt(Sigma[2,2])
cov.alpha.beta <- Sigma[1,2]

sigma.residual ~ T(dt(0, 1, 0.0001), 0, )
tau.residual <- pow(sigma.residual, -2)

for (i in 1:nPops){
  B[i,1:2] ~ dmnorm(Mu[1:2], Tau[1:2,1:2])
  alpha[i] <- B[i,1]
  beta[i] <- B[i,2]
}

# Likelihood
for (i in 1:n){
  mass[i] ~ dnorm(mu[i], tau.residual)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthN[i] #Expected value
}
} )

# Inits
inits <- function(){ list(mu.alpha = rnorm(1, 0, 1), 
  mu.beta = rnorm(1, 0, 1), sigma.residual = rlnorm(1)) }

# Parameters monitored: same as before
params <- c("mu.alpha", "mu.beta", "sigma.alpha", "sigma.beta",
            "cov.alpha.beta", "sigma.residual", "alpha", "beta")

# MCMC settings
ni <- 10000  ;  nb <- 5000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 68 sec), check convergence, summarize posteriors and save estimates
system.time( out10.5.5 <- 
    nimbleMCMC(code = model10.5.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out10.5.5)
(nsum <- nimble_summary(out10.5.5, params))   # not shown
nimble_est <- nsum[1:6,1]


#   10.5.6 Bayesian analysis with Stan
# ------------------------------------

# Summarize data set again
str(dataList)

# Write Stan model
cat(file = "model10_5_6.stan", "
data {
  int n; //Number of observations
  int nPops; //Number of populations
  vector[n] mass; //Response variable
  vector[n] lengthN; //Covariate
  array[n] int pop; //Population assignment of each obs
  cov_matrix[2] V; //Wishart scale matrix
  int Wdf; //Wishart degrees of freedom
}

parameters {
  real mu_alpha;
  real mu_beta;
  cov_matrix[2] Sigma;
  real <lower = 0> sigma_residual;
  matrix[nPops, 2] B; //Matrix of random effects
}

transformed parameters {
  vector[2] Mu;
  real sigma_alpha;
  real sigma_beta;
  real cov_alpha_beta;
  vector[nPops] alpha;
  vector[nPops] beta;

  Mu[1] = mu_alpha;
  Mu[2] = mu_beta;
  sigma_alpha = sqrt(Sigma[1,1]);
  sigma_beta = sqrt(Sigma[2,2]);
  cov_alpha_beta = Sigma[1,2];
  alpha = B[1:nPops,1];
  beta = B[1:nPops,2];
}

model {
  vector[n] mu;

  mu_alpha ~ normal(0, 1000);
  mu_beta ~ normal(0, 1000);
  Sigma ~ inv_wishart(Wdf, V);
  sigma_residual ~ cauchy(0, 100);

  for (i in 1:nPops){
    B[i,1:2] ~ multi_normal(Mu, Sigma);
  }

  for (i in 1:n){
    mu[i] = alpha[pop[i]] + beta[pop[i]] * lengthN[i];
    mass[i] ~ normal(mu[i], sigma_residual);
  }
}
")

# Function to generate starting values
inits <- function(){ list(mu_alpha = rnorm(1, 0, 1),
  mu_beta = rnorm(1, 0, 1), sigma_residual = rlnorm(1))}

# Parameters monitored: we use hyphens as separators
params <- c("mu_alpha", "mu_beta", "sigma_alpha", "sigma_beta",
            "cov_alpha_beta", "sigma_residual", "alpha", "beta")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 220 sec/130 sec)
system.time(
out10.5.6 <- stan(file = "model10_5_6.stan", data = dataList,
  pars = params, warmup = nb,iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out10.5.6) # not shown
print(out10.5.6, dig = 3) # not shown
stan_est <- summary(out10.5.6)$summary[1:6,1]


#   10.5.7 Do-it-yourself MLEs
# ----------------------------

library(mvtnorm)
f <- function(alpha.j, beta.j, Mu, Sigma,
  sigma.residual, mass.j, lengthN.j){
  mu <- alpha.j + beta.j * lengthN.j
  B.j <- c(alpha.j, beta.j)
  prod(dnorm(mass.j, mu, sigma.residual)) *
  dmvnorm(B.j, Mu, Sigma)
}

f <- function(alpha.j, beta.j, Mu, Sigma,
              sigma.residual, mass.j, lengthN.j){
  matdim <- dim(alpha.j)
  alpha.j <- as.vector(alpha.j)
  beta.j <- as.vector(beta.j)

  par.mat <- cbind(alpha.j, beta.j)
  Xmat <- cbind(1, lengthN.j)
  mu <- Xmat %*% t(par.mat)
  mass.j.mat <- matrix(mass.j, nrow = nrow(mu), ncol = ncol(mu))
  L <- dnorm(mass.j.mat, mu, sigma.residual)
  L <- apply(L, 2, prod)
  L <- L * dmvnorm(par.mat, Mu, Sigma)
  
  matrix(L, nrow = matdim[1], ncol = matdim[2])
}

library(pracma)
g <- function(Mu, Sigma, sigma.residual, mass.j, lengthN.j){
  tryCatch({
  integral2(fun = f, xmin = 150, xmax = 400, ymin = 0, ymax = 150,
            Mu = Mu, Sigma = Sigma, sigma.residual = sigma.residual,
            mass.j = mass.j, lengthN.j = lengthN.j)$Q
  }, error = function(e) return(Inf))
}

# Definition of NLL
NLL <- function(pars, data){
  mu.alpha <- pars[1]
  mu.beta <- pars[2]
  Beta <- c(mu.alpha, mu.beta) # mean vector
  
  theta <- pars[3:5]
  # Make valid covariance matrix from unconstrained theta vector
  # Log-Cholesky parameterization
  # Sigma = M'M, where M is upper tri matrix with positive diagonal
  M <- diag(exp(theta[1:2])) # diagonal elements must be positive
  M[1,2] <- theta[3]
  Sigma <- t(M)%*%M

  sigma.residual <- exp(pars[6])
  
  LL <- 0 # initialize LL at 0

  for (j in 1:data$nPops){
    mass.j <- data$mass[data$pop == j]
    lengthN.j <- data$lengthN[data$pop == j]
    L <- g(Beta, Sigma, sigma.residual,
          mass.j, lengthN.j) # likelihood for pop j
    LL <- LL + log(L)
  }
return(-LL)
}

# Minimize that NLL to find MLEs and get SEs (ART 420 sec)
inits <- c(mu_alpha = mean(dataList$mass), mu_beta = 0, theta1 = 1, theta2 = 1,
  theta3 = 0, log_sigma_residual = log(sd(dataList$mass)))

system.time(
  out10.5.7 <- optim(inits, NLL, hessian = TRUE, data = dataList,
                     control = list(trace = 1, REPORT = 5, maxit = 1000)) )
get_MLE(out10.5.7, 4)

# Recover Sigma values from theta via log-Cholesky parameterization
theta <- out10.5.7$par[3:5]
M <- matrix(0, nrow = 2, ncol = 2)
M <- diag(exp(theta[1:2]))
M[1,2] <- theta[3]
(Sigma <- t(M)%*%M)

# Save estimates
diy_est <- c(out10.5.7$par[1:2], sqrt(diag(Sigma)), Sigma[1,2],
             exp(out10.5.7$par[6]))


#   10.5.8 Likelihood analysis with TMB
# -------------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop-1 # convert pop index to 0-based
str(tmbData) # not shown

# Write TMB model file
cat(file = "model10_5_8.cpp",
"#include <TMB.hpp>
using namespace density;
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(mass); //response
  DATA_VECTOR(lengthN); //covariate
  DATA_IVECTOR(pop); //population index
  DATA_INTEGER(nPops); //Number of populations
  DATA_INTEGER(n); //Number of observations

  //Describe parameters
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);
  PARAMETER_VECTOR(theta); //Parameters associated with VC matrix
  PARAMETER(log_sigma_residual);
  PARAMETER_MATRIX(B); //Random effects matrix

  //Make valid covariance matrix from theta
  //using log-Cholesky parameterization
  matrix <Type> Sigma(2,2);
  matrix <Type> M(2,2);
  M(0,0) = exp(theta(0));
  M(1,1) = exp(theta(1));
  M(0,1) = theta(2);
  M(1,0) = 0;
  Sigma = M.transpose() * M;

  Type sigma_residual = exp(log_sigma_residual);
  Type nll = 0.0; //Initialize negative log-likelihood at 0

  //Multivariate normal prior on correlated random slopes and intercepts
  for (int i = 0; i < nPops; i++){
    nll += MVNORM(Sigma)(B.row(i)); //MVNORM returns -loglik
  }

  //Add random effects to mean values to get complete intercepts/slopes
  vector <Type> alpha = B.col(0);
  alpha += mu_alpha;
  vector <Type> beta = B.col(1);
  beta += mu_beta;
  for (int i = 0; i < n; i++){
    Type mu = alpha(pop(i)) + beta(pop(i)) * lengthN(i);
    //Calculate log-likelihood of observation and add to total
    nll -= dnorm(mass(i), mu, sigma_residual, true);
  }
  ADREPORT(Sigma);
  ADREPORT(alpha); ADREPORT(beta);
  return nll; //Return negative log likelihood
}
")

# Compile and load TMB function
compile("model10_5_8.cpp")
dyn.load(dynlib("model10_5_8"))

# Provide dimensions and starting values for parameters
params <- list(mu_alpha = mean(tmbData$mass), mu_beta = 0,
               theta = c(1, 1, 0),
               log_sigma_residual = log(sd(tmbData$mass)),
               B = matrix(0, nPops, 2))

# Create TMB object
out10.5.8 <- MakeADFun(data = tmbData, parameters = params,
             random = "B", #Identify which params(s) are random
             DLL = "model10_5_8", silent = TRUE)

# Optimize TMB object and print results
opt <- optim(out10.5.8$par, fn = out10.5.8$fn, gr = out10.5.8$gr,
  method = "L-BFGS-B", hessian = TRUE)
(tsum <- tmb_summary(out10.5.8))

# Recover VC matrix and save results
# Recover Sigma values via same parameterization as in the DIY-MLEs
theta <- opt$par[3:5]
M <- matrix(0, nrow = 2, ncol = 2)
M <- diag(exp(theta[1:2]))
M[1,2] <- theta[3]
Sigma <- t(M)%*%M

tmb_est <- c(opt$par[1:2], sqrt(diag(Sigma)), Sigma[1,2],
             exp(opt$par[6]))


#   10.5.9 Comparison of the parameter estimates
# ----------------------------------------------

# Compare all results with truth
comp <- cbind(truth = truth, gtmbML = gtmbML_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 10.6 Summary and outlook
# ------------------------
# (no code)
