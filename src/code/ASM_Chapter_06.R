
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -----------------------------------------------------
# Chapter 6  --  Comparing two groups in a normal model
# -----------------------------------------------------

# Last changes: 11 June 2024


# 6.1 Introduction
# ----------------
# (no code)

# 6.2 Comparing two groups with equal variances
# ---------------------------------------------
# (no code)

#   6.2.1 Data generation
# -----------------------

# Generate a data set
set.seed(61)
n1 <- 60                             # Number of females
n2 <- 40                             # Number of males
mu1 <- 105                           # Population mean of females
mu2 <- 77.5                          # Population mean of males
sigma <- 2.75                        # Average population SD of both
n <- n1+ n2                          # Total sample size
y1 <- rnorm(n1, mu1, sigma)          # Data for females
y2 <- rnorm(n2, mu2, sigma)          # Date for males
y <- c(y1, y2)                       # Merge both data sets
x <- rep(c(0,1), c(n1, n2))          # Indicator variable indexing a male

# Make a plot (Fig. 6-1)
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1, frame = FALSE)

# Alternative and equivalent data generation
set.seed(61)
n <- n1 + n2                         # Total sample size
alpha <- mu1                         # Mean for females serves as the intercept
beta <- mu2 - mu1                    # beta is the difference male-female
E.y <- alpha + beta*x                # Expectation (linear predictor)
y.obs <- rnorm(n = n, mean = E.y, sd = sigma) # Add random variation
boxplot(y.obs ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1, frame = FALSE)

# Save true values for later comparisons
truth <- c(alpha = mu1, beta = mu2 - mu1, sigma = sigma)

# Load required libraries
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)


#   6.2.2 Likelihood analysis with canned functions in R
# ------------------------------------------------------

summary(out6.2.2 <- lm(y ~ x))
lm_est <- c(coef(out6.2.2), sigma(out6.2.2)) # save results

t.test(y ~ x, var.equal = TRUE)

xfac <- factor(rep(c("group1", "group2"), c(n1, n2)))
head(xfac) # look at the first few values

summary(lm(y ~ xfac)) # not shown (effects or treatment contrast parameterization)
summary(lm(y ~ xfac - 1)) # not shown (means parameterization)

# Diagnostic checks of residuals/model GOF (not shown)
plot(out6.2.2)                       # Traditional Pearson residual check
simOut <- simulateResiduals(out6.2.2, n = 1000, plot = TRUE)
plotResiduals(simOut)                # same as right plot . . .
testDispersion(out6.2.2)             # Test for over- or underdispersion


#   6.2.3 Bayesian analysis with JAGS
# -----------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, x = x, n = n))

# Write JAGS model file
cat(file = "model6.2.3.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)             # Intercept (=female mean)
beta ~ dnorm(0, 0.0001)              # Slope (=diff between females and males)
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)                 # common SD

# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta * x[i]       # linear predictor
}

# Derived parameters
mu1 <- alpha                         # female mean
mu2 <- alpha + beta                  # male mean
for (i in 1:n){                      # residuals
  residual[i] <- y[i] - mu[i]
}
}")

# Function to generate starting values
inits <- function(){list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta", "sigma", "mu1", "mu2", "residual")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out6.2.3 <- jags(data = dataList, inits = inits, parameters.to.save = params,
model.file = "model6.2.3.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out6.2.3) # not shown
print(out6.2.3, 3)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- out6.2.3$summary[1:3,1]
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


# Residual plots
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(1:100, out6.2.3$mean$residual) # not shown
abline(h = 0)
boxplot(out6.2.3$mean$residual ~ x, col = "grey", xlab = "Male",
  ylab = "Wingspan residuals (cm)", las = 1) # not shown
abline(h = 0)


#   6.2.4 Bayesian analysis with NIMBLE
# -------------------------------------

library(nimble)

# Bundle and summarize data (same as for JAGS)
dataList <- list(y = y, x = x, n = n)


# Write NIMBLE model file
model6.2.4 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, 0.0001)             # Intercept (=female mean)
beta ~ dnorm(0, 0.0001)              # Slope (=diff between females and males)
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)                 # common SD

# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta*x[i]         # linear predictor
}

# Derived quantities
mu1 <- alpha                         # female mean
mu2 <- alpha + beta                  # male mean
for (i in 1:n){                      # residuals
  residual[i] <- y[i] - mu[i]
}
} )

# Function to generate starting values
inits <- function(){list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta", "sigma", "mu1", "mu2", "residual")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 30 sec), check convergence, summarize posteriors and save estimates
system.time(  
  out6.2.4 <- nimbleMCMC(code = model6.2.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out6.2.4) # not shown
(nsum <- nimble_summary(out6.2.4, params))   # not shown
nimble_est <- nsum[1:3,1]                    # Save estimates


#   6.2.5 Bayesian analysis with Stan
# -----------------------------------

# Bundle and summarize data (same as before)
dataList <- list(y = y, x = x, n = n)

# Write text file with model description in BUGS language
# Version 1 with loops over reponse
cat(file = "model6_2_5.stan",        # This line is R code
"data {                              // This is the first line of Stan code
  int<lower=0> n;                    // Define the format of all data
  vector[n] y;                       // ... including the dimension of vectors
  vector[n] x;                       //
}

parameters {                         // Define format for all parameters
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters{
  vector[n] mu;
  for (i in 1:n){
    mu[i] = alpha + beta * x[i];     // Calculate linear predictor
  }
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  // Likelihood
  for (i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
}

generated quantities {
  vector[n] residuals;
  for (i in 1:n){
    residuals[i] = y[i] - mu[i];
  }
}                                    // This is the last line of Stan code
")


# Write text file with model description in BUGS language
# Version 2 (vectorized)
cat(file = "model6_2_5.stan",        # This line is R code
"data {                              // This is the first line of Stan code
  int<lower=0> n;                    // Define the format of all data
  vector[n] y;                       // ... including the dimension of vectors
  vector[n] x;                       //
}

parameters {                         // Define format for all parameters
  real alpha;
  real beta;
  real<lower=0> sigma;
}

transformed parameters{
  vector[n] mu;
  mu = alpha + beta * x;             // Calculate linear predictor
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  // Likelihood
  y ~ normal(mu, sigma);
}

generated quantities {
  vector[n] residuals;
  residuals = y - mu;
}                                    // This is the last line of Stan code
")

# HMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 57/4 sec), check convergence, summarize posteriors and save estimates
system.time(
  out6.2.5 <- stan(file = "model6_2_5.stan", data = dataList,
    chains = nc, iter = ni, warmup = nb, thin = nt) )
rstan::traceplot(out6.2.5)           # not shown

print(out6.2.5, dig = 2)             # not shown
stan_est <- summary(out6.2.5)$summary[1:3,1]


#   6.2.6 Do-it-yourself MLEs
# ---------------------------

NLL <- function(params, y, x) {
  alpha <- params[1]
  beta <- params[2]
  sigma <- exp(params[3])            # convert from log scale
  n <- length(y)                     # number of datapoints
  mu <- LL <- numeric(n)             # empty vectors
  # You could vectorize this loop to speed things up
  for (i in 1:n){
    mu[i] <- alpha + beta * x[i]
    LL[i] <- dnorm(y[i], mean = mu[i], sd = sigma, log = TRUE)
  }
  -sum(LL)
}

# Minimize that NLL to find MLEs and also get SEs
inits <- c('alpha' = 50, 'beta' = 10, 'log_sigma' = 0)
out6.2.6 <- optim(inits, NLL, y = y, x = x, hessian=TRUE)
diy_est <- c(out6.2.6$par[1:2], exp(out6.2.6$par[3])) # save estimates
get_MLE(out6.2.6, 5)


#   6.2.7 Likelihood analysis with TMB
# ------------------------------------

# Bundle and summarize data (same as before)
str(dataList <- list(y = y, x = x, n = n))   # not shown

# Write TMB model file
cat(file = "model6_2_7.cpp",
"#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(y);                    // response
  DATA_VECTOR(x);                    // covariate
  DATA_INTEGER(n);                   // Number of obs

  //Describe parameters
  PARAMETER(alpha);                  // Intercept
  PARAMETER(beta);                   // Slope
  PARAMETER(log_sigma);              // log(residual standard deviation)

  Type sigma = exp(log_sigma);       // Type = match type of function output (double)
  Type LL = 0.0;                     // Initialize total log likelihood at 0
  Type mu;

  for (int i= 0; i<n; i++){          // Note index starts at 0 instead of 1!
    mu = alpha + beta * x(i);
    // Calculate log-likelihood of observation and add to total
    LL += dnorm(y(i), mu, sigma, true); // Add log-lik of obs i
  }

  // Derived parameters
  Type mu1 = alpha;
  Type mu2 = alpha + beta;
  ADREPORT(mu1);                     // save mu1 and mu2 to output
  ADREPORT(mu2);
  return -LL;                        // Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model6_2_7.cpp")
dyn.load(dynlib("model6_2_7"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0, log_sigma = 0)

# Create TMB object
out6.2.7 <- MakeADFun(data = dataList, parameters = params,
DLL = "model6_2_7", silent = TRUE)

# Optimize TMB object, print results and save estimates
opt <- optim(out6.2.7$par, fn = out6.2.7$fn, gr = out6.2.7$gr, method = "BFGS", hessian = TRUE)
tmb_est <- c(opt$par[1:2], exp(opt$par[3]))  # save estimates
(tsum <- tmb_summary(out6.2.7))              # look at output (not shown)


#   6.2.8 Comparison of the parameter estimates
# ---------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est, NIMBLE = nimble_est,
  Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 6.3 Comparing two groups with unequal variances
# -----------------------------------------------
# (no code)


#   6.3.1 Data generation
# -----------------------

set.seed(63)
# Generate data set
n1 <- 60                             # Number of females
n2 <- 40                             # Number of males
mu1 <- 105                           # Population mean for females
mu2 <- 77.5                          # Population mean for males
sigma1 <- 3                          # Population SD for females
sigma2 <- 2.5                        # Population SD for males
n <- n1+ n2                          # Total sample size
y1 <- rnorm(n1, mu1, sigma1)         # Data for females
y2 <- rnorm(n2, mu2, sigma2)         # Data for males
y <- c(y1, y2)                       # Merge both data sets
x <- rep(c(0,1), c(n1, n2))          # Indicator for male

truth <- c(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)

# Make a plot (Fig. 6-2)
par(mfrow = c(1, 1), mar = c(6,6,6,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
boxplot(y ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1, frame = FALSE)


#   6.3.2 Frequentist analysis with canned functions in R
# -------------------------------------------------------

# Welch test
(out6.3.2 <- t.test(y ~ x))
tt_est <- c(out6.3.2$est, NA, NA) # save estimates of means

# Diagnostic checks of residuals/model
summary(fm <- lm(y ~ as.factor(x)-1) ) # Fit model from 6.2.2
simOut <- simulateResiduals(fm, n = 1000, plot = TRUE)


#   6.3.3 Bayesian analysis with JAGS
# -----------------------------------

# Bundle and summarize data
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# Write JAGS model file
cat(file = "model6.3.3.txt", "
model {
# Priors
mu1 ~ dnorm(0, 0.0001)
mu2 ~ dnorm(0, 0.0001)
tau1 <- pow(sigma1, -2)
sigma1 ~ dunif(0, 1000)
tau2 <- pow(sigma2, -2)
sigma2 ~ dunif(0, 1000)

# Likelihood
for (i in 1:n1) {                    # First sample (females)
  y1[i] ~ dnorm(mu1, tau1)
}

for (i in 1:n2) {                    # Second sample (males)
  y2[i] ~ dnorm(mu2, tau2)
}

# Derived quantities
delta <- mu2 - mu1
}
")

# Function to generate starting values
inits <- function(){ list(mu1 = rnorm(1), mu2 = rnorm(1), sigma1 = rlnorm(1), sigma2 = rlnorm(1))}

# Parameters to estimate
params <- c("mu1", "mu2", "delta", "sigma1", "sigma2")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out6.3.3 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model6.3.3.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out6.3.3) # not shown
jags_est <- out6.3.3$summary[c(1, 2, 4, 5), 1] # save estimates
print(out6.3.3, 3)

samps <- as.matrix(out6.3.3$samples)
quantile(samps[,'sigma2'] - samps[,'sigma1'], c(0.025,0.975))


#   6.3.4 Bayesian analysis with NIMBLE
# -------------------------------------

# Bundle and summarize data
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# Write NIMBLE model file
model63.4 <- nimbleCode( {
# Priors
mu1 ~ dnorm(0, 0.0001)
mu2 ~ dnorm(0, 0.0001)
tau1 <- pow(sigma1, -2)
sigma1 ~ dunif(0, 1000)
tau2 <- pow(sigma2, -2)
sigma2 ~ dunif(0, 1000)

# Likelihood
for (i in 1:n1) {                    # First sample (females)
  y1[i] ~ dnorm(mu1, tau1) 
}

for (i in 1:n2) {                    # Second sample (males)
  y2[i] ~ dnorm(mu2, tau2) 
}

# Derived quantities
delta <- mu2 - mu1
} )


# Inits
inits <- function(){ list(mu1 = rnorm(1), mu2 = rnorm(1), sigma1 = rlnorm(1), sigma2 = rlnorm(1))}

# Parameters monitored: same as before
params <- c("mu1", "mu2", "delta", "sigma1", "sigma2")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 30 sec), check convergence, summarize posteriors and save estimates
system.time(  
  out6.3.4 <- nimbleMCMC(code = model63.4,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out6.3.4)     # not shown
(nsum <- nimble_summary(out6.3.4, params))       # not shown
nimble_est <- nsum[c(1,2,4,5),1] # save estimates


#   6.3.5 Bayesian analysis with Stan
# -----------------------------------

# Bundle and summarize data
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# Write text file with model description in BUGS language
cat(file = "model6_3_5.stan",        # This line is R code
"data {                              //This is the first line of Stan code
  int n1;                            //Define the format of all data
  int n2;
  vector[n1] y1;
  vector[n2] y2;
}

parameters {                         //Define format for all parameters
  real mu1;
  real mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

model {                              //Priors and likelihood here
  //Priors
  mu1 ~ normal(0, 100);
  mu2 ~ normal(0, 100);
  sigma1 ~ cauchy(0, 10);
  sigma2 ~ cauchy(0, 10);
  //Likelihood
  y1 ~ normal(mu1, sigma1);
  y2 ~ normal(mu2, sigma2);
}

generated quantities{
  real delta = mu2 - mu1;
}                                    //This is the last line of Stan code
" )

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 55/1 sec)
system.time(
  out6.3.5 <- stan(file = "model6_3_5.stan", data = dataList,
    chains = nc, iter = ni, warmup = nb, thin = nt) )

rstan::traceplot(out6.3.5) # not shown
stan_est <- summary(out6.3.5)$summary[1:4, 1] # save estimates
print(out6.3.5, dig = 2)             # not shown


#   6.3.6 Do-it-yourself MLEs
# ---------------------------

NLL <- function(param, y1, y2){
  mu1 <- param[1]
  mu2 <- param[2]
  sigma1 <- exp(param[3])
  sigma2 <- exp(param[4])
  LL1 <- dnorm(y1, mu1, sigma1, log = TRUE) # females
  LL2 <- dnorm(y2, mu2, sigma2, log = TRUE) # males
  LL <- sum(LL1) + sum(LL2)                 # total
  -LL                                       # negative log likelihood
}

# Minimize the NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('mu1' = 50, 'mu2' = 50, 'log_sigma1' = 1, 'log_sigma2' = 1)
out6.3.6 <- optim(inits, NLL, y1 = y1, y2 = y2, hessian = TRUE)
diy_est <- c(out6.3.6$par[1:2], exp(out6.3.6$par[3:4]))
get_MLE(out6.3.6, 5)


#   6.3.7 Likelihood analysis with TMB
# ------------------------------------

# Bundle and summarize data (same as before)
str(dataList <- list(y1 = y1, y2 = y2, n1 = n1, n2 = n2))

# Write TMB model file
cat(file = "model6_3_7.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(y1);                   //response pop1
  DATA_VECTOR(y2);                   //response pop2
  DATA_INTEGER(n1);                  //Number of obs pop1
  DATA_INTEGER(n2);                  //Number of obs pop2
  
  //Describe parameters
  PARAMETER(mu1);                    //Intercept pop1
  PARAMETER(mu2);                    //Intercept pop2
  PARAMETER(log_sigma1);             //log(residual standard deviation) pop1
  PARAMETER(log_sigma2);             //log(residual standard deviation) pop2

  Type sigma1 = exp(log_sigma1);
  Type sigma2 = exp(log_sigma2);
  Type LL = 0.0;                     //Initialize total log likelihood at 0
  
  for (int i= 0; i<n1; i++){         //Note index starts at 0 instead of 1!
    LL += dnorm(y1(i), mu1, sigma1, true); //Add log-lik of obs i
  }

  for (int i= 0; i<n2; i++){         //Note index starts at 0 instead of 1!
    LL += dnorm(y2(i), mu2, sigma2, true); //Add log-lik of obs i
  }
  
  Type delta = mu2 - mu1;
  ADREPORT(delta);

  return -LL;                        //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model6_3_7.cpp")
dyn.load(dynlib("model6_3_7"))

# Provide dimensions and starting values for parameters
params <- list(mu1 = 0, mu2 = 0, log_sigma1 = 0, log_sigma2 = 0)

# Create TMB object
out6.3.7 <- MakeADFun(data = dataList, parameters = params,
  DLL = "model6_3_7", silent = TRUE)

# Optimize TMB object, print results and save estimates
opt <- optim(out6.3.7$par, fn = out6.3.7$fn, gr = out6.3.7$gr,
method = "BFGS", hessian = TRUE)
tmb_est <- c(opt$par[1:2], exp(opt$par[3:4])) # save estimates
(tsum <- tmb_summary(out6.3.7))


#   6.3.8 Comparison of the parameter estimates
# ---------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, tt = tt_est, JAGS = jags_est, NIMBLE = nimble_est,
 Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)


# 6.4 Summary and a comment on the modeling of variances
# ------------------------------------------------------
# (no code)
