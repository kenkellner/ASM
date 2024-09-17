
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ---------------------------------
# Chapter 20  --  Integrated models
# ---------------------------------

# Last changes: 12 June 2024


# 20.1 Introduction
# -----------------
# (no code)


# 20.2 Data generation: simulating three abundance data sets 
#      with different observation/aggregation models
# ----------------------------------------------------------

set.seed(20)

# Simulate the two data sets and plot them
# Choose sample size and parameter values for both data sets
nsites1 <- 500                                         # Sample size for count data
nsites2 <- 1000                                        # Sample size for zero-truncated counts
nsites3 <- 2000                                        # Sample size for detection/nondetection data
mean.lam <- 2                                          # Average expected abundance (lambda) per site
beta <- -2                                             # Coefficient of elevation covariate on lambda
truth <- c("log.lam" = log(mean.lam),"beta" = beta)    # Save truth

# Simulate elevation covariate for all three and standardize to mean of 1000 and
# standard deviation also of 1000 m
elev1 <- sort(runif(nsites1, 200, 2000))               # Imagine 200–2000 m a.s.l.
elev2 <- sort(runif(nsites2, 200, 2000))
elev3 <- sort(runif(nsites3, 200, 2000))
selev1 <- (elev1-1000)/1000                            # Scaled elev1
selev2 <- (elev2-1000)/1000                            # Scaled elev2
selev3 <- (elev3-1000)/1000                            # Scaled elev3

# Create three regular count data sets with log-linear effects
C1 <- rpois(nsites1, exp(log(mean.lam) + beta * selev1))
C2 <- rpois(nsites2, exp(log(mean.lam) + beta * selev2))
C3 <- rpois(nsites3, exp(log(mean.lam) + beta * selev3))

table(C1)                                              # Tabulate data set 1

# Create data set 2 (C2) by zero-truncating (discard all zeroes)
ztC2 <- C2                                             # Make a copy
ztC2 <- ztC2[ztC2 > 0]                                 # Tossing out zeroes yields zero-truncated data

table(C2); table(ztC2)                                 # tabulate both original and ZT data set

# Turn count data set 3 (C3) into detection/nondetection data (y)
y <- C3                                                # Make a copy
y[y > 1] <- 1                                          # Squash to binary

table(C3) ; table(y)                                   # tabulate both original counts and DND

# Plot counts/DND data in all data sets (Fig. 20.3)
library(scales)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1), cex = 1.2, cex.lab = 1.5, cex.axis = 1.5, las = 1)
plot(elev2[C2>0], jitter(ztC2), pch = 16, xlab = 'Elevation (m)', ylab = 'Counts',
frame = FALSE, ylim = range(c(C1, ztC2)), col = alpha('grey80', 1), main = 'Data sets 1 and 2')
points(elev1, jitter(C1), pch = 16)
lines(200:2000, exp(log(2) -2 * ((200:2000)-1000)/1000 ), col = 'red', lty = 1, lwd = 2)
axis(1, at = c(250, 750, 1250, 1750), tcl = -0.25, labels = NA)
plot(elev3, jitter(y, amount = 0.04), xlab = 'Elevation (m)', ylab = 'Detection/nondetection',
axes = FALSE, pch = 16, col = alpha('grey60', 0.3), main = 'Data set 3')
axis(1)
axis(1, at = c(250, 750, 1250, 1750), tcl = -0.25, labels = NA)
axis(2, at = c(0, 1), labels = c(0, 1))

# Required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 20.3 Fitting models to the three individual data sets first
# -----------------------------------------------------------
# (no code)


#   20.3.1 Fitting a standard Poisson generalized linear model in R 
#          and JAGS to data set 1
# -----------------------------------------------------------------

# Get MLEs for individual data sets
# Data set 1: Poisson GLM with log link for counts

## Canned function in R
# ---------------------
summary(fm1 <- glm(C1 ~ selev1, family = poisson(link = "log")))
exp(coef(fm1)[1])                                      # Estimate of lambda on natural scale from counts

## DIY-MLEs
# ---------
# Definition of negative log-likelihood (NLL) for Poisson GLM
NLL1 <- function(param, y, x) {
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x)
  LL <- dpois(y, lambda, log = TRUE)                   # LL contribution for each datum
  return(-sum(LL))                                     # NLL for all data
}

# Minimize NLL
inits <- c(alpha = 0, beta = 0)                        # Need to provide initial values
sol1 <- optim(inits, NLL1, y = C1, x = selev1, hessian = TRUE, method = 'BFGS')
tmp1 <- get_MLE(sol1, 4)


## JAGS
# ------

# Bundle data
str(bdata <- list(C1 = C1, nsites1 = nsites1, selev1 = selev1))

# Write JAGS model file
cat(file = "model1.txt", "
model {
# Priors and linear models
alpha ~ dunif(-10, 10)                                 # Abundance intercept on log scale
mean.lam <- exp(alpha)                                 # Abundance intercept on natural scale
beta ~ dnorm(0, 0.0001)                                # Slope on elevation

# Likelihood for data set 1
for (i in 1:nsites1){                                  # Data set 1
  C1[i] ~ dpois(lambda1[i])
  log(lambda1[i]) <- alpha + beta * selev1[i]
}
}
")

# Initial values
inits <- function(){list(alpha = runif(1), beta = rnorm(1))}

# Parameters monitored
params <- c("mean.lam", "alpha", "beta")

# MCMC settings
ni <- 6000; nb <- 2000; nc <- 4; nt <- 4; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "model1.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out1)                                # Not shown
print(out1, 4)

# Compare truth with likelihood and Bayesian solutions
comp <- cbind(truth = truth, tmp1[,1:2], out1$summary[2:3, 1:2])
colnames(comp)[3:5] <- c("SE(MLE)", "Post.mean", "Post.sd")
print(comp, 3)


#   20.3.2 Fitting the zero-truncated Poisson generalized linear model
#          in R and JAGS to data set 2
# -------------------------------------------------------------------

## DIY-MLEs
# ---------
# Definition of negative log-likelihood (NLL) for the ztPoisson GLM
NLL2 <- function(param, y, x) {
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x)
  L <- dpois(y, lambda)/(1-ppois(0, lambda))           # L. contribution for each datum
  return(-sum(log(L)))                                 # NLL for all data
}

# Minimize NLL
inits <- c(alpha = 0, beta = 0)                        # Need to provide initial values
sol2 <- optim(inits, NLL2, y = ztC2, x = selev2[C2>0], hessian = TRUE, method = 'BFGS')
tmp2 <- get_MLE(sol2, 4)

# Minimize 'wrong' Poisson NLL which ignores the zero truncation
sol2a <- optim(inits, NLL1, y = ztC2, x = selev2[C2>0], hessian = TRUE, method = 'BFGS')
get_MLE(sol2a, 4)

# Compare intercepts of right and wrong Poisson NLL on the natural scale
exp(sol2$par[1]) ; exp(sol2a$par[1])


## JAGS
# -----

# Bundle data
str(bdata <- list(C2 = ztC2, nsites2 = length(ztC2), selev2 = selev2[C2 > 0]))

# Write JAGS model file
cat(file = "model2.txt", "
model {
# Priors and linear models
alpha ~ dunif(-10, 10)                                 # Abundance intercept on log scale
mean.lam <- exp(alpha)                                 # Abundance intercept on natural scale
beta ~ dnorm(0, 0.0001)                                # Slope on elevation

# Zero-truncated Poisson likelihood for data set 2
for (j in 1:nsites2){                                  # Data set 2
  C2[j] ~ dpois(lambda1[j])T(1,)                       # truncation is accommodated easily !
  log(lambda1[j]) <- alpha + beta * selev2[j]
}
}
")

# Initial values
inits <- function(){list(alpha = runif(1), beta = rnorm(1))}

# Parameters monitored
params <- c("mean.lam", "alpha", "beta")

# MCMC settings
ni <- 6000; nb <- 2000; nc <- 4; nt <- 4; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out2 <- jags(bdata, inits, params, "model2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out2)                                # Not shown
print(out2, 2)

# Compare truth with likelihood and Bayesian solutions
comp <- cbind(truth = truth, tmp2[,1:2], out2$summary[2:3, 1:2])
colnames(comp)[3:5] <- c("SE(MLE)", "Post.mean", "Post.sd")
print(comp, 3)


#   20.3.3 Fitting the cloglog Bernoulli generalized linear model in 
#          R and JAGS to data set 3
# ------------------------------------------------------------------

## Canned function in R
# ---------------------
# Data set 3: Bernoulli GLM with cloglog link for detection/nondetection
summary(fm3 <- glm(y ~ selev3, family = binomial(link = "cloglog")))
exp(coef(fm3)[1])                                      # Estimate of lambda on natural scale from binary data


## DIY-MLEs
# ---------
# Definition of NLL for Bernoulli GLM with cloglog link
NLL3 <- function(param, y, x) {
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x)
  psi <- 1-exp(-lambda)
  LL <- dbinom(y, 1, psi, log = TRUE)                  # L. contribution for each datum
  return(-sum(LL))                                     # NLL for all data
}

# Minimize NLL
inits <- c(alpha = 0, beta = 0)
sol3 <- optim(inits, NLL3, y = y, x = selev3, hessian = TRUE, method = 'BFGS')
tmp3 <- get_MLE(sol3, 4)


## JAGS
# -----
# Bundle data
str(bdata <- list(y = y, nsites3 = nsites3, selev3 = selev3))

# Write JAGS model file
cat(file = "model3.txt", "
model {
# Priors and linear models
alpha ~ dunif(-10, 10)                                 # Abundance intercept on log scale
mean.lam <- exp(alpha)                                 # Abundance intercept on natural scale
beta ~ dnorm(0, 0.0001)                                # Slope on elevation

# Likelihood for data set 3
for (k in 1:nsites3){                                  # Data set 3
  y[k] ~ dbern(psi[k])
  cloglog(psi[k]) <- alpha + beta * selev3[k]
  # Alternative implementation of same model for data set 3
  # y[k] ~ dbern(psi[k])
  # psi[k] <- 1 - exp(-lambda3[k])
  # log(lambda3[k]) <- alpha + beta * selev3[k]
}
}
")

# Initial values
inits <- function(){list(alpha = runif(1), beta = rnorm(1))}

# Parameters monitored
params <- c("mean.lam", "alpha", "beta")

# MCMC settings
ni <- 6000; nb <- 2000; nc <- 4; nt <- 4; na <- 1000

# Call JAGS from R (ART 1.2 min), check convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "model3.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out3)                                # Not shown
print(out3, 3)

# Compare truth with likelihood and Bayesian solutions
comp <- cbind(truth = truth, tmp3[,1:2], out3$summary[2:3, 1:2])
colnames(comp)[3:5] <- c("SE(MLE)", "Post.mean", "Post.sd")
print(comp, 3)


# 20.4 Fitting the IM to all three data sets simultaneously
# ---------------------------------------------------------
# (no code)

#   20.4.1 Fitting the IM with a DIY likelihood function
# ------------------------------------------------------

# Definition of the joint NLL for the integrated model
NLLjoint <- function(param, y1, x1, y2, x2, y3, x3) {
  # Definition of elements in param vector (shared between data sets)
  alpha <- param[1]                                    # log-linear intercept
  beta <- param[2]                                     # log-linear slope
  # Likelihood for the Poisson GLM for data set 1 (y1, x1)
  lambda1 <- exp(alpha + beta * x1)
  L1 <- dpois(y1, lambda1)
  # Likelihood for the ztPoisson GLM for data set 2 (y2, x2)
  lambda2 <- exp(alpha + beta * x2)
  L2 <- dpois(y2, lambda2)/(1-ppois(0, lambda2))
  # Likelihood for the cloglog Bernoulli GLM for data set 3 (y3, x3)
  lambda3 <- exp(alpha + beta * x3)
  psi <- 1-exp(-lambda3)
  L3 <- dbinom(y3, 1, psi)
  # Joint log-likelihood and joint NLL: here you can see that sum!
  JointLL <- sum(log(L1)) + sum(log(L2)) + sum(log(L3)) # Joint LL
  return(-JointLL)                                      # Return joint NLL
}

# Minimize NLLjoint
inits <- c(alpha = 0, beta = 0)
solJoint <- optim(inits, NLLjoint, y1 = C1, x1 = selev1, y2 = ztC2, x2 = selev2[C2> 0],
  y3 = y, x3 = selev3, hessian = TRUE, method = 'BFGS')

# Get MLE and asymptotic SE and print and save
(tmp4 <- get_MLE(solJoint, 4))
diy_est <- tmp4[,1]


#   20.4.2 Fitting the IM with JAGS
# ---------------------------------

# Bundle data
str(dataList <- list(C1 = C1, C2 = ztC2, y = y, nsites1 = nsites1, nsites2 = length(ztC2),
  nsites3 = nsites3, selev1 = selev1, selev2 = selev2[C2>0], selev3 = selev3))

# Write JAGS model file
cat(file = "model4.txt", "
model {
# Priors and linear models: shared for models of all three data sets
alpha ~ dunif(-10, 10)                                 # Abundance intercept on log scale
mean.lam <- exp(alpha)                                 # Abundance intercept on natural scale
beta ~ dnorm(0, 0.0001)                                # Slope on elevation

# Joint likelihood: Note identical alpha and beta for all data sets
# Likelihood portion for data set 1: regular counts
for (i in 1:nsites1){
  C1[i] ~ dpois(lambda1[i])
  log(lambda1[i]) <- alpha + beta * selev1[i]
}
# Likelihood portion for data set 2: zero-truncated counts
for (j in 1:nsites2){
  C2[j] ~ dpois(lambda2[j])T(1,)
  log(lambda2[j]) <- alpha + beta * selev2[j]
}
# Likelihood portion for data set 3: detection/nondetection
for (k in 1:nsites3){
  y[k] ~ dbern(psi[k])
  cloglog(psi[k]) <- alpha + beta * selev3[k]
}
}
")

# Initial values
inits <- function(){list(alpha = runif(1), beta = rnorm(1))}

# Parameters monitored
params <- c("mean.lam", "alpha", "beta")

# MCMC settings
na <- 1000; ni <- 6000; nb <- 2000; nc <- 4; nt <- 4

# Call JAGS from R (ART 170 sec), check convergence and summarize posteriors
out4 <- jags(dataList, inits, params, "model4.txt", n.iter = ni, n.burnin = nb,
n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out4)
print(out4, 2)
jags_est <- out4$summary[2:3,1]

# Compare truth with likelihood and Bayesian solutions
comp <- cbind(truth = truth, tmp4[,1:2], out4$summary[2:3, 1:2])
colnames(comp)[3:5] <- c("SE(MLE)", "Post.mean", "Post.sd")
print(comp, 3)


#   20.4.3 Fitting the IM with NIMBLE
# -----------------------------------

library(nimble)

# Bundle data
str(dataList <- list(C1 = C1, C2 = ztC2, y = y, nsites1 = nsites1, nsites2 = length(ztC2),
  nsites3 = nsites3, selev1 = selev1, selev2 = selev2[C2>0], selev3 = selev3) )

# Write NIMBLE model file
model20.4.3 <- nimbleCode( {
# Priors and linear models: shared for models of all three data sets
alpha ~ dunif(-10, 10)                                 # Abundance intercept on log scale
mean.lam <- exp(alpha)                                 # Abundance intercept on natural scale
beta ~ dnorm(0, sd = 100)                              # Slope on elevation

# Joint likelihood: Note identical alpha and beta for all data sets
# Likelihood portion for data set 1: regular counts
for (i in 1:nsites1){
  C1[i] ~ dpois(lambda1[i])
  log(lambda1[i]) <- alpha + beta * selev1[i]
}
# Likelihood portion for data set 2: zero-truncated counts
for (j in 1:nsites2){
  C2[j] ~ T(dpois(lambda2[j]), 1, )
  log(lambda2[j]) <- alpha + beta * selev2[j]
}
# Likelihood portion for data set 3: detection/nondetection
for (k in 1:nsites3){
  y[k] ~ dbern(psi[k])
  cloglog(psi[k]) <- alpha + beta * selev3[k]
}} )

# Inits
inits <- function(){list(alpha = runif(1), beta = rnorm(1))}

# Parameters monitored
params <- c("mean.lam", "alpha", "beta")

# MCMC settings
ni <- 6000; nb <- 2000; nc <- 4; nt <- 4; na <- 1000

# Call NIMBLE (ART 90 sec), check convergence, summarize posteriors and save estimates
system.time( out20.4.3 <- 
    nimbleMCMC(code = model20.4.3,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out20.4.3)          # not shown
(nsum <- nimble_summary(out20.4.3, params))            # not shown
nimble_est <- nsum[2:3, 1]                             # Save estimates


#   20.4.4 Fitting the IM with Stan
# ---------------------------------

cat(file = "model4.stan",
"data {
  int nsites1;
  int nsites2;
  int nsites3;
  array[nsites1] int C1;
  array[nsites2] int C2;
  array[nsites3] int y;
  vector[nsites1] selev1;
  vector[nsites2] selev2;
  vector[nsites3] selev3;
}
parameters {
  real alpha;
  real beta;
}
model {
  vector[nsites1] lambda1;
  vector[nsites2] lambda2;
  vector[nsites3] psi;

  // Priors
  alpha ~ uniform(-10, 10);
  beta ~ normal(0, 100);

  // Likelihood
  for (i in 1:nsites1){
    lambda1[i] = exp(alpha + beta * selev1[i]);
    C1[i] ~ poisson(lambda1[i]);
  }
  for (j in 1:nsites2){
    lambda2[j] = exp(alpha + beta * selev2[j]);
    C2[j] ~ poisson(lambda2[j]) T[1, ];
  }
  for (k in 1:nsites3){
    psi[k] = inv_cloglog(alpha + beta * selev3[k]);
    y[k] ~ bernoulli(psi[k]);
  }
}

generated quantities {
  real mean_lam = exp(alpha);
}
" )

# Parameters monitored
params <- c("mean_lam", "alpha", "beta")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 90/45 sec), assess convergence and print results table
system.time(
  out4.stan <- stan(file = "model4.stan", data = dataList,
    warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out4.stan) # not shown
print(out4.stan, dig = 3) # not shown
stan_est <- summary(out4.stan)$summary[1:2,1]


#   20.4.5 Fitting the IM with TMB
# --------------------------------

cat(file = "model4.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_INTEGER(nsites1);
  DATA_INTEGER(nsites2);
  DATA_INTEGER(nsites3);
  DATA_VECTOR(C1);
  DATA_VECTOR(C2);
  DATA_VECTOR(y);
  DATA_VECTOR(selev1);
  DATA_VECTOR(selev2);
  DATA_VECTOR(selev3);

  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);

  Type LL = 0.0;                                       //Initialize log-likelihood at 0

  vector<Type> lambda1(nsites1);
  vector<Type> lambda2(nsites2);
  vector<Type> psi(nsites3);

  for (int i= 0; i<nsites1; i++){
    lambda1(i) = exp(alpha + beta * selev1(i));
    LL += dpois(C1(i), lambda1(i), true);
  }
  for (int j = 0; j< nsites2; j++){
    lambda2(j) = exp(alpha + beta * selev2(j));
    // Truncated Poisson
    LL += log(dpois(C2(j), lambda2(j))/(1 - ppois(Type(0), lambda2(j))));
  }
  for (int k= 0; k<nsites3; k++){
    psi(k) = 1 - exp(-exp(alpha + beta * selev3(k))); //inverse cloglog
    LL += dbinom(y(k), Type(1), psi(k), true);
  }

  Type mean_lam = exp(alpha);
  ADREPORT(mean_lam);

  return -LL;
}
")

# Compile and load TMB function
compile("model4.cpp")
dyn.load(dynlib("model4"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0)

# Create TMB object
out4.tmb <- MakeADFun(data = dataList, parameters = params,
DLL = "model4", silent = TRUE)

# Optimize TMB object, print and save results
opt <- optim(out4.tmb$par, fn = out4.tmb$fn, gr = out4.tmb$gr, method = "BFGS")
(tsum <- tmb_summary(out4.tmb))
tmb_est <- tsum[1:2,1]


#   20.4.6 Comparison of the parameter estimates for the IM
# ---------------------------------------------------------

# Compare point estimates from the five engines
comp <- cbind(truth = truth, DIY = diy_est, JAGS = jags_est, NIMBLE = nimble_est,
  Stan = stan_est, TMB = tmb_est)
print(comp, 4)


# 20.5 What do we gain by analyzing the joint likelihood in our analysis?
# -----------------------------------------------------------------------

# Compare truth with MLEs only from all 4 models (stacked sideways)
print(cbind(truth = truth, "MLE(Poisson)"= tmp1[,1], "MLE(ZTPois)" = tmp2[,1],
  "MLE(cloglogBern)" = tmp3[,1], "MLE(integrated)" = tmp4[,1]), 3)

# Compare ASEs from all 4 models (stacked sideways)
print(cbind("ASE(Poisson)" = tmp1[,2], "ASE(ZTPois)" = tmp2[,2], "ASE(cloglogBern)" = tmp3[,2],
  "ASE(integrated)" = tmp4[,2]), 3)


# Compare MLEs and ASEs from all 4 models (Fig. 20.4)
par(mfrow = c(1, 2), mar = c(12, 6, 5, 3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.8)
# Plot for abundance intercept (on log scale)
all.mles <- c(tmp1[1,1], tmp2[1,1], tmp3[1,1], tmp4[1,1])
all.lower.CL <- c(tmp1[1,3], tmp2[1,3], tmp3[1,3], tmp4[1,3])
all.upper.CL <- c(tmp1[1,4], tmp2[1,4], tmp3[1,4], tmp4[1,4])
plot(1:4, all.mles, pch = 16, xlab = '', ylab = '', axes = FALSE, frame = FALSE,
  main = 'Abundance intercept', cex = 2, ylim = c(0.5, 0.9))
axis(1, at = 1:4, c('Model 1', 'Model 2', 'Model 3', 'Integrated model'), las = 2)
segments(1:4, all.lower.CL, 1:4, all.upper.CL, lwd = 1.5)
axis(2, las = 1)
abline(h = log(2), lwd = 1, lty = 3)

# Plot for abundance slope (on log scale)
all.mles <- c(tmp1[2,1], tmp2[2,1], tmp3[2,1], tmp4[2,1])
all.lower.CL <- c(tmp1[2,3], tmp2[2,3], tmp3[2,3], tmp4[2,3])
all.upper.CL <- c(tmp1[2,4], tmp2[2,4], tmp3[2,4], tmp4[2,4])
plot(1:4, all.mles, pch = 16, xlab = '', ylab = '', axes = FALSE, frame = FALSE,
  main = 'Abundance slope', cex = 2, ylim = c(-2.2, -1.8))
axis(1, at = 1:4, c('Model 1', 'Model 2', 'Model 3', 'Integrated model'), las = 2)
segments(1:4, all.lower.CL, 1:4, all.upper.CL, lwd = 1.5)
axis(2, las = 1)
abline(h = -2, lwd = 1, lty = 3)


# 20.6 Summary and outlook
# ------------------------
# (no code)
