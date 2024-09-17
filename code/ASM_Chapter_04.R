
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -----------------------------------------------------------------------
# Chapter 4  --  Introduction to general-purpose model fitting engines 
#                and the "model of the mean"
# -----------------------------------------------------------------------

# Last changes: 11 June 2024


# 4.1 Introduction
# ----------------
# (no code)


# 4.2 Data generation
# -------------------

# Generate two samples of body mass measurements of male peregrines
set.seed(4)
y10 <- rnorm(n = 10, mean = 600, sd = 30) # Sample of 10 birds
y1000 <- rnorm(n = 1000, mean = 600, sd = 30) # Sample of 1000 birds
# Save the data-generating values of the parameters for later comparisons
truth <- c(mean = 600, sd = 30)
# Plot data (Fig. 4-2)
xlim = c(450, 750)
par(mfrow = c(1, 2), mar = c(6, 6, 6, 3), cex.lab = 1.5, cex.axis = 1.5)
hist(y10, col = 'grey ', xlim = xlim, main = 'Body mass (g) of 10 male peregrines')
hist(y1000, col = 'grey', xlim = xlim, main = 'Body mass (g) of 1000 male peregrines')


# 4.3 Analysis using canned functions in R
# ----------------------------------------

summary(out4.3 <- lm(y10 ~ 1)) # small data set: not shown
summary(out4.3 <- lm(y1000 ~ 1)) # large data set

# Save estimates
lm_est <- c(coef(out4.3), sigma = sigma(out4.3))


# 4.4 JAGS
# --------
# (no code)

#   4.4.1 Introduction to JAGS
# ----------------------------
# (no code)


#   4.4.2 Fit the model with JAGS
# -------------------------------

library(jagsUI)

# Bundle and summarize data
str(dataList <- list(mass = y1000, n = length(y1000)))

# Write JAGS model file
cat(file = "model4.4.txt", " # This code line is R
model { # Starting here, we have BUGS code
# Priors
pop.mean ~ dunif(0, 5000) # Population mean
precision <- 1 / pop.var # Precision = 1/variance
pop.var <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# Likelihood
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
} # This is the last line of BUGS code
") # . . . and this is R again

# Function to generate starting values
inits <- function(){
  list(pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))
}

# Parameters monitored
params <- c("pop.mean", "pop.sd", "pop.var")

# MCMC settings
na <- 1000 # Number of iterations in the adaptive phase
ni <- 12000 # Number of draws from the posterior (in each chain)
nb <- 2000 # Number of draws to discard as burn-in
nc <- 4 # Number of chains
nt <- 1 # Thinning rate (nt = 1 means we do not thin)

# Call JAGS (ART 1 min) and marvel at JAGS' progress bar
out4.4 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model4.4.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = FALSE)

# Call JAGS in parallel (ART <1 min) and check convergence
out4.4 <- jags(data = dataList, inits = inits, parameters.to.save = params,
  model.file = "model4.4.txt", n.iter = ni, n.burnin = nb, n.chains = nc,
  n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out4.4) # Produce Fig. 4-3

print(out4.4, 3) # Produce a summary of the fitted model object

names(out4.4)

hist(out4.4$summary[,8]) # Rhat values in column 8 of the summary
which(out4.4$summary[,8] > 1.1) # check for non-convergence: none here

summary(out4.4$samples) # not shown

par(mfrow = c(3,1), mar = c(5,5,4,2)) # not shown
pars <- colnames(out4.4$samples[[1]])[1:3]
for (i in pars){
  matplot(as.array(out4.4$samples)[,i,], type = 'l', lty = 1,
          col = c("green", "red", "blue", "orange"),
          xlab = "Iteration", ylab = "Value", main = i, frame = FALSE)
}

par(mfrow = c(3, 2), mar = c(5,5,4,2)) # Fig. 4-4
for (i in pars){
  samps <- as.matrix(out4.4$samples[,i])
  hist(samps, col = 'grey', breaks = 100, xlab = i, main = "")
  plot(density(samps), type = 'l', lwd = 2, col = 'gray40', main = "",
       xlab = i, frame = FALSE)
}

sims <- as.matrix(out4.4$samples)
head(sims[,"pop.mean"] < 597) # show first of MANY logical tests!
mean(sims[,"pop.mean"] < 597) # Prob(mu < 597)

# Fig. 4–5: compute 2d probability game and plot bivariate posterior
test.true <- sims[,"pop.mean"] < 600 & sims[,"pop.sd"] > 30
mean(test.true)
par(mfrow = c(1,1))
plot(sims[,"pop.mean"], sims[,"pop.sd"], pch = 16, col = rgb(0,0,0,0.3),
cex = 0.8, frame = FALSE)
points(sims[test.true,"pop.mean"], sims[test.true,"pop.sd"], pch = 16,
col = rgb(1,0,0,0.6), cex = 0.8)
abline(h = 30, col = 'red')
abline(v = 600, col = 'red')

# alternatively do this to see bivariate posterior plots
pairs(sims[,1:3])

apply(as.matrix(out4.4$samples[,1:2]), 2, summary)

apply(as.matrix(out4.4$samples[,1:2]), 2, sd)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- out4.4$summary[1:2,1]
print(cbind(truth = truth, lm = lm_est, JAGS = jags_est), 5)


# 4.5 NIMBLE
# ----------
# (no code)

#   4.5.1 Introduction to NIMBLE
# ------------------------------
# (no code)


#   4.5.2 Fit the model with NIMBLE
# ---------------------------------

# Load NIMBLE
library(nimble)

# Bundle data (same as for JAGS)
str(dataList <- list(mass = y1000, n = length(y1000))) # not shown

# Write Nimble model file
model4.5 <- nimbleCode( {
# Priors and linear models
pop.mean ~ dunif(0, 5000) # Normal parameterized by precision
precision <- 1 / pop.variance # Precision = 1/variance
pop.variance <- pop.sd * pop.sd
pop.sd ~ dunif(0, 100)

# Likelihood
for(i in 1:n){
  mass[i] ~ dnorm(pop.mean, precision)
}
} )

# Can use same function to generate starting values as for JAGS
inits <- function()
  list (pop.mean = rnorm(1, 600), pop.sd = runif(1, 1, 30))

# Parameters monitored: same as before
params <- c("pop.mean", "pop.sd", "pop.variance")

# MCMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call NIMBLE (ART 30 sec)
system.time(
out4.5 <- nimbleMCMC(code = model4.5, constants = dataList,
  inits = inits, monitors = params, nburnin = nb, niter = ni,
  thin = nt, nchains = nc, samplesAsCodaMCMC = TRUE, summary = TRUE) )

lapply(out4.5, head) # print first 5 values of each chain (not shown)

par(mfrow = c(1, 3));
  coda::traceplot(out4.5$samples) # not shown
out4.5$summary$all.chains # Posterior summaries from all 4 chains

library(ASMbook)
nsum <- nimble_summary(out4.5$samples, params) # Summary table
round(nsum, 3)
nimble_est <- nsum[1:2,1] # save estimates

# Create a NIMBLE model from BUGS code
rawModel <- nimbleModel(code = model4.5, constants = dataList, inits = inits())

# Configure the MCMC algorithm: create a default MCMC configuration
# and monitor our selected list of parameters
mcmcConfig <- configureMCMC(rawModel, monitors = params)

# Build the MCMC algorithm function
rawMCMC <- buildMCMC(mcmcConfig)

system.time(rawMCMC$run(10)) # Take 10 samples, takes about 10 sec
as.matrix(rawMCMC$mvSamples) # View samples, not shown

compModel <- compileNimble(rawModel)

compMCMC <- compileNimble(rawMCMC, project = rawModel)

system.time(compMCMC$run(10)) # not shown, takes about 0.02 sec

system.time(
samples <- runMCMC(compMCMC, niter = ni, nburnin = nb, thin = nt,
  nchains = nc, samplesAsCodaMCMC = TRUE) )

# Peak at samples
lapply(samples, head) # not shown

# Produce marginal posterior summaries
(nsum <- nimble_summary(samples)) # not shown

# Traceplots
par(mfrow = c(2,2)); coda::traceplot(samples) # not shown


# 4.6 Stan
# --------
# (no code)


#   4.6.1 Introduction to Stan
# ----------------------------
# (no code)


#   4.6.2 Fit the model with Stan
# -------------------------------

# Load Stan R package
library(rstan)

# Bundle and summarize data
str(dataList <- list(n = length(y1000), mass = y1000))

# Write text file with model description in Stan language
cat(file = "model4_6.stan", # This line is still R code
"data { // This is the first line of Stan code
  int n; // Define the format of all data
  vector[n] mass; // . . .including the dimension of vectors
}

parameters { // Same for the parameters
  real pop_mean;
  real <lower = 0> pop_sd;
}

model {
  // Priors
  pop_mean ~ normal(0, 1000);
  pop_sd ~ cauchy(0, 10);

  // Likelihood
  for (i in 1:n){
    mass[i] ~ normal(pop_mean, pop_sd);
  }
} // This is the last line of Stan code
" )

# HMC settings
ni <- 1200 ; nb <- 200 ; nc <- 4 ; nt <- 1

# Call STAN (ART 33 / 3 sec)
system.time(
out4.6 <- stan(file = "model4_6.stan", data = dataList,
  chains = nc, iter = ni, warmup = nb, thin = nt) )

# Check convergence
rstan::traceplot(out4.6) # not shown

# Print posterior summaries
print(out4.6)

# Save estimates
stan_est <- summary(out4.6)$summary[1:2,1]

str(out4.6) # not shown

sims <- extract(out4.6)
lapply(sims, head) # not shown

sims <- As.mcmc.list(out4.6) # Note capital 'A' !
lapply(sims, head) # not shown
par(mfrow = c(2,2)); coda::traceplot(sims) # not shown


# 4.7 Maximum likelihood in R
# ---------------------------
# (no code)


#   4.7.1 Introduction to maximum likelihood in R
# -----------------------------------------------
# (no code)


#   4.7.2 Fit the model using maximum likelihood in R
# ---------------------------------------------------

dnorm(y10[1], mean = 600, sd = 30)
dnorm(y10[1], mean = 550, sd = 25)

prod(dnorm(y10, mean = 600, sd = 30))

prod(dnorm(y10, mean = 550, sd = 25))

sum(dnorm(y10, mean = 600, sd = 30, log = TRUE))

# Definition of LL for a normal linear model with intercept only
LL <- function(param, y) {
  pop.mean <- param[1] # First parameter is mean (mu)
  pop.sd <- param[2] # Second is SD (sigma)
  sum(dnorm(y, mean = pop.mean, sd = pop.sd, log = TRUE))
}

LL(c(600, 30), y10)

starts <- c(pop.mean = 550, pop.sd = 25)

# Definition of NLL for a normal linear model with constant mean
NLL <- function(param, y) {
  pop.mean <- param[1] # First parameter is mean (mu)
  pop.sd <- param[2] # Second is SD (sigma)
  -sum(dnorm(y, mean = pop.mean,
    sd = pop.sd, log = TRUE)) # Note minus sign
}

out4.7 <- optim(starts, fn = NLL, hessian = TRUE, y = y1000)

out4.7$convergence

out4.7$par
diy_est <- out4.7$par # save estimates

vcov <- solve(out4.7$hessian) # invert Hessian
ASE <- sqrt(diag(vcov)) # Get asymptotic SEs
ASE

get_MLE(out4.7, 5)


# 4.8 Maximum likelihood using Template Model Builder
# ---------------------------------------------------
# (no code)


#   4.8.1 Introduction to Template Model Builder
# ----------------------------------------------
# (no code)


#   4.8.2 Fit the model using TMB
# -------------------------------

library(TMB)

# Bundle data
str(dataList <- list(y = y1000, n = length(y1000)))

# Write TMB model file
cat(file = "model4_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type> ::operator() ()
{
  // Describe input data
  DATA_INTEGER(n); // Sample size
  DATA_VECTOR(y); // observations
  
  // Describe parameters
  PARAMETER(pop_mean); // Mean
  PARAMETER(log_pop_sd); // log(standard deviation)
  Type pop_sd = exp(log_pop_sd); // Type = match type of function output

  Type LL = 0.0; // Initialize total log likelihood at 0
  for (int i = 0; i < n; i ++ ){ // Note index starts at 0 instead of 1!
    // Calculate log-likelihood of observation
    // value of true = function should return log(lik)
    LL += dnorm(y(i), pop_mean, pop_sd, true);
  }
  return -LL; // Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model4_8.cpp") # Produces gibberish . . . have no fear !
dyn.load(dynlib("model4_8"))

# Provide dimensions and starting values for parameters
params <- list(pop_mean = 0, log_pop_sd = 0)

# Create TMB object
out4.8 <- MakeADFun(data = dataList, parameters = params, random = NULL,
  DLL = "model4_8", silent = TRUE)

# Optimize TMB object
opt <- optim(out4.8$par, fn = out4.8$fn, gr = out4.8$gr,
  method = "BFGS", hessian = TRUE)

opt$par

summary(sdreport(out4.8))

tsum <- tmb_summary(out4.8) # not shown
tmb_est <- tsum[,1] # save results
tmb_est[2] <- exp(tmb_est[2]) # convert SD from log scale


# 4.9 Comparison of engines and concluding remarks
# ------------------------------------------------

comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 5)
