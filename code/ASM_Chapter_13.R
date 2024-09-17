
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Chapter 13  --  Poisson GLM with continuous and categorical 
#                 explanatory variables
# ---------------------------------------------------------------

# Last changes: 11 June 2024


# 13.1 Introduction
# -----------------
# (no code)


# 13.2 Data generation
# --------------------

set.seed(13)
nPops <- 3 # Three populations
nSample <- 100 # . . .with 100 individuals each
n <- nPops * nSample # Total sample size
x <- rep(1:nPops, rep(nSample, nPops)) # Population indicator
pop <- factor(x, labels = c("Pyrenees", "Massif Central", "Jura"))
orig.length <- runif(n, 4.5, 7.0) # Wing length (cm)
wing.length <- orig.length - mean(orig.length) # Center values

Xmat <- model.matrix(~ pop * wing.length)
head(Xmat, 10) # Look at first 10 rows of design matrix

# Save truth for comparisons
truth <- beta.vec <- c(-2, 1, 2, 4, -2, -5)

lin.pred <- Xmat[,] %*% beta.vec # Value of lin.predictor
lambda <- exp(lin.pred) # Poisson mean: expected count
load <- rpois(n = n, lambda = lambda) # Add Poisson noise

# Inspect what we’ve created (Fig. 13.2)
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.axis = 1.5, cex.lab = 1.5)
hist(load, col = "grey", breaks = 30, xlab = "Parasite load", main = "", las = 1)
plot(wing.length, load, pch = rep(c("P", "M", "J"), each = nSample), las = 1,
col = rep(c("Red", "Green", "Blue"), each = nSample), ylab = "Parasite load",
xlab = "Wing length", cex = 1.2, frame = FALSE)

# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 13.3 Likelihood analysis with canned functions in R
# ---------------------------------------------------

# Fit model using least-squares and save estimates
summary(out13.3 <- glm(load ~ pop * wing.length, family = poisson))
glm_est <- coef(out13.3)

# Model assessment by checks of residuals (not shown)
plot(out13.3) # Check of traditional (Pearson) residuals
library(DHARMa) # Compute quantile residuals based on simulation
simOut <- simulateResiduals(out13.3, n = 1000, plot = TRUE) # Check


# 13.4 Bayesian analysis with JAGS
# --------------------------------
# (no code)


#   13.4.1 Fitting the model
# --------------------------

# Bundle and summarize data
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops,
wing.length = wing.length, n = n) )

# Write JAGS model file
cat(file = "model13.4.txt", "
model {
# Priors
for (i in 1:nPops){ # Loop over 3 populations
  alpha[i] ~ dnorm(0, 0.0001) # Intercepts
  beta[i] ~ dnorm(0, 0.0001) # Slopes
}
# Likelihood
for (i in 1:n) { # Loop over all 300 data points
  load[i] ~ dpois(lambda[i]) # The response variable (C above)
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]]* wing.length[i])
} # Note nested indexing: alpha[pop[i]] and beta[pop[i]]
# Derived quantities
# Recover effects relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1] # Intercept Massif Central vs. Pyr.
a.effe3 <- alpha[3] - alpha[1] # Intercept Jura vs. Pyr.
b.effe2 <- beta[2] - beta[1] # Slope Massif Central vs. Pyr.
b.effe3 <- beta[3] - beta[1] # Slope Jura vs. Pyr.
# Custom test
test1 <- beta[3] - beta[2] # Slope Jura vs. Massif Central
}
")

# Function to generate starting values
inits <- function(){list(alpha = rlnorm(nPops, 3, 1), beta = rlnorm(nPops, 2, 1))}

# Parameters to estimate
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")
# MCMC settings
na <- 2000; ni <- 25000; nb <- 5000; nc <- 4; nt <- 5
# Call JAGS (ART <1 min), check convergence and summarize posteriors
out13.4 <- jags(dataList, inits, params, "model13.4.txt", n.iter = ni,
n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out13.4) # not shown
print(out13.4, 3)

# Compare glm() estimes with Bayesian estimates and with truth
jags_est <- unlist(out13.4$mean)[c(1, 7, 8, 4, 9, 10)]
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est)
print(comp, 4)


#   13.4.2 Forming predictions
# ----------------------------

# Create a vector with 100 wing lengths
orig.wlength <- sort(orig.length)
wlength <- orig.wlength - mean(orig.length)
# Create matrices to contain prediction for each winglength and MCMC iteration
nsamp <- out13.4$mcmc.info$n.samples # Get size of posterior sample
mite.load.Pyr <- mite.load.MC <- mite.load.Ju <- array(dim = c(nsamp, 300))
# Fill in these vectors: this is clumsy, but it works
alpha <- out13.4$sims.list$alpha # posterior of alpha
beta <- out13.4$sims.list$beta # posterior of beta
for (i in 1:300){
  mite.load.Pyr[,i] <- exp(alpha[,1] + beta[,1] * wlength[i])
  mite.load.MC[,i] <- exp(alpha[,2] + beta[,2] * wlength[i])
  mite.load.Ju[,i] <- exp(alpha[,3] + beta[,3] * wlength[i])
}

# Compute 95% Bayesian credible intervals around predictions
LCB.Pyr <- apply(mite.load.Pyr, 2, quantile, prob = 0.025)
UCB.Pyr <- apply(mite.load.Pyr, 2, quantile, prob = 0.975)
LCB.MC <- apply(mite.load.MC, 2, quantile, prob = 0.025)
UCB.MC <- apply(mite.load.MC, 2, quantile, prob = 0.975)
LCB.Ju <- apply(mite.load.Ju, 2, quantile, prob = 0.025)
UCB.Ju <- apply(mite.load.Ju, 2, quantile, prob = 0.975)

# Compute posterior means
mean.rel <- cbind(exp(out13.4$mean$alpha[1] + out13.4$mean$beta[1] * wlength),
exp(out13.4$mean$alpha[2] + out13.4$mean$beta[2] * wlength),
exp(out13.4$mean$alpha[3] + out13.4$mean$beta[3] * wlength))
covar <- cbind(orig.wlength, orig.wlength, orig.wlength)

# Plot (Fig. 13.3)
par(mar = c(6, 6, 5, 3), cex.lab = 1.5, cex.axis = 1.5)
matplot(orig.wlength, mean.rel, col = c("red", "green", "blue"), type = "l",
lty = 1, lwd = 2, las = 1, ylab = "Expected mite load",
xlab = "Wing length (cm)", frame = FALSE, ylim = c(0, 25))
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.Pyr, rev(UCB.Pyr)),
col = rgb(1,0,0,0.2), border = NA)
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.MC, rev(UCB.MC)),
col = rgb(0,1,0,0.2), border = NA)
polygon(c(orig.wlength, rev(orig.wlength)), c(LCB.Ju, rev(UCB.Ju)),
col = rgb(0,0,1,0.2), border = NA)
matplot(orig.wlength, mean.rel, col = c("red", "green", "blue"),
type = "l", lty = 1, lwd = 3, add = TRUE)


# 13.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops, 
                     wing.length = wing.length, n = n) )

# Write NIMBLE model file
model13.5 <- nimbleCode( {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, sd = 100)  # Intercepts
  beta[i] ~ dnorm(0, sd = 100)   # Slopes
}

# Likelihood
for (i in 1:n) {
  load[i] ~ dpois(lambda[i])     # The response variable
  lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]]* wing.length[i])
}   # Note nested indexing: alpha[pop[i]] and beta[pop[i]]

# Derived quantities
# Recover effects relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1]   # Intercept Massif Central vs. Pyr.
a.effe3 <- alpha[3] - alpha[1]   # Intercept Jura vs. Pyr.
b.effe2 <- beta[2] - beta[1]     # Slope Massif Central vs. Pyr.
b.effe3 <- beta[3] - beta[1]     # Slope Jura vs. Pyr.

# Custom comparison
test1 <- beta[3] - beta[2]       # Slope Jura vs. Massif Central
} )

# Inits
inits <- function(){list(alpha=rlnorm(nPops, 3, 1), beta = rlnorm(nPops, 2, 1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMC settings
ni <- 25000  ;  nb <- 5000  ; nc <- 4  ; nt <- 5

# Call NIMBLE (ART 40 sec), check convergence and summarize posteriors
system.time( out13.5 <- 
    nimbleMCMC(code = model13.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out13.5)
(nsum <- nimble_summary(out13.5, params))
nimble_est <- nsum[c(1,7,8,4,9,10),1]    # Save param estimates


# 13.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(load = load, pop = as.numeric(pop), nPops = nPops,
wing_length = wing.length, n = n) )

# Write Stan model
cat(file = "model13_6.stan", "
data {
  int n;
  int nPops;
  array[n] int load;
  vector[n] wing_length;
  array[n] int pop;
}

parameters {
  vector[nPops] alpha;
  vector[nPops] beta;
}

model {
  vector[n] lambda;

  for (i in 1:nPops){
    alpha[i] ~ normal(0, 100);
    beta[i] ~ normal(0, 100);
  }

  for (i in 1:n) {
    lambda[i] = exp(alpha[pop[i]] + beta[pop[i]] * wing_length[i]);
    load[i] ~ poisson(lambda[i]);
  }
}

generated quantities {
  real a_effe2 = alpha[2] - alpha[1];
  real a_effe3 = alpha[3] - alpha[1];
  real b_effe2 = beta[2] - beta[1];
  real b_effe3 = beta[3] - beta[1];
  real test1 = beta[3] - beta[2];
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 60/10 sec), assess convergence and print results table
system.time(
out13.6 <- stan(file = "model13_6.stan", data = dataList,
warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out13.6) # not shown
print(out13.6, dig = 3) # not shown
stan_est <- summary(out13.6)$summary[c(1,7,8,4,9,10),1] # Save estimates

# 13.7 Do-it-yourself maximum likelihood estimates

head(Xmat, 10) # Look at first 10 rows of design matrix

# Define NLL for general Poisson regression
NLL <- function(beta, y, Xmat) {
  lambda <- exp(Xmat %*% beta) # Multiply design matrix by beta
  LL <- dpois(y, lambda, log = TRUE) # Log-likelihood of each obs
  NLL <- -sum(LL) # NLL for all observations in data set
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
# Need to use method = "BFGS" here
inits <- rep(0, 6)
names(inits) <- colnames(Xmat)
out13.7 <- optim(inits, NLL, y = load, Xmat = Xmat, hessian = TRUE, method = "BFGS")
get_MLE(out13.7, 4)
diy_est <- out13.7$par


# 13.8 Likelihood analysis with Template Model Builder
# ----------------------------------------------------

# Bundle and summarize data
str(tmbData <- list(load = load, pop = as.numeric(pop) - 1,
nPops = nPops, wing_length = wing.length, n = n) ) # not shown

# Write TMB model file
cat(file = "model13_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type> ::operator() ()
{
  //Describe input data
  DATA_VECTOR(load); //response
  DATA_VECTOR(wing_length); //covariate
  DATA_IVECTOR(pop); //population index (of integers; IVECTOR)
  DATA_INTEGER(n); //Number of observations

  //Describe parameters
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type LL = 0; //Initialize log-likelihood at 0
  for (int i = 0; i < n; i ++ ){
    Type lambda = exp(alpha(pop(i)) + beta(pop(i)) * wing_length(i));
    //Calculate log-likelihood of observation and add to total
    LL += dpois(load(i), lambda, true);
  }

  //Derived parameters
  Type a_effe2 = alpha(1) - alpha(0);
  Type a_effe3 = alpha(2) - alpha(0);
  Type b_effe2 = beta(1) - beta(0);
  Type b_effe3 = beta(2) - beta(0);
  Type test1 = beta(2) - beta(1);
  ADREPORT(a_effe2);
  ADREPORT(a_effe3);
  ADREPORT(b_effe2);
  ADREPORT(b_effe3);
  ADREPORT(test1);
  return -LL; //Return negative log likelihood
}
")

# Compile and load TMB function
compile("model13_8.cpp")
dyn.load(dynlib("model13_8"))

# Provide dimensions and starting values for parameters
params <- list(alpha = rep(0, tmbData$nPops), beta = rep(0, tmbData$nPops))

# Create TMB object
out13.8 <- MakeADFun(data = tmbData,
                     parameters = params,
                     DLL = "model13_8", silent = TRUE)

# Optimize TMB object and print results
opt <- optim(out13.8$par, fn = out13.8$fn, gr = out13.8$gr, method = "BFGS")
(tsum <- tmb_summary(out13.8))
tmb_est <- tsum[c(1, 7, 8, 4, 9, 10), 1] # Save parameter estimates


# 13.9 Comparison of the parameter estimates
# ------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = beta.vec, glm = glm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)

# 13.10 Summary
# -------------
# (no code)
