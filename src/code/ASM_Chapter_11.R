
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ------------------------------------------------------
# Chapter 11  --  Introduction to the generalized linear
#                 model (GLM): comparing two groups in a
#                 Poisson regression
# ------------------------------------------------------

# Last changes: 11 June 2024


# 11.1 Introduction
# -----------------
# (no code)


# 11.2 An important but often forgotten issue with count data
# -----------------------------------------------------------
# (no code)


# 11.3 How to deal with missing values in our data
# ------------------------------------------------
# (no code)


# 11.4 Data generation
# --------------------

set.seed(11)
nSites <- 30
x <- gl(n = 2, k = nSites, labels = c("grassland", "arable"))
n <- 2 * nSites

lambda <- exp(0.69 + 0.92*(as.numeric(x)-1))           # x has levels 1 and 2, not 0 and 1

# Save true values of parameters for later comparison
truth <- c(intercept = 0.69, arable = 0.92)

y <- rpois(n = n, lambda = lambda)                     # Add Poisson noise
y[c(1, 10, 35)] <- NA                                  # Make some observations NA
y                                                      # Print the counts (not shown)

aggregate(y, by = list(x), FUN = mean, na.rm = TRUE)   # Observed means
boxplot(y ~ x, col = "grey", xlab = "Land-use",
  ylab = "Hare count", las = 1, frame = FALSE)         # Fig. 11–2  

# Required packages
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)


# 11.5 Likelihood analysis with canned functions in R
# ---------------------------------------------------

out11.5 <- glm(y ~ x, family = poisson)                # Fit the model
summary(out11.5)                                       # Two-group comparison
anova(out11.5, test = "Chisq")                         # Likelihood ratio test (LRT)

# Save estimates
glm_est <- coef(out11.5)

# Check goodness-of-fit using traditional and quantile residuals
plot(out11.5)                                          # Traditional residual check for comparison
simOut <- simulateResiduals(out11.5, n = 1000, plot = TRUE)


# 11.6 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# Write JAGS model file
cat(file = "model11.6.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:n) {
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha + beta *x[i]
  # lambda[i]) <- exp(alpha + beta *x[i])              # Same
}

# Fit assessments
for (i in 1:n) {
  Presi[i] <- (y[i]-lambda[i]) / sqrt(lambda[i])       # Pearson residuals
  y.new[i] ~ dpois(lambda[i])                          # Replicate data set
  Presi.new[i] <- (y.new[i]-lambda[i]) / sqrt(lambda[i]) # Pearson resi
  D[i] <- pow(Presi[i], 2)
  D.new[i] <- pow(Presi.new[i], 2)
}
# Add up discrepancy measures
fit <- sum(D[])
fit.new <- sum(D.new[])
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit.new", "y")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out11.6 <- jags(dataList, inits, params, "model11.6.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out11.6) # not shown (yet)
print(out11.6, 3) # not shown


#   11.6.1 Assessment of model adequacy
# -------------------------------------

# Fig. 11-3
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(out11.6$mean$Presi, ylab = "Pearson residual", las = 1, pch = 16,
  cex = 1.5, col = rgb(0,0,0,0.5), frame = FALSE)
abline(h = 0)
plot(out11.6$sims.list$fit, out11.6$sims.list$fit.new, main = "",
  xlab = "Discrepancy actual data set", ylab = "Discrepancy perfect data sets",
  pch = 16, cex = 1.5, col = rgb(0,0,0,0.2), frame = FALSE)
abline(0, 1, lwd = 2, col = "black")

# Bayesian p-value
mean(out11.6$sims.list$fit.new > out11.6$sims.list$fit)

# Fig. 11.4
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(out11.6$sims.list$beta, col = "grey", las = 1,
  xlab = "Coefficient", main = "Effect of 'arable'",
  breaks = 50, freq = FALSE, xlim= c(0,2.2))
abline(v=0, col='red', lty=2, lwd = 3)
hist(exp(out11.6$sims.list$alpha), main= "Comparison of land covers",
  col = "red", xlab = "Expected hare count", xlim = c(0,7),
  breaks = 20, freq = FALSE)
hist(exp(out11.6$sims.list$alpha + out11.6$sims.list$beta),
  col = "blue", breaks = 20, freq = FALSE, add = TRUE)
  legend('topright', legend=c('Grassland', 'Arable'), 
  fill = c("red", "blue"), cex = 1.5)


#   11.6.2 Inference under the model
# ----------------------------------

print(out11.6, 2)

# Compare likelihood with Bayesian estimates and with truth
jags_est <- unlist(out11.6$mean[1:2])
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est)
print(comp, 4)


# 11.7 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data (same as for JAGS)
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# Write NIMBLE model file
model11.7 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, sd = 100)
beta ~ dnorm(0, sd = 100)

# Likelihood
for (i in 1:n) {
  y[i] ~ dpois(lambda[i]) 
  log(lambda[i]) <- alpha + beta * x[i]
}

# Fit assessments
for (i in 1:n) {
  Presi[i] <- (y[i]-lambda[i]) / sqrt(lambda[i])   # Pearson residuals
  y.new[i] ~ dpois(lambda[i])                      # Replicate data set
  Presi.new[i] <- (y.new[i]-lambda[i]) / sqrt(lambda[i]) # Pearson resi
  D[i] <- pow(Presi[i], 2)
  D.new[i] <- pow(Presi.new[i], 2)
}

# Add up discrepancy measures
 fit <- sum(D[1:20])
 fit.new <- sum(D.new[1:20])
} )

# Inits
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit.new")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 28 sec), check convergence and summarize posteriors
system.time( out11.7 <- 
    nimbleMCMC(code = model11.7,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out11.7)
(nsum <- nimble_summary(out11.7, params))

# Save estimates from NIMBLE
nimble_est <- nsum[1:2,1]


# 11.8 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
not_na <- which(!is.na(y))
x_stan <- as.numeric(x) - 1
x_stan <- x_stan[not_na]
str(dataList <- list(y = y[not_na], x = x_stan, n = length(x_stan)))

# Write Stan model
cat(file = "model11_8.stan", "
data {
  int n;
  array[n] int y;
  vector[n] x;
}

parameters {
  real alpha;
  real beta;
}

transformed parameters {
  vector[n] lambda;
  for (i in 1:n) {
    lambda[i] = exp(alpha + beta * x[i]);
  }
}

model {

  //Priors
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  
  //Likelihood
  for (i in 1:n) {
    y[i] ~ poisson(lambda[i]);
  }
}

generated quantities {
  vector[n] Presi;
  vector[n] D;
  array[n] int y_new;
  vector[n] Presi_new;
  vector[n] D_new;
  real fit;
  real fit_new;
  for (i in 1:n) {
    Presi[i] = (y[i] - lambda[i])/sqrt(lambda[i]);
    D[i] = pow(Presi[i], 2);
    y_new[i] = poisson_rng(lambda[i]);
    Presi_new[i] = (y_new[i] - lambda[i])/sqrt(lambda[i]);
    D_new[i] = pow(Presi_new[i], 2);
  }
  fit = sum(D);
  fit_new = sum(D_new);
}
")

# Parameters monitored: same as before
params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit_new")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 32/3 sec), assess convergence, print estimates and save point estimates
system.time(
  out11.8 <- stan(file = "model11_8.stan", data = dataList, pars = params,
                  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out11.8) # not shown
print(out11.8, dig = 3) # not shown
stan_est <- summary(out11.8)$summary[1:2,1]


# 11.9 Do-it-yourself MLEs
# ------------------------

# Definition of NLL for a Poisson GLM with 1 covariate
NLL <- function(param, y, x) {
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x) # Calculate lambda for each datum
  LL <- dpois(y, lambda, log = TRUE) # Log-likelihood for each datum
  -sum(LL, na.rm = TRUE) # NLL for all observations
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('Intercept' = 0, 'beta(arable)' = 0)
out11.9 <- optim(inits, NLL, y = y, x = as.numeric(x)-1, hessian = TRUE)
get_MLE(out11.9, 5)
diy_est <- out11.9$par


# 11.10 Likelihood analysis with TMB
# ----------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, x = as.numeric(x)-1, n = length(x)))

# Write TMB model file
cat(file = "model11_10.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(y); //response
  DATA_VECTOR(x); //covariate
  DATA_INTEGER(n); //Number of observations

  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  
  Type LL = 0.0; //Initialize log-likelihood at 0
  Type lambda; //Initialize lambda

  for (int i= 0; i<n; i++){
    if(R_IsNA(asDouble(y(i)))) continue; // Here accommodate NAs
    lambda = exp(alpha + beta * x(i)); //Expected count

    //Calculate log-likelihood of observation and add to total
    LL += dpois(y(i), lambda, true);
  }

  return -LL; //Return negative log likelihood
}
")


# Compile and load TMB function
compile("model11_10.cpp")
dyn.load(dynlib("model11_10"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0)

# Create TMB object
out11.10 <- MakeADFun(data = dataList,
                      parameters = params,
                      DLL = "model11_10", silent = TRUE)

# Optimize TMB object and print results
opt <- optim(out11.10$par, fn = out11.10$fn, gr = out11.10$gr, method = "BFGS")
(tsum <- tmb_summary(out11.10)) # not shown
tmb_est <- c(opt$par[1:2]) # save estimates


# 11.11 Comparison of the parameter estimates
# -------------------------------------------

# Compare estimates with truth
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 11.12 Summary and outlook
# -------------------------
# (no code)
