
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -------------------------------------------------------------------
# Chapter 15  --  Comparing two groups in a logistic regression model
# -------------------------------------------------------------------

# Last changes: 11 June 2024


# 15.1 Introduction
# -----------------
# (no code)


# 15.2 Data generation
# --------------------

# Simulate Bernoulli variant of the data set first
set.seed(15)
N <- 50 # Total number of sites (binomial total)
theta.cr <- 12/50 # Success probability Cross-leaved
theta.ch <- 38/50 # Success probability Chiltern gentian
# Simulate 50 'coin flips' for each species
y.cr <- rbinom(N, 1, prob = theta.cr); y.cr # Det/nondet data for cr
y.ch <- rbinom(N, 1, prob = theta.ch); y.ch # ditto for ch
y <- c(y.cr, y.ch) # Merge the two binary vectors
species.long <- factor(rep(c(0,1), each = N), labels = c("Cross-leaved", "Chiltern"))
data.frame('species' = species.long, 'det.nondet' = y) # not shown
# Aggregate the binary data to become two binomial counts
C <- c(sum(y.cr), sum(y.ch)) # Tally up detections
species <- factor(c(0,1), labels = c("Cross-leaved", "Chiltern"))
# Save true parameter values
truth <- c(Intercept = log(theta.cr/(1-theta.cr)),
chiltern = (log(theta.ch/(1-theta.ch)) - log(theta.cr/(1-theta.cr))))

# Required packages
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 15.3 Likelihood analysis with canned functions in R
# ---------------------------------------------------

# Fit the Bernoulli model to the binary variant of the data set
out15.3X <- glm(y ~ species.long, family = binomial)
summary(out15.3X)
predict(out15.3X, type = "response")[c(1, 51)] # Pred for two species

# Fit the Binomial model to the aggregated counts and save estimates
out15.3 <- glm(cbind(C, N - C) ~ species, family = binomial)
summary(out15.3)
predict(out15.3, type = "response")
glm_est <- out15.3$coef

# 15.4 Bayesian analysis with JAGS

# Bundle and summarize data
str(dataList <- list(C = C, species = c(0,1), N = 50, n = 2) )

# Write JAGS model file
cat(file = "model15.4.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:n) {
  C[i] ~ dbin(theta[i], N) # Note theta before N in JAGS
  logit(theta[i]) <- alpha + beta * species[i] # species Chiltern = 1
  # theta[i] <- ilogit(alpha + beta * species[i]) # same (see also STAN)
}

# Derived quantities
Occ.cross <- exp(alpha)/(1 + exp(alpha))
Occ.chiltern <- exp(alpha + beta)/(1 + exp(alpha + beta))
Occ.Diff <- Occ.chiltern - Occ.cross # Test quantity
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1))}
# Parameters to estimate
params <- c("alpha", "beta", "Occ.cross", "Occ.chiltern", "Occ.Diff")
# MCMC settings
na <- 1000 ; ni <- 6000 ; nb <- 1000 ; nc <- 4 ; nt <- 5

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save estimates
out15.4 <- jags(dataList, inits, params, "model15.4.txt",
n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out15.4) # not shown
print(out15.4, 3)
jags_est <- unlist(out15.4$mean[1:2])

# Draw Fig. 15-2
par(mfrow = c(1, 2), mar = c(5, 5, 5, 3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(out15.4$sims.list$Occ.cross, main = "", col = "red", xlab = "Observed occupancy",
  xlim = c(0,1), breaks = 20, freq = FALSE, ylim = c(0, 8))
hist(out15.4$sims.list$Occ.chiltern, col = "blue", breaks = 20, freq = FALSE, add = TRUE)
legend('topright', legend = c('Cross-leaved', 'Chiltern'),
      fill = c("red", "blue"), cex = 1.5)

hist(out15.4$sims.list$Occ.Diff, col = "grey", las = 1,
xlab = "Difference in occupancy", main = "", breaks = 30,
freq = FALSE, xlim = c(-0.2, 0.6), ylim = c(0, 5))
abline(v = 0, col = 'grey', lty = 2, lwd = 3)

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est)
print(comp, 4)


# 15.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(C = C, species = c(0, 1), N = 50, n = 2) )


# Write NIMBLE model file
model15.5 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, sd = 100)
beta ~ dnorm(0, sd = 100)

# Likelihood
for (i in 1:2) {
  C[i] ~ dbin(p[i], N)       # Note again p before N
  logit(p[i]) <- alpha + beta *species[i]
}

# Derived quantities
Occ.cross <- exp(alpha) / (1 + exp(alpha))
Occ.chiltern <- exp(alpha + beta) / (1 + exp(alpha + beta))
Occ.Diff <- Occ.chiltern - Occ.cross     # Test quantity
} )


# Inits
inits <- function(){ list(alpha = rlnorm(1), beta = rlnorm(1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "Occ.cross", "Occ.chiltern", "Occ.Diff")

# MCMC settings
na <- 1000  ;  ni <- 6000  ;  nb <- 1000  ; nc <- 4  ; nt <- 5

# Call NIMBLE (ART 20 sec), check convergence, summarize posteriors and save estimates
system.time( out15.5 <- 
    nimbleMCMC(code = model15.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out15.5)
(nsum <- nimble_summary(out15.5, params))
nimble_est <- nsum[1:2,1]     # Save param estimates


# 15.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, species = c(0,1), N = 50, n = 2) )

# Write Stan model
cat(file = "model15_6.stan", "
data{
  array[2] int C;
  array[2] int species;
  int N;
  int n;
}

parameters{
  real alpha;
  real beta;
}

model{
  vector[2] theta;

  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  
  for (i in 1:n){
    theta[i] = inv_logit(alpha + beta * species[i]);
    C[i] ~ binomial(N, theta[i]);
  }
}

generated quantities{
  real Occ_cross = inv_logit(alpha);
  real Occ_chiltern = inv_logit(alpha + beta);
  real Occ_Diff = Occ_chiltern - Occ_cross;
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 42/1 sec), assess convergence, print and save results
system.time(
out15.6 <- stan(file = "model15_6.stan", data = dataList,
warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out15.6) # not shown
print(out15.6, dig = 3) # not shown
stan_est <- summary(out15.6)$summary[1:2,1] # save estimates


# 15.7 Do-it-yourself MLEs
# ------------------------

# Definition of NLL for logistic regression with 1 factor
NLL <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  LL <- theta <- numeric(2)
  for (i in 1:2){
    theta[i] <- plogis(alpha + beta * data$species[i])
    LL[i] <- dbinom(data$C[i], data$N, theta[i], log = TRUE)
  }
  NLL <- -sum(LL) # NLL for all observations
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- c('Intercept' = 0, 'beta.Chiltern' = 0)
out15.7 <- optim(inits, NLL, data = dataList, hessian = TRUE)
get_MLE(out15.7, 4)
diy_est <- out15.7$par # Save estimates


# 15.8 Likelihood analysis with TMB
# ---------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, species = c(0, 1), N = 50, n = 2) )

# Write TMB model file
cat(file = "model15_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(C);
  DATA_VECTOR(species);
  DATA_INTEGER(N);
  DATA_INTEGER(n);
  
  //Describe parameters
  PARAMETER(alpha);
  PARAMETER(beta);
  
  Type LL = 0.0; //Initialize log-likelihood at 0
  
  for (int i = 0; i < n; i++){
    Type theta = invlogit(alpha + beta * species(i));
    LL += dbinom(C(i), Type(N), theta, true);
  }

  Type Occ_cross = invlogit(alpha);
  Type Occ_chiltern = invlogit(alpha + beta);
  Type Occ_Diff = Occ_chiltern - Occ_cross;
  ADREPORT(Occ_cross);
  ADREPORT(Occ_chiltern);
  ADREPORT(Occ_Diff);
  return -LL;
}
")

# Compile and load TMB function
compile("model15_8.cpp")
dyn.load(dynlib("model15_8"))

# Provide dimensions and starting values for parameters
params <- list(alpha = 0, beta = 0)

# Create TMB object
out15.8 <- MakeADFun(data = dataList,
                     parameters = params,
                     DLL = "model15_8", silent = TRUE)

# Optimize TMB object and print results
starts <- rep(0, 2)
opt <- optim(starts, fn = out15.8$fn, gr = out15.8$gr, method = "BFGS")
(tsum <- tmb_summary(out15.8)) # not shown
tmb_est <- tsum[1:2,1] # save estimates


# 15.9 Comparison of the parameter estimates
# ------------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, glm = glm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)
