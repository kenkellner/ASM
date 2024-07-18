
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ---------------------------------
# Chapter 19  --  Occupancy models
# ---------------------------------

# Last changes: 11 June 2024


# 19.1 Introduction
# -----------------
# (no code)


# 19.2 Data generation
# --------------------

set.seed(19)
nSites <- 150 # 150 sites visited

humidity <- runif(n = nSites, min = -1, max = 1)

alpha.occ <- 0 # Logit-linear intercept for humidity on occupancy
beta.occ <- 2 # Logit-linear slope for humidity on occupancy
occ.prob <- plogis(alpha.occ + beta.occ * humidity)

# Fig. 19–2 left
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(humidity, occ.prob, ylim = c(0,1), xlab = "Humidity index", ylab = "Occupancy probability",
main = "State process", las = 1, pch = 16, cex = 2, col = rgb(0,0,0,0.3), frame = FALSE)
z <- rbinom(n = nSites, size = 1, prob = occ.prob)
z # Look at the true occupancy state of each site
(true_Nz <- sum(z)) # Number of occupied sites among the visited sites

alpha.p <- 0 # Logit-linear intercept for humidity on detection
beta.p <- -3 # Logit-linear slope for humidity on detection
lp <- alpha.p + beta.p * humidity # The linear predictor for detection
p <- plogis(lp) # Get p on the probability scale

# Save true parameter values into a vector
truth <- c(alpha.occ = alpha.occ, beta.occ = beta.occ,
alpha.p = alpha.p, beta.p = beta.p)

# Fig. 19–2 right
plot(humidity, p, ylim = c(0,1), xlab = "Humidity index", ylab = "Detection probability",
     main = "Observation process", las = 1, pch = 16, cex = 2, col = rgb(0,0,0,0.3), frame = FALSE)

eff.p <- z * p # Zeroes out detection at sites where gentian is absent

nVisits <- 3
y <- array(dim = c(nSites, nVisits))

# Simulate results of first through last surveys
for(t in 1:nVisits){
  y[,t] <- rbinom(n = nSites, size = 1, prob = eff.p)
}

y # Look at the data
# Apparent number of occupied sites among the 150 visited sites
(obs_Nz <- sum(apply(y, 1, sum) > 0)) # 39
cbind('True state' = z, 'Obs Visit 1' = y[,1],
  'Obs Visit 2' = y[,2], 'Obs Visit 3' = y[,3])
# Look true state and at three
# observed states at each site

obsZ <- as.numeric(apply(y, 1, sum) > 0) # 'Observed presence/absence'
naive.analysis <- glm(obsZ ~ humidity, family = binomial)
summary(naive.analysis)
lpred.naive <- predict(naive.analysis, type = 'link', se = TRUE)
pred.naive <- plogis(lpred.naive$fit)
LCL.naive <- plogis(lpred.naive$fit-2*lpred.naive$se)
UCL.naive <- plogis(lpred.naive$fit + 2*lpred.naive$se)

# Fig. 19–3
par(mfrow = c(1, 1), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(humidity, pred.naive, ylim = c(0, 0.6), xlab = "Humidity index",
  ylab = "Apparent occupancy prob.", main = "Confounding of state and observation processes",
  las = 1, pch = 16, cex = 2, col = rgb(0,0,0,0.4), frame = FALSE)
polygon(c(sort(humidity), rev(humidity[order(humidity)])), c(LCL.naive[order(humidity)],
  rev(UCL.naive[order(humidity)])), col = rgb(0,0,0, 0.2), border = NA)

library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 19.3 Likelihood analysis with canned functions in the R package unmarked
# ------------------------------------------------------------------------

# Load unmarked, format data into unmarked data frame and summarize
library(unmarked)
summary( umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(humidity = humidity)) )

# Fit model and extract estimates
# Detection covariates follow first tilde, then occupancy covariates
summary(out19.3 <- occu(~humidity ~humidity, data = umf))
unm_est <- coef(out19.3) # Save estimates

# Get estimates of latent occurrence from unmarked output with ranef()
unm_Nz <- round(sum(ranef(out19.3)@post[,2,1]), 2)
comp <- c(truth = true_Nz, observed = obs_Nz, unmarked = unm_Nz)
print(comp, 2) # not shown

# Make predictions
state.pred <- predict(out19.3, type = 'state')
det.pred <- predict(out19.3, type = 'det')
p.pred <- matrix(det.pred[,1], nrow = nSites, byrow = TRUE) # reformat
p.LCL <- matrix(det.pred[,3], nrow = nSites, byrow = TRUE) # reformat
p.UCL <- matrix(det.pred[,4], nrow = nSites, byrow = TRUE) # reformat

# Predict humidity relation for state and observation processes
# (Fig. 19.4)
ooo <- order(humidity) # Get the order of the humidity values
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)

# Occupancy
plot(humidity[ooo], state.pred[ooo,1], xlab = 'Humidity index', ylab = 'Occupancy probability',
  frame = FALSE, col = 'blue', lwd = 3, main = 'State process', type = 'l', ylim = c(0, 1))
lines(humidity[ooo], occ.prob[ooo], lwd = 3, col = 'red')
polygon(c(humidity[ooo], rev(humidity[ooo])), c(state.pred[ooo,3], rev(state.pred[ooo,4])),
  col = rgb(0,0,1, 0.1), border = NA)
legend('bottomright', legend = c('Estimate', 'Truth'), lwd = 2, col = c('blue', 'red'),
  bty = 'n', cex = 1.2)

# Detection
plot(humidity[ooo], p.pred[ooo,1], xlab = 'Humidity index', ylab = 'Detection probability',
  frame = FALSE, col = 'blue', lwd = 3, main = 'Observation process', type = 'l', ylim = c(0, 1))
lines(humidity[ooo], p[ooo], lwd = 3, col = 'red')
  polygon(c(humidity[ooo], rev(humidity[ooo])), c(p.LCL[ooo,1], rev(p.UCL[ooo,1])),
  col = rgb(0,0,1, 0.1), border = NA)


# 19.4 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, humidity = humidity, nSites = nSites, nVisits = nVisits) )

# Write JAGS model file
cat(file = "model19.4.txt", "
model {
# Priors
occ.int ~ dunif(0, 1) # Occupancy intercept on prob. scale
alpha.occ <- logit(occ.int)
beta.occ ~ dnorm(0, 0.0001)
p.int ~ dunif(0, 1) # Detection intercept on prob. scale
alpha.p <- logit(p.int)
beta.p ~ dnorm(0, 0.0001)

# Likelihood
for (i in 1:nSites) { # start initial loop over sites
# True state model for the partially observed true state
  z[i] ~ dbern(psi[i]) # True occupancy z at site i
  logit(psi[i]) <- alpha.occ + beta.occ * humidity[i]
  for (t in 1:nVisits) { # start a second loop over visits
    # Observation model for the actual observations
    y[i,t] ~ dbern(eff.p[i,t]) # Detection-nondetection at i and t
    eff.p[i,t] <- z[i] * p[i,t]
    logit(p[i,t]) <- alpha.p + beta.p * humidity[i]
  }
}

# Derived quantities
# Finite-sample inference on number occ. sites in our sample of 150
occ.fs <- sum(z[])
}
")

# Function to generate starting values
zst <- apply(y, 1, max) # Get *observed* presence/absence status
inits <- function(){list(z = zst, occ.int = runif(1), beta.occ = runif(1, -3, 3),
  p.int = runif(1), beta.p = runif(1, -3, 3))}

# Parameters to estimate
params <- c("alpha.occ","beta.occ", "alpha.p", "beta.p", "occ.fs", "occ.int", "p.int")

# MCMC settings
na <- 5000 ; ni <- 50000 ; nb <- 10000 ; nc <- 4 ; nt <- 10

# Call JAGS (ART 60 sec), check convergence, summarize posteriors and save results
out19.4 <- jags(dataList, inits, params, "model19.4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out19.4) # not shown
print(out19.4, 3)
jags_est <- out19.4$summary[1:4,1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, unmarked = unm_est, JAGS = jags_est)
print(comp, 4)

# Back-transform the intercept estimates
plogis(unm_est[c(1,3)]) # Estimates from unmarked (MLEs)
plogis(jags_est[c(1,3)]) # Estimates from JAGS (posterior means)

# Get estimate of the number of occupied sites
jags_Nz <- round(out19.4$mean$occ.fs,2)
comp <- c(truth = true_Nz, observed = obs_Nz, unmarked = unm_Nz, JAGS = jags_Nz)
print(comp, 2)

# Get predictions
str(post.draws <- out19.4$sims.list) # Grab posterior draws and look
nsamp <- length(post.draws[[1]]) # Check number of posterior draws
pred.occ <- array(NA, dim = c(length(humidity), nsamp))
for(i in 1:length(humidity)){ # Posterior predictive distribution
  pred.occ[i,] <- plogis(post.draws$alpha.occ + post.draws$beta.occ * humidity[i])
}
pm <- apply(pred.occ, 1, mean) # Posterior mean
CRI <- apply(pred.occ, 1,
function(x) quantile(x, c(0.025, 0.975))) # Central 95% percentiles

# Fig. 19–5
par(mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(humidity, pred.naive, ylim = c(0, 1), xlab = "Humidity index",
  ylab = "Occupancy probability", main = "", las = 1, pch = 16, cex = 1.6,
col = rgb(0, 0, 0, 0.4), frame = FALSE) # Detection-naive analysis
polygon(c(humidity[ooo], rev(humidity[ooo])), c(LCL.naive[ooo], rev(UCL.naive[ooo])),
  col = rgb(0, 0, 0, 0.2), border = NA)
points(humidity, occ.prob, pch = 16, cex = 1.6, col = rgb(1, 0, 0, 0.4)) # Truth
points(humidity, pm, pch = 16, cex = 1.6, col = rgb(0, 0, 1, 0.4)) # Esti
polygon(c(humidity[ooo], rev(humidity[ooo])), c(CRI[1,ooo], rev(CRI[2,ooo])),
  col = rgb(0, 0, 1, 0.2), border = NA)
legend('topleft', legend = c('Site-occ SDM', 'Truth', 'Detection-naive SDM'),
pch = 16, col = c('blue', 'red', 'black'), bty = 'n', cex = 1.4)


# 19.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data (same as for JAGS)
str(dataList <- list(y = y, humidity = humidity, nSites = nSites, nVisits = nVisits) )

# Write NIMBLE model file
model19.5 <- nimbleCode( {
# Priors
occ.int ~ dunif(0, 1)
alpha.occ <- logit(occ.int)
beta.occ ~ dnorm(0, sd = 100)
p.int ~ dunif(0, 1)
alpha.p <- logit(p.int)
beta.p ~ dnorm(0, , sd = 100)

# Likelihood
for (i in 1:nSites) {           # start initial loop over sites
# True state model for the partially observed true state
  z[i] ~ dbern(psi[i])          # True occupancy z at site i
  logit(psi[i]) <- alpha.occ + beta.occ * humidity[i]

  for (t in 1:nVisits) {        # start a second loop over visits
    # Observation model for the actual observations
    y[i,t] ~ dbern(eff.p[i,t])  # Detection-nondetection at i and t
    eff.p[i,t] <- z[i] * p[i,t]
    logit(p[i,t]) <- alpha.p + beta.p * humidity[i]
  }
}

# Derived quantities
occ.fs <- sum(z[1:nSites])      # Number of occupied sites among 150
} )

# Inits
zst <- apply(y, 1, max)  # Good starting values essential !
inits <- function(){list(z = zst, occ.int=runif(1), beta.occ = runif(1, -3, 3), p.int = runif(1), beta.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.occ","beta.occ", "alpha.p", "beta.p", "occ.fs", "occ.int", "p.int")

# MCMC settings
ni <- 50000  ;  nb <- 10000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 95 sec), check convergence, summarize posteriors and save results
system.time( out19.5 <- 
    nimbleMCMC(code = model19.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out19.5)   # not shown
(nsum <- nimble_summary(out19.5, params))     # not shown
nimble_est <- nsum[1:4,1]    # Save estimates
nimble_Nz <- round(nsum[5,1],2)


# 19.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(y = y, humidity = humidity, nSites = nSites, nVisits = nVisits) )

# Add new variable indicating if there were 0 detections at a site
dataList$no_detect <- 1 - apply(dataList$y, 1, max)

# Add new variable for total # of detections at each site
dataList$nd <- apply(dataList$y, 1, sum)
str(dataList)

# Write Stan model
cat(file = "model19_6.stan", "
data{
  int nSites;
  int nVisits;
  array[nSites, nVisits] int y;
  vector[nSites] humidity;
  array[nSites] int no_detect;
  array[nSites] int nd;
}

parameters{
  real occ_int;
  real beta_occ;
  real p_int;
  real beta_p;
}

transformed parameters{
  vector[nSites] psi;
  vector[nSites] p;
  real alpha_occ = logit(occ_int);
  real alpha_p = logit(p_int);
  
  for (i in 1:nSites){
    psi[i] = inv_logit(alpha_occ + beta_occ * humidity[i]);
    p[i] = inv_logit(alpha_p + beta_p * humidity[i]);
  }
}

model{
  vector[nSites] lik; //Likelihood for each site
  
  //priors
  occ_int ~ uniform(0, 1);
  beta_occ ~ normal(0, 100);
  p_int ~ uniform(0, 1);
  beta_p ~ normal(0, 100);
  
  for (i in 1:nSites){
    //Calculate occupancy likelihood
    lik[i] = exp(binomial_lpmf(nd[i] | nVisits, p[i])) * psi[i] +
    no_detect[i] * (1-psi[i]);
    target += log(lik[i]); //Add to total likelihood
  }
}

generated quantities{
  real qT;
  vector[nSites] psi_con; //P(z = 1 | Number of detects = 0)
  array[nSites] int zpost;
  real occ_fs; //Number of occupied sites

  for (i in 1:nSites){
    if(no_detect[i] == 0){
      zpost[i] = 1; //If species was detected at least once, z = 1
    } else {
      qT = pow(1-p[i], nVisits); //P(never detected | z = 1)
      psi_con[i] = psi[i] * qT/(psi[i] * qT + (1-psi[i]));
      zpost[i] = bernoulli_rng(psi_con[i]);
  }
  }
  occ_fs = sum(zpost);
}
")

# Inits
inits <- function(){ list(occ_int = runif(1,0,1), beta_occ = runif(1, -5, 5),
  p_int = runif(1,0,1), beta_p = runif(1, -5, 5))}

# Parameters to estimate
params <- c("alpha_occ", "beta_occ", "alpha_p", "beta_p", "occ_fs")

# HMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 75/26 sec), assess convergence and print results table
system.time(
out19.6 <- stan(file = "model19_6.stan", data = dataList, pars = params,
                warmup = nb, iter = ni, chains = nc, thin = nt,
                init = inits, control = list(adapt_delta = 0.99)) )
rstan::traceplot(out19.6) # not shown
print(out19.6, dig = 3) # not shown
stan_est <- summary(out19.6)$summary[1:4,1] # Save estimates
stan_Nz <- round(summary(out19.6)$summary[5,1],2)


# 19.7 Bayesian analysis with canned functions in the R package ubms
# ------------------------------------------------------------------

library(ubms) # will also load unmarked as a dependency

# Format data and summarize data (exactly the same as in Section 19.3)
summary( umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(humidity = humidity)) )

# Fit model in parallel
options(mc.cores = 4) # number of parallel cores to use
system.time( # ART 28 sec
  out19.7 <- stan_occu(~ humidity ~ humidity, umf, chains = 4,
                       prior_coef_state = normal(0, 10),
                       prior_coef_det = normal(0, 10))
)
out19.7
ubms_est <- coef(out19.7) # Save estimates

# Get estimate of the number of occupied sites and save
postz_ubms <- posterior_predict(out19.7, "z") # Get posterior of latent z
post_Nz_ubms <- apply(postz_ubms, 1, sum) # Calculate sumN for each draw
ubms_Nz <- mean(post_Nz_ubms)


# 19.8 Do-it-yourself MLEs
# ------------------------

# Same data as section 19.6
str(dataList) # not shown

# Definition of NLL for a simple static site-occupancy model
# with one single site covariate humidity and fit to detection frequency data
NLL <- function(param, data) {
  alpha.lpsi <- param[1] # Occupancy intercept (logit scale)
  beta.lpsi <- param[2] # Occupancy slope on humidity
  alpha.lp <- param[3] # Detection intercept (logit scale)
  beta.lp <- param[4] # Detection slope on humidity
  
  # Linear predictor for occupancy
  psi <- plogis(alpha.lpsi + beta.lpsi * data$humidity)
  # Linear predictor for detection
  p <- plogis(alpha.lp + beta.lp * data$humidity)
  
  L <- numeric(data$nSites)
  
  for (i in 1:data$nSites){
    # Likelihood contribution for 1 observation
    L[i] <- psi[i] * dbinom(data$nd[i], data$nVisits, p[i]) +
              data$no_detect[i] * (1-psi[i])
  }
  
  LL <- log(L) # Log-likelihood contr. for each observation
  NLL <- -sum(LL) # NLL for all observations
  return(NLL)
}

# Minimize that NLL to find MLEs and also get SEs
inits <- rep(0, 4)
names(inits) <- names(truth)
out19.8 <- optim(inits, NLL, data = dataList, hessian = TRUE, method = "BFGS")
get_MLE(out19.8, 5)
diy_est <- out19.8$par


# 19.9 Likelihood analysis with TMB
# ---------------------------------

# Same data as section 19.6
str(dataList) # not shown

# Write TMB model file
cat(file = "model19_9.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_VECTOR(nd);
  DATA_VECTOR(no_detect);
  DATA_VECTOR(humidity);
  DATA_INTEGER(nSites);
  DATA_INTEGER(nVisits);
  
  //Describe parameters
  PARAMETER(alpha_occ);
  PARAMETER(beta_occ);
  PARAMETER(alpha_p);
  PARAMETER(beta_p);

  Type LL = 0.0; //Initialize log-likelihood at 0

  vector <Type> p(nSites);
  vector <Type> psi(nSites);
  vector <Type> L(nSites);

  for (int i = 0; i < nSites; i++){
    psi(i) = invlogit(alpha_occ + beta_occ * humidity(i));
    p(i) = invlogit(alpha_p + beta_p * humidity(i));
  
    L(i) = dbinom(nd(i), Type(nVisits), p(i), false) * psi(i) +
            no_detect(i) * (1 - psi[i]);
    LL += log(L(i));
  
  }
  return -LL;
}
")

# Compile and load TMB function
compile("model19_9.cpp")
dyn.load(dynlib("model19_9"))

# Provide dimensions and starting values for parameters
params <- list(alpha_occ = 0, beta_occ = 0, alpha_p = 0, beta_p = 0)

# Create TMB object
out19.9 <- MakeADFun(data = dataList, parameters = params,
                     DLL = "model19_9", silent = TRUE)

# Optimize TMB object and print results
starts <- rep(0, 4)
opt <- optim(starts, fn = out19.9$fn, gr = out19.9$gr, method = "BFGS")
(tsum <- tmb_summary(out19.9))
tmb_est <- tsum[1:4,1] # Save estimates

# Estimate number of occupied sites with a bootstrap

# Generate samples of parameter vector
nsims <- 10000
library(MASS)# for mvrnorm()
sdr <- sdreport(out19.9)
Beta <- sdr$par.fixed
Sigma <- sdr$cov.fixed
param_samples <- mvrnorm(nsims, Beta, Sigma)

# Calculate estimated sum of sites occupied for each parameter vector sample
# conditional on observed data
occ_fs_sims <- rep(NA, nsims)
for (i in 1:nsims){
  beta <- param_samples[i,]
  psi <- plogis(beta[1] + beta[2] * dataList$humidity)
  p <- plogis(beta[3] + beta[4] * dataList$humidity)
  qT <- (1-p)^dataList$nVisits

  zpost <- rep(NA, dataList$nSites)
  for (j in 1:dataList$nSites){
    if(dataList$no_detect[j] == 0){ #If detected at least once at site j
      zpost[j] <- 1
    } else{ # If never detected at site j
      psi_con <- psi[j] * qT[j]/(psi[j] * qT[j] + (1-psi[j]))
      zpost[j] <- rbinom(1, 1, psi_con)
    }
  }
  occ_fs_sims[i] <- sum(zpost)
}

# Plot distribution of bootstrap samples and compare truth and estimate
# (not shown)
par(mfrow = c(1, 1))
hist(occ_fs_sims, main = "Distribution of occ_fs", xlab = "occ_fs")
abline(v = mean(occ_fs_sims), col = 'blue', lwd = 4)
abline(v = true_Nz, col = 'red', lwd = 4)
legend("topleft", legend = c("Estimate", "Truth"), col = c("blue","red"), lwd = 3,
       bty = 'n')
tmb_Nz <- mean(occ_fs_sims)


# 19.10 Comparison of parameter estimates
# ---------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(truth = truth, unmarked = unm_est, JAGS = jags_est,
  NIMBLE = nimble_est, Stan = stan_est, ubms = ubms_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)

# Compare estimates of the number of occupied sites
comp <- c(truth = true_Nz, observed = obs_Nz, unmarked = unm_Nz, JAGS = jags_Nz,
  NIMBLE = nimble_Nz, Stan = stan_Nz, ubms = ubms_Nz, TMB = tmb_Nz)
print(comp, 2)

# 19.11 Summary and outlook
# -------------------------
# (no code)
