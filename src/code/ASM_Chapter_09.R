
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ----------------------------------------------------------------
# Chapter 9  --  General linear model for a normal response with 
#                continuous andcategorical explanatory variables
# ----------------------------------------------------------------

# Last changes: 11 June 2024


# 9.1 Introduction
# ----------------
# (no code)


# 9.2 Data generation
# -------------------

set.seed(9)
nPops <- 3
nSample <- 10
n <- nPops * nSample # Total number of data points
x <- rep(1:nPops, rep(nSample, nPops)) # Indicator for population
pop <- factor(x, labels = c("Pyrenees", "Massif Central", "Jura"))
length <- runif(n, 45, 70) # ad. body length rarely <45 cm
lengthC <- length-mean(length) # Use centered length

Xmat <- model.matrix(~ pop * lengthC)
print(Xmat, dig = 2) # not shown, but make sure to understand this !
beta.vec <- c(80, -30, -20, 6, -3, -4)

sigma <- 10 # Choose residual SD
lin.pred <- Xmat[,] %*% beta.vec # Value of lin.predictor
eps <- rnorm(n = n, mean = 0, sd = sigma) # residuals
mass <- as.numeric(lin.pred + eps) # response = lin.pred + residual
hist(mass) # Inspect what we’ve created (not shown)

par(mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
matplot(cbind(length[1:10], length[11:20], length[21:30]), cbind(mass[1:10],
mass[11:20], mass[21:30]), ylim = c(0, max(mass)), ylab = "Body mass (g)",
xlab = "Body length (cm)", col = c("Red","Green","Blue"), pch = c("P", "M", "J"),
las = 1, cex = 1.6, cex.lab = 1.5, frame = FALSE) # Fig. 9–2

# Save chosen values for later comparisons
truth <- c(beta.vec, sigma)

# Load required libraries
library(ASMbook); library(DHARMa); library(jagsUI); library(rstan); library(TMB)


# 9.3 Likelihood analysis with canned functions in R
# --------------------------------------------------

summary(out9.3 <- lm(mass ~ pop * lengthC))

# Save least-squares estimates
lm_est <- c(coef(out9.3), sigma = summary(out9.3)$sigma)

# Check goodness-of-fit using quantile residuals
plot(out9.3) # Traditional residual check for comparison
simOut <- simulateResiduals(out9.3, n = 1000, plot = TRUE)

par(mfrow = c(1,2))
# Plot effect of length, holding pop constant (left plot)
# make sequence of 100 length values
length.seq <- seq(min(length), max(length), length.out = 100)
# Center the length sequence using the original mean
lengthC.seq <- length.seq - mean(length)
# make newdata data.frame, predict and compute 95% CI around prediction
newdata <- data.frame(lengthC = lengthC.seq, pop = "Pyrenees")
pr.length <- predict(out9.3, newdata = newdata, se.fit = TRUE)
LCL <- pr.length$fit - 1.96 * pr.length$se.fit # Wald 95% LCL
UCL <- pr.length$fit + 1.96 * pr.length$se.fit # Wald 95% UCL

# Draw the plot (Fig. 9-3 left)
plot(length[pop == "Pyrenees"], mass[pop == "Pyrenees"], pch = 16,
ylab = "Mass", xlab = "Length", main = "Effect of length on mass for Pyrenees",
frame = FALSE, ylim = c(0, 180), cex = 1.5, col = rgb(0,0,0,0.5))
lines(length.seq, pr.length$fit, col = 'red', lwd = 3)
polygon(c(length.seq, rev(length.seq)), c(LCL, rev(UCL)), col = rgb(1,0,0,0.2), border = NA)

# Plot effect of pop, holding length constant at its median (right plot)
newdata <- data.frame(lengthC = median(lengthC), pop = levels(pop))
pr.pop <- predict(out9.3, newdata = newdata, se.fit = TRUE)
LCL <- pr.pop$fit - 1.96 * pr.pop$se.fit
UCL <- pr.pop$fit + 1.96 * pr.pop$se.fit
plot(as.numeric(pop), mass, pch = 16, xaxt = 'n', xlab = "Pop",
ylab = "Mass", main = "Effect of pop on mass at median length",
frame = FALSE, cex = 1.5, col = rgb(0,0,0,0.4)) # Draw the plot (Fig. 9.3 right)
axis(1, at = 1:3, labels = levels(pop))
points(1:3, pr.pop$fit, pch = 16, col = "red", cex = 1.5)
segments(1:3, LCL, 1:3, UCL, lwd = 2, col = "red")
par(mfrow = c(1,1))


# 9.4 Bayesian analysis with JAGS
# -------------------------------

# Bundle and summarize data
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC,
    nPops = nPops, n = n) )

# Write JAGS model file
cat(file = "model9.4.txt", "
model {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, 0.0001) # Intercepts
  beta[i] ~ dnorm(0, 0.0001) # Slopes
}
sigma ~ dunif(0, 100) # Residual standard deviation
tau <- pow(sigma, -2)

# Likelihood
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthC[i]
}

# Define effects relative to baseline level
a.effe2 <- alpha[2] - alpha[1] # Intercept Massif Central vs. Pyr.
a.effe3 <- alpha[3] - alpha[1] # Intercept Jura vs. Pyr.
b.effe2 <- beta[2] - beta[1] # Slope Massif Central vs. Pyr.
b.effe3 <- beta[3] - beta[1] # Slope Jura vs. Pyr.
# Custom comparison
test1 <- beta[3] - beta[2] # Slope Jura vs. Massif Central
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(nPops, 0, 2), beta = rnorm(nPops, 1, 1), sigma = runif(1))}

# Parameters to estimate
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMC settings
na <- 1000 ; ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call JAGS (ART 2 sec), check convergence, summarize posteriors and save estimates
out9.4 <- jags(dataList, inits, params, "model9.4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out9.4) # not shown
print(out9.4, 3)
jags_est <- out9.4$summary[c(1, 8, 9, 4, 10, 11, 7), 1]

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


# 9.5 Bayesian analysis with NIMBLE
# ---------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC, nPops = nPops, n = n) )

# Write NIMBLE model file
model9.5 <- nimbleCode( {
# Priors
for (i in 1:nPops){
  alpha[i] ~ dnorm(0, 0.0001)     # Intercepts
  beta[i] ~ dnorm(0, 0.0001)      # Slopes
}
sigma ~ dunif(0, 100)             # Residual standard deviation
tau <- pow(sigma, -2)

# Likelihood
for (i in 1:n) {
  mass[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha[pop[i]] + beta[pop[i]] * lengthC[i]
}

# Define effects relative to baseline level
a.effe2 <- alpha[2] - alpha[1]    # Intercept Massif Central vs. Pyr.
a.effe3 <- alpha[3] - alpha[1]    # Intercept Jura vs. Pyr.
b.effe2 <- beta[2] - beta[1]      # Slope Massif Central vs. Pyr.
b.effe3 <- beta[3] - beta[1]      # Slope Jura vs. Pyr.

# Custom comparison
test1 <- beta[3] - beta[2]        # Slope Jura vs. Massif Central
} )

# Inits
inits <- function(){ list(alpha = rnorm(nPops, 0, 2), beta = rnorm(nPops, 1, 1), sigma = runif(1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# MCMC settings
ni <- 3000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 25 sec), check convergence, summarize posteriors and save estimates
system.time( out9.5 <- 
    nimbleMCMC(code = model9.5,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE) )
par(mfrow=c(2,2)); coda::traceplot(out9.5)
(nsum <- nimble_summary(out9.5, params))
nimble_est <- nsum[c(1, 8, 9, 4, 10, 11, 7), 1]


# 9.6 Bayesian analysis with Stan
# -------------------------------

# Bundle and summarize data
str(dataList <- list(mass = as.numeric(mass), pop = as.numeric(pop), lengthC = lengthC,
  nPops = nPops, n = n) )

# Write Stan model
cat(file = "model9_6.stan", "
data {
  int n; //sample size
  int nPops; //number of populations
  vector[n] mass; //response
  vector[n] lengthC; //covariate
  array[n] int pop; //population of each observation
}

parameters {
  vector[nPops] alpha; //intercepts
  vector[nPops] beta; //slopes
  real <lower = 0> sigma; //residual standard deviation
}

model {
  vector[n] mu; //expected value of observations
  
  //Priors
  for (i in 1:nPops){
    alpha[i] ~ normal(0, 100);
    beta[i] ~ normal(0, 100);
  }
  sigma ~ uniform(0, 100);

  //Likelihood
  for (i in 1:n){
    mu[i] = alpha[pop[i]] + beta[pop[i]] * lengthC[i];
    mass[i] ~ normal(mu[i], sigma);
  }
}

generated quantities {
  real a_effe2;
  real a_effe3;
  real b_effe2;
  real b_effe3;
  real test1;
  a_effe2 = alpha[2] - alpha[1]; //Intercept Massif Central vs. Pyr.
  a_effe3 = alpha[3] - alpha[1]; //Intercept Jura vs. Pyr.
  b_effe2 = beta[2] - beta[1]; //Slope Massif Central vs. Pyr.
  b_effe3 = beta[3] - beta[1]; //Slope Jura vs. Pyr.
  test1 = beta[3] - beta[2]; //Slope Jura vs. Massif Central
}
")

# Parameters to estimate
params <- c("alpha", "beta", "sigma", "a.effe2", "a.effe3", "b.effe2", "b.effe3", "test1")

# HMC settings
ni <- 1000 ; nb <- 500 ; nc <- 4 ; nt <- 1

# Call STAN (ART 40/2 sec), print output and save estimates
system.time(
out9.6 <- stan(file = "model9_6.stan", data = dataList,
  warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out9.6) # not shown
print(out9.6, dig = 3) # not shown
stan_est <- summary(out9.6)$summary[c(1, 8, 9, 4, 10, 11, 7), 1]


# 9.7 Do-it-yourself MLEs
# -----------------------

NLL <- function(param, y, Xmat) {
  beta <- param[1:6] # Intercept/slopes matching columns of Xmat
  sigma <- exp(param[7]) # Residual SD
  mu <- Xmat %*% beta # Matrix-multiply beta and design matrix
  LL <- dnorm(y, mu, sigma, log = TRUE) # Log-lik for each obs
  -sum(LL) # NLL for all observations (whole data set)
}

# Minimize NLL to find MLEs, get SEs and 95% CIs and save estimates
inits <- c(mean(mass), rep(0, 5), log(sd(mass)))
names(inits) <- c(colnames(Xmat), 'log-sigma')
out9.7 <- optim(inits, NLL, y = mass, Xmat = Xmat, hessian = TRUE, method = "BFGS")
get_MLE(out9.7, 4)
diy_est <- c(out9.7$par[1:6], exp(out9.7$par[7]))


# 9.8 Likelihood analysis with TMB
# --------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$pop <- tmbData$pop-1
str(tmbData)

# Write TMB model file
cat(file = "model9_8.cpp",
"#include <TMB.hpp>
template <class Type>
Type objective_function <Type>::operator() ()
{
  //Describe input data
  DATA_INTEGER(n); //number of obs
  DATA_VECTOR(mass); //response
  DATA_VECTOR(lengthC); //covariate
  DATA_IVECTOR(pop); //population index

  //Describe parameters
  PARAMETER_VECTOR(alpha); //Intercepts
  PARAMETER_VECTOR(beta); //Slopes
  PARAMETER(log_sigma); //Residual sd on log scale

  Type sigma = exp(log_sigma); //Residual SD
  
  vector <Type> mu(n); //Expected value of observations
  
  Type LL = 0; // Initialize log-likelihood at 0

  //Iterate over observations
  for (int i = 0; i < n; i++){ //Note index starts at 0 instead of 1
    mu(i) = alpha(pop(i)) + beta(pop(i)) * lengthC(i);
    //Calculate log-likelihood of observation and add to total loglik
    LL += dnorm(mass(i), mu(i), sigma, true);
  }

  //Calculate derived parameters
  Type a_effe2 = alpha(1) - alpha(0); //Intercept M. Central vs. Pyr.
  Type a_effe3 = alpha(2) - alpha(0); //Intercept Jura vs. Pyr.
  Type b_effe2 = beta(1) - beta(0); //Slope Massif Central vs. Pyr.
  Type b_effe3 = beta(2) - beta(0); //Slope Jura vs. Pyr.
  Type test1 = beta(2) - beta(1); //Slope Jura vs. Massif Central

  //Tell TMB to report these derived parameters
  ADREPORT(a_effe2);
  ADREPORT(a_effe3);
  ADREPORT(b_effe2);
  ADREPORT(b_effe3);
  ADREPORT(test1);
  
  return -LL; //Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model9_8.cpp")
dyn.load(dynlib("model9_8"))

# Provide dimensions and starting values for parameters
# Starting values need to be semi-informative (use same ones as DIY)
params <- list(alpha = c(mean(mass), 0, 0), beta = rep(0,3),
log_sigma = log(sd(mass)))

# Create TMB object
out9.8 <- MakeADFun(data = tmbData, parameters = params,
  DLL = "model9_8", silent = TRUE)

# Optimize TMB object and print and save results
opt <- optim(out9.8$par, fn = out9.8$fn, gr = out9.8$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out9.8))
tmb_est <- c(tsum[c(1, 8, 9, 4, 10, 11), 1], exp(tsum[7, 1]))


# 9.9 Comparison of the parameter estimates
# -----------------------------------------

# Compare results with truth and previous estimates
comp <- cbind(cbind(truth = truth, lm = lm_est, JAGS = jags_est,
NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est))
print(comp, 3)


# 9.10 Summary and outlook
# ------------------------
# (no code)
