
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# --------------------------------------------------
# Chapter 16  --  Binomial GLM with continuous and
#                 categorical explanatory variables
# --------------------------------------------------

# Last changes: 11 June 2024


# 16.1 Introduction
# -----------------
# (no code)


# 16.2 Data generation
# --------------------

set.seed(16)
nRegion <- 3
nSite <- 10
n <- nRegion * nSite
x <- rep(1:nRegion, rep(nSite, nRegion))
region <- factor(x, labels = c("Jura", "Black Forest", "Alps"))

wetness.Jura <- sort(runif(nSite, 0, 1))
wetness.BlackF <- sort(runif(nSite, 0, 1))
wetness.Alps <- sort(runif(nSite, 0, 1))
wetness <- c(wetness.Jura, wetness.BlackF, wetness.Alps)

N <- round(runif(n, 10, 50) )        # Get discrete uniform values for N

Xmat <- model.matrix(~ region * wetness)

print(Xmat, dig = 2)                 # not shown, but make sure you understand this!
truth <- beta.vec <- c(-4, 1, 2, 6, 2, -5)

lin.pred <- Xmat[,] %*% beta.vec     # Value of lin.predictor
exp.p <- exp(lin.pred) / (1 + exp(lin.pred)) # Expected proportion, note this is
                                             # same as plogis(lin.pred) in R
C <- rbinom(n = n, size = N, prob = exp.p)   # Add binomial noise
hist(C)                              # Inspect simulated binomial counts

# Draw Fig. 16.2
par(mfrow = c(1,2), mar = c(5,5,3,1))
matplot(cbind(wetness[1:10], wetness[11:20], wetness[21:30]), cbind(exp.p[1:10],
  exp.p[11:20], exp.p[21:30]), ylab = "Expected proportion black", xlab = "Wetness index",
  col = c("red","green","blue"), pch = c("J","B","A"), lty = "solid", type = "b", las = 1,
  cex = 1.2, main = "Expected proportion", lwd = 2, frame = FALSE)
matplot(cbind(wetness[1:10], wetness[11:20], wetness[21:30]), cbind(C[1:10]/N[1:10],
  C[11:20]/N[11:20], C[21:30]/N[21:30]), ylab = "Observed proportion black",
  xlab = "Wetness index", col = c("red","green","blue"), pch = c("J","B","A"), las = 1,
  cex = 1.2, main = "Realized proportion", frame = FALSE)

# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 16.3 Likelihood analysis with canned functions in R
# ---------------------------------------------------

summary(out16.3 <- glm(cbind(C, N-C) ~ region * wetness, family = binomial))
glm_est <- coef(out16.3)             # Save estimates

# Diagnostic checks of residuals/model GOF (not shown)
plot(out16.3)                        # Traditional Pearson residual check

# Compute quantile residuals based on simulation
library(DHARMa)
simOut <- simulateResiduals(out16.3, n = 1000, plot = TRUE)
plotResiduals(simOut, form = wetness)# same as right plot ...
testDispersion(out16.3)              # Test for over- or underdispersion


# 16.4 Bayesian analysis with JAGS
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, N = N, nRegion = nRegion,
  region = as.numeric(region),wetness = wetness, n = n) )

# Write JAGS model file
cat(file = "model16.4.txt", "
model {
# Priors
for (i in 1:nRegion){
  alpha[i] ~ dnorm(0, 0.0001)        # Intercepts
  beta[i] ~ dnorm(0, 0.0001)         # Slopes
}

# Likelihood
for (i in 1:n) {
  C[i] ~ dbin(p[i], N[i])
  logit(p[i]) <- alpha[region[i]] + beta[region[i]]*wetness[i]    # Jura is baseline
  # p[i] <- ilogit(alpha[region[i]] + beta[region[i]]*wetness[i]) # same !

  # Fit assessments: Pearson residuals and posterior predictive check
  Presi[i] <- (C[i]-N[i]*p[i])/sqrt(N[i]*p[i]*(1-p[i]))           # Pearson resi
  C.new[i] ~ dbin(p[i], N[i])                                     # Create replicate data set
  Presi.new[i] <- (C.new[i]-N[i]*p[i])/sqrt(N[i]*p[i]*(1-p[i]))
}

# Add up squared residual as our discrepancy measures
fit <- sum(Presi[]^2)
fit.new <- sum(Presi.new[]^2)

# Derived quantities
# Recover the effects relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1]       # Intercept Black Forest vs. Jura
a.effe3 <- alpha[3] - alpha[1]       # Intercept Alps vs. Jura
b.effe2 <- beta[2] - beta[1]         # Slope Black Forest vs. Jura
b.effe3 <- beta[3] - beta[1]         # Slope Alps vs. Jura
# Custom comparison
test1 <- beta[3] - beta[2]           # Difference slope Alps -Black Forest
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(nRegion, 3, 1), beta = rnorm(nRegion, 2, 1))}

# Parameters to estimate (add "C.new" for DHARMa GoF test)
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", "b.effe3",
  "test1", "Presi", "fit", "fit.new", "C.new")

# MCMC settings
na <- 5000 ; ni <- 50000 ; nb <- 10000 ; nc <- 4 ; nt <- 40

# Call JAGS (ART <1 min), check convergence, summarize posteriors and save results
out16.4 <- jags(dataList, inits, params, "model16.4.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out16.4)           # not shown
print(out16.4, 3)                    # shown partially
jags_est <- unlist(out16.4$mean)[c(1, 7, 8, 4, 9, 10)]

mean(out16.4$sims.list$fit.new > out16.4$sims.list$fit)

# Draw Fig. 16.3
par(mfrow = c(1, 3), cex = 1.5)
plot(out16.4$mean$Presi, ylab = "Residual", las = 1, pch = 16, frame = FALSE)
abline(h = 0)
plot(wetness, out16.4$mean$Presi, ylab = "Residual", las = 1, pch = 16, frame = FALSE)
abline(h = 0)
plot(out16.4$sims.list$fit, out16.4$sims.list$fit.new, main = "",
  xlab = "Discrepancy actual data", ylab = "Discrepancy ideal data",
col = rgb(0,0,0,0.3), pch = 16, cex = 1, frame = FALSE)
abline(0,1, lwd = 2, col = "black")

# Do quantile residual assessments (not shown)
C.new <- out16.4$sims.list$C.new
sim <- createDHARMa(simulatedResponse = t(C.new), observedResponse = C,
  fittedPredictedResponse = apply(C.new, 2, median), integerResponse = T)
plot(sim)

# Compare likelihood with Bayesian estimates and with truth
comp <- cbind(truth = beta.vec, glm = glm_est, JAGS = jags_est)
print(comp, 4)


# 16.5 Bayesian analysis with NIMBLE
# ----------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(C = C, N = N, nRegion = nRegion,
  region = as.numeric(region), wetness = wetness, n = n) )

# Write NIMBLE model file
model16.5 <- nimbleCode( {
# Priors
for (i in 1:nRegion){
  alpha[i] ~ dnorm(0, sd = 100)      # Intercepts
  beta[i] ~ dnorm(0, sd = 100)       # Slopes
}

# Likelihood
for (i in 1:n) {
  C[i] ~ dbin(p[i], N[i])
  logit(p[i]) <- alpha[region[i]] + beta[region[i]]* wetness[i] # Jura is baseline
# p[i] <- ilogit(alpha[region[i]] + beta[region[i]]*wetness[i]) # same

# Fit assessments: Pearson residuals and posterior predictive check
  Presi[i] <- (C[i]-N[i]*p[i]) / sqrt(N[i]*p[i]*(1-p[i]))       # Pearson resi
  C.new[i] ~ dbin(p[i], N[i])                                   # Create replicate data set
  Presi.new[i] <- (C.new[i]-N[i]*p[i]) / sqrt(N[i]*p[i]*(1-p[i]))
}

# Add up squared Pearson residuals as our discrepancy measures
fit <- sum(Presi[1:n]^2)
fit.new <- sum(Presi.new[1:n]^2)

# Derived quantities
# Recover the effects relative to baseline level (no. 1)
a.effe2 <- alpha[2] - alpha[1]       # Intercept Black Forest vs. Jura
a.effe3 <- alpha[3] - alpha[1]       # Intercept Alps vs. Jura
b.effe2 <- beta[2] - beta[1]         # Slope Black Forest vs. Jura
b.effe3 <- beta[3] - beta[1]         # Slope Alps vs. Jura

# Custom comparison
test1 <- beta[3] - beta[2]           # Difference slope Alps -Black Forest
} )


# Inits
inits <- function(){ list(alpha = rnorm(nRegion, 3, 1), 
  beta = rnorm(nRegion, 2, 1))}

# Parameters monitored: same as before (except for C.new)
params <- c("alpha", "beta", "a.effe2", "a.effe3", "b.effe2", 
  "b.effe3", "test1", "Presi", "fit", "fit.new")

# MCMC settings
ni <- 50000  ;  nb <- 10000  ; nc <- 4  ; nt <- 40

# Call NIMBLE (ART 45 sec), check convergence, summarize and save posteriors
system.time( out16.5 <- 
    nimbleMCMC(code = model16.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out16.5)
(nsum <- nimble_summary(out16.5, params))     # not shown
nimble_est <- nsum[c(1, 7, 8, 4, 9, 10), 1]   # save estimates


# 16.6 Bayesian analysis with Stan
# --------------------------------

# Bundle and summarize data
str(dataList <- list(C = C, N = N, nRegion = nRegion,
  region = as.numeric(region), wetness = wetness, n = n) )

# Write Stan model
cat(file = "model16_6.stan", "
data{
  int n;                             // Number of samples
  int nRegion;                       // Number of regions
  int C[n];                          // Counts for each sample
  int N[n];                          // Sizes of each sample
  int region[n];                     // Region indices
  vector[n] wetness;                 // Covariate
}
parameters{
  real alpha[nRegion];               // Intercepts for each region
  real beta[nRegion];                // Slopes for each region
}
transformed parameters{
  vector[n] p; // Estimated probabilities
  for (i in 1:n){
    p[i] = inv_logit(alpha[region[i]] + beta[region[i]] * wetness[i]);
  }
}

model{
  for (i in 1:nRegion){
    alpha[i] ~ normal(0, 100);
    beta[i] ~ normal(0, 100);
  }
  for (i in 1:n){
    C[i] ~ binomial(N[i], p[i]);
  }
}

generated quantities{
  real a_effe2 = alpha[2] - alpha[1]; // Intercept Black Forest vs. Jura
  real a_effe3 = alpha[3] - alpha[1]; // Intercept Alps vs. Jura
  real b_effe2 = beta[2] - beta[1]; // Slope Black Forest vs. Jura
  real b_effe3 = beta[3] - beta[1]; // Slope Alps vs. Jura
  real test1 = beta[3] - beta[2];   // Difference slope Alps-Black Forest
  int C_new[n];                     // New simulated dataset
  vector[n] Presi;                  // Pearson residuals for real dataset
  vector[n] Presi_new;              // Pearson residuals for simulated dataset
  vector[n] Presi2;                 // Squared Pearson resi. for real dataset
  vector[n] Presi2_new;             // Squared Pearson resi. for simulated dataset
  real fit;                         // Sum of Pearson residuals
  real fit_new;
  for (i in 1:n){
    Presi[i] = (C[i]-N[i]*p[i])/sqrt(N[i]*p[i]*(1-p[i]));
    Presi2[i] = pow(Presi[i], 2);
    C_new[i] = binomial_rng(N[i], p[i]);
    Presi_new[i] = (C_new[i]-N[i]*p[i])/sqrt(N[i]*p[i]*(1-p[i]));
    Presi2_new[i] = pow(Presi_new[i], 2);
  }
fit = sum(Presi2);                  // Add up squared discrepancies
fit_new = sum(Presi2_new);
}
")

# HMC settings
ni <- 2000 ; nb <- 1000 ; nc <- 4 ; nt <- 1

# Call STAN (ART 53/4 sec), assess convergence and print and save result
system.time(
out16.6 <- stan(file = "model16_6.stan", data = dataList,
warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out16.6)           # not shown
print(out16.6, dig = 3)             # not shown
stan_est <- summary(out16.6)$summary[c(1,37,38,4,39,40),1] # Save esti


# 16.7 Do-it-yourself MLEs
# ------------------------

# Define NLL for general logistic regression with Binomial response
NLL <- function(beta, y, N, Xmat) {
  p <- plogis(Xmat %*% beta)
  LL <- dbinom(y, N, p, log = TRUE) # log-likelihood contr. for each obs
  NLL <- -sum(LL)                   # NLL for all observations in data set
  return(NLL)
}

# Minimize that NLL to find MLEs, get SEs and CIs and save estimates
inits <- rep(0, 6)
names(inits) <- names(coef(out16.3))
out16.7 <- optim(inits, NLL, y = C, N = N, Xmat = Xmat, hessian = TRUE, method = "BFGS")
get_MLE(out16.7, 4)
diy_est <- out16.7$par              # Save estimates


# 16.8 Likelihood analysis with TMB
# ---------------------------------

# Bundle and summarize data
tmbData <- dataList
tmbData$region <- tmbData$region-1  # Indices start at 0 in TMB/C++
str(tmbData)

# Write TMB model file
cat(file = "model16_8.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Describe input data
  DATA_VECTOR(C);
  DATA_VECTOR(N);
  DATA_IVECTOR(region);
  DATA_VECTOR(wetness);
  DATA_INTEGER(n);

  // Describe parameters
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type LL = 0.0;                    // Initialize log-likelihood at 0

  vector<Type> p(n);
  vector<Type> Presi(n);

  for (int i= 0; i<n; i++){
    p(i) = invlogit(alpha(region(i)) + beta(region(i)) * wetness(i));
    LL += dbinom(C(i), N(i), p(i), true);
    Presi(i) = (C(i)-N(i)*p(i))/sqrt(N(i)*p(i)*(1-p(i)));
  }

  Type a_effe2 = alpha(1) - alpha(0); // Intercept Black Forest vs. Jura
  Type a_effe3 = alpha(2) - alpha(0); // Intercept Alps vs. Jura
  Type b_effe2 = beta(1) - beta(0);   // Slope Black Forest vs. Jura
  Type b_effe3 = beta(2) - beta(0);   // Slope Alps vs. Jura
  Type test1 = beta(2) - beta(1);     // Difference slope Alps-Black Forest
  ADREPORT(a_effe2);
  ADREPORT(a_effe3);
  ADREPORT(b_effe2);
  ADREPORT(b_effe3);
  ADREPORT(test1);
  ADREPORT(Presi);
  return - LL;
}
")

# Compile and load TMB function
compile("model16_8.cpp")
dyn.load(dynlib("model16_8"))

# Provide dimensions and starting values for parameters
params <- list(alpha = rep(0, tmbData$nRegion), beta = rep(0, tmbData$nRegion))

# Create TMB object
out16.8 <- MakeADFun(data = tmbData,
  parameters = params,
  DLL = "model16_8", silent = TRUE)

# Optimize TMB object and print and save results
opt <- optim(out16.8$par, fn = out16.8$fn, gr = out16.8$gr, method = "BFGS")
(tsum <- tmb_summary(out16.8))          # not shown
tmb_est <- tsum[c(1, 7, 8, 4, 9, 10),1] # save estimates


# 16.9 Comparison of the estimates
# --------------------------------

# Compare estimates with truth
comp <- cbind(truth = truth, glm = glm_est, JAGS =jags_est,
  NIMBLE = nimble_est, Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 3)


# 16.10 Summary
# -------------
# (no code)
