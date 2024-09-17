
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# ------------------------------------------
# Chapter 5  --  Normal linear regression
# ------------------------------------------

# Last changes: 11 June 2024


# 5.1 Introduction
# ----------------
# (no code)


# 5.2 Data generation
# -------------------

set.seed(5)
n <- 16 # Number of years
a <- 40 # Intercept
b <- -0.5 # Slope
sigma2 <- 25 # Residual variance
# Save true values for later comparisons
truth <- c(alpha = a, beta = b, sigma = sqrt(sigma2))

# See Fig. 5–3 (left) for a variant of this plot
x <- 1:16 # Values of covariate year
eps <- rnorm(n, mean = 0, sd = sqrt(sigma2)) # Residuals
y <- a + b*x + eps # Assemble data set
plot((x + 1989), y, xlab = "Year", las = 1,
ylab = "Prop. occupied (%)", cex = 1.5,
pch = 16, col = rgb(0,0,0,0.6), frame = FALSE) # not shown

# Load required libraries
library(ASMbook); library(jagsUI); library(rstan); library(TMB)


# 5.3 Analysis with canned functions in R
# ---------------------------------------
# (no code)


#   5.3.1 Fitting the model
# -------------------------

# Fit model and save lm estimates
summary(out5.3 <- lm(y ~ x))
lm_est <- c(coef(out5.3), sigma = sigma(out5.3))


#   5.3.2 Goodness-of-fit assessment using traditional and quantile residuals
# ---------------------------------------------------------------------------

# Print out traditional residuals (=raw residuals)
print(residuals(out5.3), 4)

# Visually assess adequacy of fitted linear regression model
hist(residuals(out5.3), breaks = 20) # not shown
plot(out5.3) # not shown

# Computation and goodness-of-fit assessment using quantile residuals
library(DHARMa)
simOut <- simulateResiduals(out5.3, n = 1000)
residuals(simOut) # print out scaled quantile residuals

# Normal-transformed quantile residuals
print(residuals(simOut, quantileFunction = qnorm), 4)

par(mfrow = c(1,2))
plotQQunif(simOut)
plotResiduals(simOut, rank = FALSE) # See comment above
par(mfrow = c(1,1))
testDispersion(simOut) # Test for over- or underdispersion
# (not shown, and note this test really
# makes sense only for nonnormal GLMs)


#   5.3.3 Computing predictions
# -----------------------------

xpred <- seq(1, 16, length.out = 100) # Predict for 100 values of x
newdata <- list(x = xpred)
pred <- predict(out5.3, newdata = newdata, se.fit = TRUE) # also get SEs
LCL <- pred$fit - 1.96 * pred$se.fit # Lower 95% confidence limit
UCL <- pred$fit + 1.96 * pred$se.fit # Upper 95% confidence limit
pred.table <- cbind(x = newdata$x, "point_est" =
pred$fit, "se" = pred$se.fit, "LCL" = LCL, "UCL" = UCL) # Table with x and predictions
print(pred.table) # Shown only partially,
# but see also Fig. 5.4


# Section 5.4 Bayesian analysis with JAGS
# ---------------------------------------
# (no code)


#   5.4.1 Fitting the model
# -------------------------

# Bundle and summarize data
str(dataList <- list(y = y, x = x, n = n))

# Write JAGS model file
cat(file = "model5.4.txt", "
model {
# Priors
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
sigma ~ dunif(0, 100) # May consider smaller range
tau <- pow(sigma, -2)
# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], tau) # Response y distributed normally
  mu[i] <- alpha + beta * x[i] # Linear predictor
}
# Assess model fit using a sums-of-squares-type discrepancy
for (i in 1:n) {
  residual[i] <- y[i]-mu[i] # Raw residuals for observed data
  predicted[i] <- mu[i] # Predicted, or fitted, values
  sq[i] <- residual[i]^2 # Squared residuals for observed data
  # Posterior predictive check
  # based on posterior predictive distribution of the data
  # Generate replicate data and compute fit stats for them
  y.new[i] ~ dnorm(mu[i], tau) # one new data set at each MCMC iteration:
  # this is the posterior predictive distribution of y
  sq.new[i] <- (y.new[i]-predicted[i])^2 # Squared resi. for new data
}
fit <- sum(sq[]) # Sum of squared residuals for actual data set
fit.new <- sum(sq.new[]) # Sum of squared residuals for new data set
}
")

# Function to generate starting values
inits <- function(){ list(alpha = rnorm(1),
beta = rnorm(1), sigma = rlnorm(1))}
# Note lognormal init for sigma
# to avoid values < 0
# Parameters to estimate
params <- c("alpha","beta", "sigma", "fit", "fit.new", "residual", "predicted", "y.new")
# MCMC settings
na <- 1000 ; ni <- 6000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# Call JAGS (ART <1 min), check convergence and summarize posteriors
out5.4 <- jags(dataList, inits, params, "model5.4.txt", n.iter = ni, n.burnin = nb,
n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
jagsUI::traceplot(out5.4) # not shown
print(out5.4, 2)

# Save JAGS estimates
jags_est <- out5.4$summary[1:3, 1]

# Compare results of lm() and JAGS
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est)
print(comp, 4)


#   5.4.2 Goodness-of-fit assessment in Bayesian analyses
# -------------------------------------------------------

# Fig. 5-3 (left)
par(mfrow = c(1, 2), mar = c(5,5,2,3))
plot(out5.4$mean$predicted, out5.4$mean$residual, main = "", las = 1,
xlab = "Predicted values", ylab = "Residuals", frame = FALSE,
pch = 16, col = rgb(0,0,0,0.5), cex = 1.5)
abline(h = 0)

# Fig. 5-3 (right)
lim <- c(0, 4300)
plot(out5.4$sims.list$fit, out5.4$sims.list$fit.new, main = "", las = 1,
xlab = "SSQ for actual data", ylab = "SSQ for new data", xlim = lim,
ylim = lim, frame = FALSE, pch = 16, col = rgb(0,0,0,0.3), cex = 1.5)
abline(0, 1)
mean(out5.4$sims.list$fit.new > out5.4$sims.list$fit) # Bayes. p-value

# How to use DHARMa for quantile residuals when fitting model in JAGS
y.new <- out5.4$sims.list$y.new
sim <- createDHARMa(simulatedResponse = t(y.new), observedResponse = y,
fittedPredictedResponse = apply(y.new, 2, median))
par(mfrow = c(1,2)) # Graphical results as in Fig. 5.2 (not shown)
plotQQunif(sim)
plotResiduals(sim, rank = FALSE)


#   Section 5.4.3 Computing predictions
# -------------------------------------

# Get predictions for Bayesian model fit
sims <- out5.4$sims.list # Grab posterior draws produced by JAGS
predictions <- array(dim = c(length(xpred), length(sims$alpha)))
for(i in 1:length(xpred)){
predictions[i,] <- sims$alpha + sims$beta * xpred[i]
}
BPE <- rowMeans(predictions) # Bayesian point estimates (post. means)
LPB <- apply(predictions, 1, quantile, probs = 0.025) # Lower bound
UPB <- apply(predictions, 1, quantile, probs = 0.975) # Upper bound

# Plot predictions with 95% uncertainty intervals from both model fits
# Fig. 5.4
plot(1990:2005, y, xlab = "Year", las = 1, ylab = "Prop. occupied (%)", pch = 16,
ylim = c(25, 50), col = rgb(0,0,0,0.5), frame = FALSE, cex = 1.5)
# CIs from frequentist inference with lm() (Section 5.3.3)
lines(xpred + 1989, pred.table[,2], lwd = 2, col = 'blue')
polygon(c(newdata$x + 1989, rev(newdata$x + 1989)), c(pred.table[,4],
rev(pred.table[,5])), col = rgb(0, 0, 1, 0.1), border = NA)
# CRIs from Bayesian inference (here)
lines(xpred + 1989, BPE, lwd = 2, col = 'red')
polygon(c(xpred+ 1989, rev(xpred+ 1989)), c(LPB, rev(UPB)), col = rgb(1, 0, 0, 0.2), border = NA)
legend('bottomleft', legend = c("Maximum likelihood", "Posterior inference"),
cex = 1.2, bty = 'n', lty = 1, col = c('blue', 'red'), lwd = 3)


#   5.4.4 Interpretation of confidence intervals versus credible intervals
# ------------------------------------------------------------------------

par(mfrow = c(1, 1), mar = c(5,5,2,3)) # Fig. 5–5
hist(out5.4$sims.list$beta, main = "", col = "grey", xlab = "Trend estimate",
xlim = c(-2, 1), breaks = 100, freq = FALSE)
abline(v = 0, col = "black", lwd = 2, lty = 3)
mean(out5.4$sims.list$beta < 0) # Probability of population decline


#   5.4.5 And what about prediction intervals?
# --------------------------------------------

# Get posterior predictive distribution of y at fine resolution
sims <- out5.4$sims.list # Grab all posterior draws
str(sims) # Remind ourselves of MCMC output format
xpred <- seq(1, 16, length.out = 1000) # Predict for 1000 values of x
y_tilde <- array(NA, dim = c(1000, 20000)) # Post. predictive dist.
for(k in 1:1000){ # Loop over all 1000 values of xpred
mu_k <- sims$alpha + sims$beta * xpred[k]
y_tilde[k,] <- rnorm(20000, mean = mu_k, sd = sims$sigma)
}
# Get 95% prediction interval for the data y_tilde
LPI <- apply(y_tilde, 1, quantile, probs = 0.025) # Lower bound
UPI <- apply(y_tilde, 1, quantile, probs = 0.975) # Upper bound
# Plot 95% prediction interval for the data (Fig. 5.6)
plot(1990:2005, y, xlab = "Year", las = 1, ylab = "Prop. occupied (%)",
pch = 16, ylim = c(18, 58), col = rgb(0,0,0,0.5), frame = FALSE,cex = 1.5)
polygon(c(xpred + 1989, rev(xpred + 1989)), c(LPI, rev(UPI)),
col = rgb(1, 0, 0, 0.2), border = NA)


# 5.5 Bayesian analysis with NIMBLE
# ---------------------------------

library(nimble)

# Bundle and summarize data
str(dataList <- list(y = y, x = x, n = n)) # not shown

# Write NIMBLE model file
model5.5 <- nimbleCode( {
# Priors
alpha ~ dnorm(0, sd = 100)
beta ~ dnorm(0, sd = 100)
sigma ~ dunif(0, 100)

# Likelihood
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], sd = sigma) 
  mu[i] <- alpha + beta * x[i]
}

# Assess model fit using a sums-of-squares-type discrepancy
for (i in 1:16) {
  residual[i] <- y[i]-mu[i]     # Residuals for observed data
  predicted[i] <- mu[i]         # Predicted values
  sq[i] <- residual[i]^2        # Squared residuals for observed data

# Generate replicate data and compute fit stats for them
  y.new[i] ~ dnorm(mu[i], sd = sigma) # one new data set at each MCMC iteration
  sq.new[i] <- (y.new[i]-predicted[i])^2 # Sq. residuals for new data
}
fit <- sum(sq[1:16])     # Sum of squared residuals for actual data set
fit.new <- sum(sq.new[1:16])# Sum of squared residuals for new data set
} )

# Inits
inits <- function(){ list(alpha = rnorm(1), beta = rnorm(1), sigma = rlnorm(1))}

# Parameters monitored: same as before
params <- c("alpha", "beta", "sigma", "fit", "fit.new",  "residual", "predicted")

# MCMC settings
#  The number of samples returned will be floor((niter-nburnin)/thin).
ni <- 6000  ;  nb <- 1000  ; nc <- 4  ; nt <- 1

# Call NIMBLE (ART 25 sec), check convergence and summarize posteriors
system.time(  
  out5.5 <- nimbleMCMC(code = model5.5,
    constants = dataList, inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2)); coda::traceplot(out5.5)     # not shown
(nsum <- nimble_summary(out5.5, params))       # not shown

# Bayesian p-value
mean(as.matrix(out5.5[,"fit.new"]) > as.matrix(out5.5[,"fit"]))

# Save NIMBLE estimates
nimble_est <- nsum[1:3,1]


# 5.6 Bayesian analysis with Stan
# -------------------------------

# Bundle and summarize data (same as before)
str(dataList <- list(y = y, x = x, n = n)) # not shown

# Write text file with model description in Stan language
cat(file = "model5_6.stan", # This line is R code
"data { // This is the first line of Stan code
  int <lower = 0> n; // Define the format of all data
  vector[n] y; //. . . including the dimension of vectors
  vector[n] x; //
}
parameters { // Define format for all parameters
  real alpha;
  real beta;
  real <lower = 0> sigma; // sigma (sd) cannot be negative
}
transformed parameters {
  vector[n] mu;
  for (i in 1:n){
    mu[i] = alpha + beta * x[i]; // Calculate linear predictor
  }
}
model {
  // Priors
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 10);
  // Likelihood (could be vectorized to increase speed)
  for (i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
}
generated quantities {
  vector[n] residuals;
  vector[n] sq;
  vector[n] sq_new;
  vector[n] y_new;
  real fit;
  real fit_new;
  for (i in 1:n){
    residuals[i] = y[i] - mu[i];
    sq[i] = residuals[i]^2;
    y_new[i] = normal_rng(mu[i], sigma);
    sq_new[i] = (y_new[i] - mu[i])^2;
  }
  fit = sum(sq);
  fit_new = sum(sq_new);
} // This is the last line of Stan code
" )

# HMC settings
ni <- 3000 ; nb <- 1000 ; nc <- 4 ; nt <- 1
# Call STAN (ART 30/3 sec), assess convergence, save estimates
system.time(
out5.6 <- stan(file = "model5_6.stan", data = dataList,
chains = nc, iter = ni, warmup = nb, thin = nt))
rstan::traceplot(out5.6) # not shown
print(out5.6, dig = 2) # not shown
stan_est <- summary(out5.6)$summary[1:3,1]


# Section 5.7 Do-it-yourself maximum likelihood estimation
# --------------------------------------------------------

NLL <- function(params, y, x) {
alpha <- params[1]
beta <- params[2]
sigma <- exp(params[3]) # convert from log scale
n <- length(y) # number of datapoints
mu <- LL <- numeric(n) # empty vectors
# You could vectorize this loop to speed things up
for (i in 1:n){
mu[i] <- alpha + beta * x[i]
LL[i] <- dnorm(y[i], mean = mu[i], sd = sigma, log = TRUE)
}
-sum(LL)
}

# Minimize that NLL to find MLEs and also get SEs
inits <- c('alpha' = 0, 'beta' = 0, 'log.sigma' = 0)
out5.7 <- optim(inits, NLL, y = y, x = x, hessian = TRUE, method = 'BFGS')
get_MLE(out5.7, 4)
# Save DIY estimates
diy_est <- c(out5.7$par[1:2], exp(out5.7$par[3]))


# 5.8 Likelihood analysis with TMB
# --------------------------------

# Bundle and summarize data (same as before)
str(tmbData <- list(y = y, x = x, n = n)) # not shown

# Write TMB model file
cat(file = "model5_8.cpp",
"#include <TMB.hpp> // (1) First some boilerplate code

template <class Type>
Type objective_function <Type>::operator() ()
{
  // (2) Describe input data
  DATA_VECTOR(y); // response
  DATA_VECTOR(x); // covariate
  DATA_INTEGER(n); // Number of obs

  // (3) Describe parameters
  PARAMETER(alpha); // Intercept
  PARAMETER(beta); // Slope
  PARAMETER(log_sigma); // log(residual standard deviation)
  Type sigma = exp(log_sigma); // Type = match type of function output

  // (4) Calculate the log-likelihood
  Type LL = 0.0; // Initialize total log likelihood at 0
  Type mu; //Initialize the expected value

  for (int i = 0; i < n; i++){ // Note index starts at 0 instead of 1!
    mu = alpha + beta * x(i);
    LL += dnorm(y(i), mu, sigma, true);
    //Equivalent to LL = LL + dnorm(. . .)
  }

  return -LL; // Return negative of total log likelihood
}
")

# Compile and load TMB function
compile("model5_8.cpp")
dyn.load(dynlib("model5_8"))

# Provide dimensions, names and starting values for parameters
params <- list(alpha = 0, beta = 0, log_sigma = 0)

# Create TMB object
out5.8 <- MakeADFun(data = tmbData, parameters = params,
DLL = "model5_8", silent = TRUE)

# Optimize TMB object and print results
starts <- rep(0, length(unlist(params)))
opt <- optim(starts, fn = out5.8$fn, gr = out5.8$gr, method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out5.8))

# Save TMB estimates
tmb_est <- c(opt$par[1:2], exp(opt$par[3]))


# 5.9 Comparison of parameter estimates
# -------------------------------------

# Compare results of all engines with truth
comp <- cbind(truth = truth, lm = lm_est, JAGS = jags_est, NIMBLE = nimble_est,
Stan = stan_est, DIY = diy_est, TMB = tmb_est)
print(comp, 4)


# 5.10 Summary and outlook
# ------------------------
# (no code)
