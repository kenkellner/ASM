
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -------------------------------------------------------------------
# Chapter 18  --  Model building, model checking, and model selection
# -------------------------------------------------------------------

# This is the complete code for this book chapter. It contains
# all the code in the printed book, plus some extended code,
# the results of which are described and discussed in the chapter,
# but which is not or only incompletely shown in the printed book.
# This extra code appears in Sections 18.6 and 18.8.


# Last changes: 8 July 2024

# Load required packages
# ----------------------

library(ASMbook)
library(DHARMa)
library(jagsUI)


# 18.1 Introduction
# ---------------
# (no code)


# 18.2 Why do we build a statistical model?
# -----------------------------------------
# (no code)


# 18.3 Ticklish rattlesnakes 
# --------------------------
# (no code)


# 18.4 How do we build a statistical model? 
# -----------------------------------------
# (no code)


#   18.4.1 Model expansion and iterating on a single model
# --------------------------------------------------------
# (no code)


#   18.4.2 Causal, path, and structural equation models
# -----------------------------------------------------
# (no code)


#   18.4.3 How to build a predictive model
# ----------------------------------------
# (no code)


# 18.5 Model checking and goodness-of-fit testing
# -----------------------------------------------
# (no code)


#   18.5.1 Traditional residual diagnostics
# -----------------------------------------

dat <- simDat9()                     # Simulate data set like in Chapter 9
fm1 <- lm(dat$mass ~ dat$pop)        # Fit model with population factor only
fm2 <- lm(dat$mass ~ dat$lengthC)    # Fit model with length only
plot(fm1)                            # Produce plots with residual diagnostics
plot(fm2)

set.seed(18)
str(dat <- simDat18(nSites = 200, beta1.vec = c(2.5, 0.2, 0.5, 1, -1),
  ncov2 = 1, beta2.vec = rnorm(1, 0, 0)))
fm <- glm(dat$C ~ dat$rock + dat$oak + dat$chip1 + dat$chip2 + dat$Xrest, family = poisson)
plot(fm)


#   18.5.2 Bootstrapped and posterior predictive distributions
# ------------------------------------------------------------

# Create data set with major covariates plus one noise covariate (ncov2)
set.seed(18)
str(dat <- simDat18(nSites = 200, beta1.vec = c(2.5, 0.2, 0.5, 1, -1),
  ncov2 = 1, beta2.vec = rnorm(1, 0, 0.2)))

# Fit the mis-specified model without chip2
summary(fm <- glm(C ~ rock + oak + chip1, family = poisson, data = as.data.frame(dat)))

# Compute expected values and residual variation of data around them
mu.lambda <- predict(fm, type = 'response')            # Expected values
resi2.obs <- ((dat$C - mu.lambda) / sqrt(mu.lambda))^2 # Residual variation
(fit.obs <- sum(resi2.obs))          # Sum over data set
sum(residuals(fm, "pearson")^2)      # Same


# Simulate data conditional on observed covariate values, model and MLEs
simrep <- 100000                     # Number of bootstrapped
                                     # (replicate) data sets

# Create R objects ... to hold replicate (= bootstrapped) data sets
yrep <- array(NA, dim = c(dat$nSites, simrep))
# ... to hold squared Pearson residuals from replicated data
resi2.rep <- numeric(200)

# ...to hold sum of squared Pearson residuals
fit.rep <- numeric(simrep)

# Launch parametric bootstrap: produce simulated data sets
for(k in 1:simrep) {
  if(k %% 1000 == 0) cat(paste('\n*** iter', k))    # Counter
  Crep <- rpois(n = dat$nSites, lambda = mu.lambda) # Draw Poisson RVs
  yrep[, k] <- Crep                                 # ... and save them

  # Compute squared Pearson residuals and residual variation
  resi2.rep <- ((Crep - mu.lambda) / sqrt(mu.lambda))^2
  fit.rep[k] <- sum(resi2.rep)
}

# Fig. 18.4
par(mfrow = c(1,2), mar = c(6,5,4,2), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
hist(fit.rep, xlim = c(0, 1000), col = 'grey', xlab = 'Sum of squared Pearson residuals',
  main = 'Fit statistic for replicate data (grey)\n and for observed data (blue line)')
abline(v = fit.obs, lwd = 3, col = 'blue')

# Plot mean simulated count vs. observed counts for entire data set
plot(dat$C, apply(yrep, 1, mean), pch = 16, col = rgb(0,0,0,0.3), frame = FALSE,
  xlab = 'Observed counts', ylab = 'Mean of replicate counts',
  main = 'Comparison of expected vs. observed')
legend('bottomright', lwd = 3, lty = 3, col = c('red', 'blue'),
legend = c("1:1 line", "Regression of y on x"), bty = "n")
abline(0, 1, col = 'red', lwd = 3, lty = 3)
abline(lm(apply(yrep, 1, mean) ~ dat$C), col = 'blue', lwd = 3, lty = 3)

# Fig. 18.5
boxplot(t(yrep[1:100,]), outline = FALSE, xlab = 'Data number 1–100',
  frame = FALSE, main = 'Bootstrapped predictive distributions (boxes) and observed data (red)')
points(1:100, dat$C[1:100], pch = 16, cex = 0.8, col = 'red')

yobs <- dat$C                        # copy data
u <- numeric(dat$nSites)             # Vector to hold quantile residual
for(i in 1:dat$nSites){              # Loop over all data points
  u[i] <- mean((yrep[i,] < yobs[i]) + runif(1) * (yrep[i,] == yobs[i]))
}
par(mfrow = c(2,2), mar = c(6,5,4,2), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

# Fig. 18.6
hist(u, col = 'grey', main = 'Distribution of quantile residuals (u)')
plot(dat$oak, u, xlab = 'oak', pch = 16, frame = FALSE)
plot(dat$chip1, u, xlab = 'chip1', pch = 16, frame = FALSE)
plot(dat$chip2, u, xlab = 'chip2', pch = 16, frame = FALSE)


# Bayesian posterior predictive check

# Bundle and summarize data
str(dataList <- list(C = dat$C, rock = dat$rock, oak = dat$oak, chip = dat$chip1, n = length(dat$C)) )

# Write JAGS model file
cat(file = "model18.1.txt", "
  model {
  # Priors
  alpha ~ dunif(-10, 10)
  for(k in 1:3){
    beta[k] ~ dunif(-5, 5)
  }

  # Likelihood
  for (i in 1:n) { # Loop over all data points
    C[i] ~ dpois(lambda[i]) # The response variable (C above)
    lambda[i] <- exp(alpha + beta[1]* rock[i] + beta[2]* oak[i] + beta[3] * chip[i])
    resi2.obs[i] <- ((C[i] - lambda[i])/sqrt(lambda[i]))^2
    # Squared Pearson residuals for the observed data
  }
  # Create replicate data under the same model for each data point
  for (i in 1:n) {
    Crep[i] ~ dpois(lambda[i])
    resi2.rep[i] <- ((Crep[i] - lambda[i])/sqrt(lambda[i]))^2
    # Squared Pearson residuals for the replicate data
  }
  fit.obs <- sum(resi2.obs) # Sum over all data
  fit.rep <- sum(resi2.rep) # ditto
  }
")

# Function to generate starting values
inits <- function(){list(alpha = rnorm(1), beta = rnorm(3))}

# Parameters to estimate
params <- c("alpha", "beta", "fit.obs", "fit.rep", "Crep")

# MCMC settings
na <- 1000 ; ni <- 6000 ; nb <- 2000 ; nc <- 4 ; nt <- 4

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out18.1 <- jags(dataList, inits, params, "model18.1.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
traceplot(out18.1) # not shown
print(out18.1, 3)

par(mfrow = c(1,3), mar = c(6,5,4,2), cex.axis = 1.5, cex.lab = 1.5,
cex.main = 1.5) # Fig. 18.7
# Traditional plot for posterior predictive check (replicated vs. observed)
xylim <- c(0, 1000)
plot(out18.1$sims.list$fit.obs, out18.1$sims.list$fit.rep, xlim = xylim,
  ylim = xylim,xlab = 'Fit statistic observed data',
  ylab = 'Fit statistic replicate data', frame = FALSE)
abline(0, 1)
# Compute a Bayesian p-value for goodness-of-fit
mean(out18.1$sims.list$fit.rep > out18.1$sims.list$fit.obs) # [1] 0
# Non-traditional plot for posterior predictive check
hist(out18.1$sims.list$fit.rep, xlim = xylim, col = 'grey', xlab = 'Sum of squared Pearson residuals',
  main = 'Fit statistic for replicate data (grey)\n and observed data (blue)', freq = FALSE, breaks = 20)
hist(out18.1$sims.list$fit.obs, col = 'blue', freq = FALSE, add = TRUE, breaks = 20)

# Plot mean simulated count vs. observed counts for entire data set
yobs <- dat$C
plot(yobs, out18.1$mean$Crep, pch = 16, col = rgb(0,0,0,0.3), frame = FALSE,
  xlab = 'Observed counts', ylab = 'Mean of replicate counts',
  main = 'Comparison of expected vs. observed counts')
legend('bottomright', lwd = 3, lty = 3, col = c('red', 'blue'),
legend = c("1:1 line", "Regression of y on x"), bty = "n")
abline(0, 1, col = 'red', lwd = 3, lty = 3)
abline(lm(out18.1$mean$Crep ~ dat$C), col = 'blue', lwd = 3, lty = 3)

# Compute quantile residuals from Bayesian posterior predictive distributions and
# plot distribution (not shown)
yrepJ <- out18.1$sims.list$Crep # J for JAGS
u <- numeric(dat$nSites) # Vector to hold quantile residual
for(i in 1:dat$nSites){
  u[i] <- mean((yrepJ[,i]<yobs[i]) + runif(1) * (yrepJ[,i] == yobs[i]))
}
hist(u, col = 'grey', main = 'Distribution of quantile residuals (u)')


# Use of DHARMa for our own bootstrapped preditive simulations
library(DHARMa)
sim <- createDHARMa(simulatedResponse = yrep, observedResponse = yobs,
  fittedPredictedResponse = apply(yrep, 1, median))
par(mfrow = c(1,2)) # not shown
plotQQunif(sim)
plotResiduals(sim, rank = FALSE)
# Use of DHARMa for our own Bayesian posterior preditive simulations
sim <- createDHARMa(simulatedResponse = t(yrepJ), observedResponse = yobs,
  fittedPredictedResponse = apply(yrepJ, 2, median))
par(mfrow = c(1,2)) # not shown
plotQQunif(sim)
plotResiduals(sim, rank = FALSE)


#   18.5.3 Cross-validation for goodness-of-fit assessments
# ---------------------------------------------------------

set.seed(18)
str(dat <- simDat18(nSites = 100, beta1.vec = c(2.5, 0.2, 0.5, 1, -1),
ncov2 = 50, beta2.vec = rnorm(50, 0, 0)))
summary(fm <- glm(C ~ rock + oak + chip1, family = poisson, data = dat))

simrep <- 1000 # Number of samples of predictive distribution
YrepCV <- array(NA, dim = c(dat$nSite, simrep)) # Array to hold CV-replicated data
for(i in 1:dat$nSites){ # Loop over all data points = sites
  if(i %% 5 == 0) cat(paste('\n*** site', i))
  # Re-fit model to all data points minus 1 (i.e., minus site i)
  summary(fm_tmp <- glm(C[-i] ~ rock[-i] + oak[-i] + chip1[-i], family = poisson, data = dat))
  # Produce 'simrep' samples of the predictive distribution for the left-out datum at site i
  for(k in 1:simrep){ # Loop over simreps of predictive distribution
    lam_i <- exp(as.numeric(cbind(1, dat$rock[i], dat$oak[i], dat$chip1[i]) %*% coef(fm_tmp)))
    YrepCV[i, k] <- rpois(n = 1, lambda = lam_i)
  }
}
# Compute quantile residuals
yobs <- dat$C
u <- numeric(dat$nSites)
for(i in 1:dat$nSites){
  u[i] <- mean((YrepCV[i,] < yobs[i]) + runif(1) * (YrepCV[i,] == yobs[i]))
}

# Summarize
par(mfrow = c(1, 2), mar = c(6,5,4,2), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.2) # Plots not shown
# Plot mean simulated datum vs. observed datum for entire data set
meanYCV <- apply(YrepCV, 1, mean)
plot(dat$C, meanYCV, pch = 16, col = rgb(0,0,0,0.3), frame = FALSE,
  main = 'Mean simulated data vs. observed data')
abline(0, 1, col = 'red', lwd = 3, lty = 3)
abline(lm(meanYCV ~ dat$C), col = 'blue', lty = 2, lwd = 3)
# Plot frequency distribution of quantile residuals
hist(u, col = 'grey', main = 'Distribution of quantile residuals')


# Bayesian LOO-CV
# Select number of samples from the CV-based predictive distribution
simrep <- 2000 # Number of MCMC draws produced below
YrepCV <- array(NA, dim = c(dat$nSite, simrep)) # Array to hold CV-replicated data
# Function to generate starting values
inits <- function(){list(alpha = rnorm(1), beta = rnorm(3))}

# Parameters to estimate
params <- c("Crep")
# MCMC settings
na <- 100 ; ni <- 1200 ; nb <- 200 ; nc <- 4 ; nt <- 2
for(i in 1:dat$nSites){ # Loop over all data points = sites
  if(i %% 10 == 0) cat(paste('\n*** site', i))
  # Bundle and summarize data of length n-1 (each time with one datum used as external data !)
  dataList <- list(C = dat$C[-i], rock = dat$rock[-i],
    oak = dat$oak[-i], chip = dat$chip1[-i], n = length(dat$C)-1,
    pred.covs = cbind(dat$rock, dat$oak, dat$chip1)[i,])

  # Write JAGS model file
  cat(file = "model18.txt", "
    model {
    # Priors
    alpha ~ dunif(-10, 10)
    for(k in 1:3){
      beta[k] ~ dunif(-5, 5)
    }
    # Likelihood
    for (i in 1:n) { # Loop over all data points
      C[i] ~ dpois(lambda[i]) # The response variable (C above)
      lambda[i] <- exp(alpha + beta[1]* rock[i] + beta[2]* oak[i] + beta[3] * chip[i])
    }

    # Predict left-out (independent) datum
    lam_pred <- exp(alpha + beta[1]* pred.covs[1] + beta[2]* pred.covs[2] + beta[3] * pred.covs[3])
    Crep ~ dpois(lam_pred)
  }
  ")

  # Call JAGS (ART <1 min)
  out18 <- jags(dataList, inits, params, "model18.txt", n.iter = ni, n.burnin = nb,
    n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)

  # Save samples from CV-posterior predictive distribution
  YrepCV[i,] <- out18$sims.list$Crep
}

# Compute quantile residuals from Bayesian LOO-CV and produce Fig. 18.8
yobs <- dat$C
u <- numeric(dat$nSites) # Vector to hold quantile residual
for(i in 1:dat$nSites){
  u[i] <- mean((YrepCV[i,]<yobs[i]) + runif(1) * (YrepCV[i,]== yobs[i]))
}
hist(u, col = 'grey', main = 'Distribution of quantile residuals')


#   18.5.4 Identifying and correcting for overdispersion
# ------------------------------------------------------

# Create data set with major covariates plus one noise covariate (ncov2)
set.seed(18)
str(dat <- simDat18(nSites = 200, beta1.vec = c(2.5, 0.2, 0.5, 1, -1), ncov2 = 1,
  beta2.vec = rnorm(1, 0, 0.2)))

# Fit the mis-specified model without chip2
summary(fm <- glm(C ~ rock + oak + chip1, family = poisson, data = as.data.frame(dat)))

# Create data set with major covariates plus 45 noise covariates
set.seed(18)
str(dat <- simDat18(nSites = 200, beta1.vec = c(2.5, 0.2, 0.5, 1, -1), ncov2 = 45,
  beta2.vec = rnorm(45, 0, 0.1)))

# Fit the correctly specified model (in terms of the main covariates !)
data <- data.frame(C = dat$C, rock = dat$rock, oak = dat$oak, chip1 = dat$chip1, chip2 = dat$chip2)
summary(fm <- glm(C ~ rock + oak + chip1 + chip2, family = poisson, data = data))


#   18.5.5 Measures of the absolute fit of a model
# ------------------------------------------------
# (no code)


#   18.5.6 Parameter identifiability and robustness to assumptions
# ----------------------------------------------------------------

library(unmarked)

# Choose simreps and create arrays to hold true values and estimates
simrep <- 1000
esti <- true.vals <- array(NA, dim = c(simrep, 2), dimnames = list(NULL, c("lambda", "p")))
# Launch simulation to study parameter identifiability
system.time( # Set timer
  for(k in 1:simrep){ # Loop over simreps
    cat(paste("*** iter", k, "\n")) # Counter
    lam <- runif(1, 0, 10) # Pick a value for lambda within some range
    p <- runif(1, 0, 1) # Same for p
    true.vals[k,] <- c(lam, p) # Save true values for both
    N <- rpois(n = 100, lambda = lam) # Simulate latent abundances (N)
    C <- matrix(rbinom(n = 100, size = N, prob = p), ncol = 1) # Simulate observed counts (C)
    umf <- unmarkedFramePCount(y = C) # Bundle data for unmarked
    fm <- pcount(~1 ~1, data = umf, se = FALSE) # Fit model
    esti[k, ] <- c(exp(coef(fm)[1]), plogis(coef(fm)[2])) # Save estimates
  }
)

# Test for identifiability: compare estimates with truth
par(mfrow = c(1, 2)) # Plots not shown
plot(true.vals[,1], esti[,1], xlab = 'True', ylab = 'Estimated',
  main = 'Abundance (lambda)', frame = FALSE)
abline(0, 1, col = 'red', lwd = 2) # Red signifies truth
plot(true.vals[,2], esti[,2], xlab = 'True', ylab = 'Estimated',
  main = 'Detection (p)', frame = FALSE)
abline(0, 1, col = 'red', lwd = 2)


# 18.6 Model selection
# --------------------
# (no code)

#   18.6.1 When to do model selection ... and when not
# ----------------------------------------------------
# (no code)

#   18.6.2 Measuring the predictive performance of a model
# --------------------------------------------------------

set.seed(18)
beta2.vec <- rnorm(10, 0, 0.1)
trainDat <- simDat18(nSites = 50, beta1.vec = c(2, 0.2, 0.5, 1, -1),
  ncov2 = 10, beta2.vec = beta2.vec, show.plot = TRUE) # Training data
testDat <- simDat18(nSites = 50, beta1.vec = c(2, 0.2, 0.5, 1, -1),
  ncov2 = 10, beta2.vec = beta2.vec, show.plot = TRUE) # Testing data



#   18.6.3 Fully external model validation using “out-of-sample” data
# -------------------------------------------------------------------

# -----------------------------------------------------
# Fit 5 different models to the training data (using ML)
# ------------------------------------------------------
# The models differ in terms of their complexity. 
# They are only examples and we could in fact combine our 3 covariates 
# to an even larger set of models. But we do this here for 
# illustration only and don't want to show all possible subsets

# Note also that in this section and the next we work with maximum likelihood
# (or IRLS) rather than Bayes. Conceptually, there is not much difference, 
# other than having to deal with a whole distribution for each parameter
# and therefore also of lambda and of the predictive density of the data.

# Bundle data
dat <- data.frame(C = trainDat$C, rock = trainDat$rock, oak = trainDat$oak,
                  chip1 = trainDat$chip1, chip2 = trainDat$chip2)

# Fit 5 models of increasing complexity
fm1 <- glm(C ~ 1, family = poisson, data = dat)                          # Intercept only
fm2 <- glm(C ~ rock, family = poisson, data = dat)                       # + rock 
fm3 <- glm(C ~ rock + oak, family = poisson, data = dat)                 # + oak
fm4 <- glm(C ~ rock + oak + chip1, family = poisson, data = dat)         # + chip1
fm5 <- glm(C ~ rock + oak + chip1 + chip2, family = poisson, data = dat) # + chip2


# Now we want to assess the predictive performance of these 5 models
# (we could of course just peek at the AIC's, but we won't do so).
cbind(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5) )

# We compute the LPPD once for the training data and then for the testing data. 

# Compute the predictions for each data point
predlam1 <- predict(fm1, type = 'response')
predlam2 <- predict(fm2, type = 'response')
predlam3 <- predict(fm3, type = 'response')
predlam4 <- predict(fm4, type = 'response')
predlam5 <- predict(fm5, type = 'response')

# Model 'validation' with the same training set, i.e. "in-sample" or "is"
lpd1is <- sum(dpois(dat$C, predlam1, log = TRUE))
lpd2is <- sum(dpois(dat$C, predlam2, log = TRUE))
lpd3is <- sum(dpois(dat$C, predlam3, log = TRUE))
lpd4is <- sum(dpois(dat$C, predlam4, log = TRUE))
lpd5is <- sum(dpois(dat$C, predlam5, log = TRUE))


# Bundle the testing data ('OOS' stands for 'out-of-sample')
datoos <- data.frame(C = testDat$C, rock = testDat$rock, oak = testDat$oak,
                  chip1 = testDat$chip1, chip2 = testDat$chip2)

# Compute the predictions for each data point in the test data set
predlam1oos <- predict(fm1, type = 'response', newdata = datoos)
predlam2oos <- predict(fm2, type = 'response', newdata = datoos)
predlam3oos <- predict(fm3, type = 'response', newdata = datoos)
predlam4oos <- predict(fm4, type = 'response', newdata = datoos)
predlam5oos <- predict(fm5, type = 'response', newdata = datoos)

# True model validation with a completely independent test data set
lpd1oos <- sum(dpois(datoos$C, predlam1oos, log = TRUE))
lpd2oos <- sum(dpois(datoos$C, predlam2oos, log = TRUE))
lpd3oos <- sum(dpois(datoos$C, predlam3oos, log = TRUE))
lpd4oos <- sum(dpois(datoos$C, predlam4oos, log = TRUE))
lpd5oos <- sum(dpois(datoos$C, predlam5oos, log = TRUE))

# Compare the 'in-sample' lpd, the 'out-of-sample' lpd and the AIC
comp <- cbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
             c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  c(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5) )) 
dimnames(comp) <- list(c("Null", "+rock", "+oak", "+chip1", "+chip2"),
                    c("lpd_is", "lpd_oos", "AIC"))
print(comp, 4)

       # lpd_is lpd_oos   AIC
# Null   -201.8  -251.5 405.7
# +rock  -182.7  -262.7 369.4
# +oak   -179.5  -234.0 365.0
# +chip1 -162.9  -210.1 333.7
# +chip2 -113.3  -181.9 236.7

apply(comp[,1:2], 1, diff)
# > apply(comp[,1:2], 1, diff)
     # Null     +rock      +oak    +chip1    +chip2 
# -49.65371 -79.99393 -54.53418 -47.25606 -68.55099 

# Can see that in-sample lpd greatly overestimates predictive performance
# on a data set on which the parameters were not estimated from
# (and this is what the AIC attempts to correct for)


#   18.6.4 Approximating external validation for “in-sample” data with cross-validation
# -------------------------------------------------------------------------------------
# Create array for log-density
n <- 50               # Select sample size

logdensCV <- array(NA, dim = c(n, 5))   # number of data x number of models 

for(i in 1:n){                          # Loop over all data
  cat(paste("\n** Left-out datum", i, ""))  # Counter

  # Bundle the new training data set with one fewer cases
  datCV <- data.frame(C = trainDat$C[-i], rock = trainDat$rock[-i], oak = trainDat$oak[-i],
                  chip1 = trainDat$chip1[-i], chip2 = trainDat$chip2[-i])

  #datCV <- data.frame(C = trainDat$C, rock = trainDat$rock, oak = trainDat$oak,
  #                chip1 = trainDat$chip1, chip2 = trainDat$chip2)[-i,]   # Seems to work, too

  # Fit the five models to the 'CV-modified' training data set with n = 99
  fm1CV <- glm(C ~ 1, family = poisson, data = datCV)
  fm2CV <- glm(C ~ rock, family = poisson, data = datCV)
  fm3CV <- glm(C ~ rock + oak, family = poisson, data = datCV)
  fm4CV <- glm(C ~ rock + oak + chip1, family = poisson, data = datCV)
  fm5CV <- glm(C ~ rock + oak + chip1 + chip2, family = poisson, data = datCV)

  # Predict the single left-out datum and evaluate log-density for it
  loDat <- data.frame(C = trainDat$C[i], rock = trainDat$rock[i], oak = trainDat$oak[i],
                  chip1 = trainDat$chip1[i], chip2 = trainDat$chip2[i])
  predlam1CV <- predict(fm1CV, type = 'response', newdata = loDat)
  predlam2CV <- predict(fm2CV, type = 'response', newdata = loDat)
  predlam3CV <- predict(fm3CV, type = 'response', newdata = loDat)
  predlam4CV <- predict(fm4CV, type = 'response', newdata = loDat)
  predlam5CV <- predict(fm5CV, type = 'response', newdata = loDat)

  # True model validation with a completely independent test data set
  lpd1CV <- dpois(loDat$C, predlam1CV, log = TRUE)
  lpd2CV <- dpois(loDat$C, predlam2CV, log = TRUE)
  lpd3CV <- dpois(loDat$C, predlam3CV, log = TRUE)
  lpd4CV <- dpois(loDat$C, predlam4CV, log = TRUE)
  lpd5CV <- dpois(loDat$C, predlam5CV, log = TRUE)

  # Save log-densities for left-out datum
  logdensCV[i, ] <- c(lpd1CV, lpd2CV, lpd3CV, lpd4CV, lpd5CV)  
}

# Compute elppd by summing the elpd over all 100 points
lpd_loocv <- apply(logdensCV, 2, sum)


# Compare the 'in-sample' lpd, the 'out-of-sample' lpd, the lpd_loocv, and the AIC
comp <- cbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
              c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  lpd_loocv,
			  c(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5) )) 
dimnames(comp) <- list(c("Null", "+rock", "+oak", "+chip1", "+chip2"),
                    c("lpd_is", "lpd_oos", "lpd_loocv", "AIC"))
print(comp, 5)

        # lpd_is lpd_oos lpd_loocv    AIC
# Null   -201.83 -251.49   -208.00 405.67
# +rock  -182.68 -262.67   -190.25 369.36
# +oak   -179.51 -234.05   -191.96 365.02
# +chip1 -162.87 -210.13   -179.59 333.74
# +chip2 -113.33 -181.88   -128.81 236.67

# Same without the AIC
comp <- cbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
              c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  lpd_loocv) 
dimnames(comp) <- list(c("Null", "+rock", "+oak", "+chip1", "+chip2"),
                    c("lpd_is", "lpd_oos", "lpd_loocv"))
print(comp, 5)

        # lpd_is lpd_oos lpd_loocv
# Null   -201.83 -251.49   -208.00
# +rock  -182.68 -262.67   -190.25
# +oak   -179.51 -234.05   -191.96
# +chip1 -162.87 -210.13   -179.59
# +chip2 -113.33 -181.88   -128.81


par(mar = c(6, 4, 4, 2))
matplot(comp[,1:3], type = 'b', lwd = 3, lty = 1, frame = FALSE, axes = FALSE,
  pch = 16, cex = 2, col = c('black', 'red', 'blue'), xlab = '', ylab = '',
  main = 'log predictive densities (lpd)')
axis(1, at = 1:5, labels = c("Null", "+rock", "+oak", "+chip1", "+chip2"), las = 2)
axis(2, ylab = 'lppd')
legend('topleft', pch = 16, cex = 1.5, col = c('black', 'blue', 'red'),
  legend = c("in-sample naive", "in-sample LOO-CV", "true out-of-sample"), bty = "n")

# Note that the cross-validated estimate of ELPD is only slightly better here
# than the naive estimate that evaluates the log-density directly in-sample :(



# Next, we try 5-fold CV

# Create array for log-density
n <- 50
logdensCV <- array(NA, dim = c(n, 5))   # number of data x number of models 

# Create index for the 5 folds
fold <- rep(1:5, each = 10)

for(k in 1:5){                          # Loop over the k folds

  hold.out <- which(fold == k)

  # Bundle the new training data set without the data in one fold
  datCV <- data.frame(C = trainDat$C[-hold.out], rock = trainDat$rock[-hold.out],
            oak = trainDat$oak[-hold.out], chip1 = trainDat$chip1[-hold.out],
			chip2 = trainDat$chip2[-hold.out])

  #datCV <- data.frame(C = trainDat$C, rock = trainDat$rock, oak = trainDat$oak,
  #                chip1 = trainDat$chip1, chip2 = trainDat$chip2)[-i,]   # Seems to work, too

  # Fit the five models to the 'CV-modified' training data set with n = 80
  fm1CV <- glm(C ~ 1, family = poisson, data = datCV)
  fm2CV <- glm(C ~ rock, family = poisson, data = datCV)
  fm3CV <- glm(C ~ rock + oak, family = poisson, data = datCV)
  fm4CV <- glm(C ~ rock + oak + chip1, family = poisson, data = datCV)
  fm5CV <- glm(C ~ rock + oak + chip1 + chip2, family = poisson, data = datCV)

  # Predict the 10 left-out data and evaluate their log-density
  loDat <- data.frame(C = trainDat$C[hold.out], rock = trainDat$rock[hold.out],
     oak = trainDat$oak[hold.out], chip1 = trainDat$chip1[hold.out], chip2 = trainDat$chip2[hold.out])
  predlam1CV <- predict(fm1CV, type = 'response', newdata = loDat)
  predlam2CV <- predict(fm2CV, type = 'response', newdata = loDat)
  predlam3CV <- predict(fm3CV, type = 'response', newdata = loDat)
  predlam4CV <- predict(fm4CV, type = 'response', newdata = loDat)
  predlam5CV <- predict(fm5CV, type = 'response', newdata = loDat)

  # True model validation with a completely independent test data set
  lpd1CV <- dpois(loDat$C, predlam1CV, log = TRUE)
  lpd2CV <- dpois(loDat$C, predlam2CV, log = TRUE)
  lpd3CV <- dpois(loDat$C, predlam3CV, log = TRUE)
  lpd4CV <- dpois(loDat$C, predlam4CV, log = TRUE)
  lpd5CV <- dpois(loDat$C, predlam5CV, log = TRUE)

  # Save log-densities for data in left-out fold
  logdensCV[hold.out, ] <- cbind(lpd1CV, lpd2CV, lpd3CV, lpd4CV, lpd5CV)  
}

# Compute elppd by summing the elpd over all 50 points
lpd_5fold <- apply(logdensCV, 2, sum)

# Compare the 'in-sample' lpd, the 'out-of-sample' lpd, the lpd_loocv, lpd_5foldCV, and the AIC
comp <- cbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
              c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  lpd_loocv, lpd_5fold,
			  c(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5) )) 
dimnames(comp) <- list(c("Null", "+rock", "+oak", "+chip1", "+chip2"),
                    c("lpd_is", "lpd_oos", "lpd_loocv", "lpd_5fold", "AIC"))
print(comp, 5)

        # lpd_is lpd_oos lpd_loocv lpd_5fold    AIC
# Null   -201.83 -251.49   -208.00   -206.54 405.67
# +rock  -182.68 -262.67   -190.25   -192.61 369.36
# +oak   -179.51 -234.05   -191.96   -195.78 365.02
# +chip1 -162.87 -210.13   -179.59   -181.06 333.74
# +chip2 -113.33 -181.88   -128.81   -130.23 236.67

# Same without AIC
comp <- cbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
              c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  lpd_loocv, lpd_5fold) 
dimnames(comp) <- list(c("Null", "+rock", "+oak", "+chip1", "+chip2"),
                    c("lpd_is", "lpd_oos", "lpd_loocv", "lpd_5fold"))
print(comp, 5)


par(mar = c(6, 4, 4, 2))
matplot(comp[,1:4], type = 'b', lwd = 3, lty = 1, frame = FALSE, axes = FALSE,
  pch = 16, cex = 2, col = c('black', 'red', 'blue', 'green'), xlab = '', ylab = '',
  main = 'log pointwise predictive densities (lppd)')
axis(1, at = 1:5, labels = c("Null", "+rock", "+oak", "+chip1", "+chip2"), las = 2)
axis(2, ylab = 'lppd')
legend('topleft', pch = 16, cex = 1.5, col = c('black', 'green', 'blue', 'red'),
  legend = c("in-sample naive", "in-sample 5fold-CV", "in-sample LOO-CV", "true out-of-sample"), bty = "n")

# Note that the cross-validated estimate of ELPPD is only slightly better here
# than the naive estimate that evaluates the log-density directly in-sample :(


#   18.6.5 Approximating external validation for “in-sample” data with 
#          information criteria: AIC, DIC, WAIC, and LOO-IC
# --------------------------------------------------------------------

# Note that for this we partly move over to Bayes (for DIC and WAIC)
# However, we use our own likelihood optimisation software in R for the AIC

# ----------------------------
# (7.1) AIC and AICc (with ML)
# ----------------------------

# Define likelihood function
Xmat <- cbind(1, dat$rock, dat$oak, dat$chip1, dat$chip2)

# Define NLL for Poisson regression with intercept only (model 1)
NLL0 <- function(beta, C) {
  lambda <- exp(beta)
  LL <- dpois(C, lambda, log=TRUE)
  -sum(LL)
}

# Define NLL for Poisson regression with covariates (models 2-5)
NLL <- function(beta, C, Xmat) {
  lambda <- exp(Xmat %*% beta)
  LL <- dpois(C, lambda, log=TRUE)
  -sum(LL)
}

# Minimize NLL for data set to find MLEs, get SEs and 95% CIs
fit1 <- optim(1, NLL0, C = dat$C, hessian=TRUE, method = "BFGS")
get_MLE(fit1, 3)

fit2 <- optim(c(1, 0), NLL, C = dat$C, Xmat = Xmat[,1:2], hessian=TRUE, method = "BFGS")
get_MLE(fit2, 3)

fit3 <- optim(c(1, 0, 0), NLL, C = dat$C, Xmat = Xmat[,1:3], hessian=TRUE, method = "BFGS")
get_MLE(fit3, 3)

fit4 <- optim(c(1, 0, 0, 0), NLL, C = dat$C, Xmat = Xmat[,1:4], hessian=TRUE, method = "BFGS")
get_MLE(fit4, 3)

fit5 <- optim(c(1, 0, 0, 0, 0), NLL, C = dat$C, Xmat = Xmat[,1:5], hessian=TRUE, method = "BFGS")
get_MLE(fit5, 3)


# Compute AIC: deviance + 2 times number of params = 2 * NLL + 2 * np
aic1 <- 2 * fit1$value + 2 * length(fit1$par)
aic2 <- 2 * fit2$value + 2 * length(fit2$par)
aic3 <- 2 * fit3$value + 2 * length(fit3$par)
aic4 <- 2 * fit4$value + 2 * length(fit4$par)
aic5 <- 2 * fit5$value + 2 * length(fit5$par)
tt1 <- rbind(aic1, aic2, aic3, aic4, aic5)

# Compare with 'official' R solution
tt2 <- rbind(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5) )
comp <- data.frame('home-grown' = tt1, 'official' = tt2)
print(comp)

     # home.grown official
# aic1   405.6669 405.6669
# aic2   369.3568 369.3568
# aic3   365.0232 365.0232
# aic4   333.7411 333.7411
# aic5   236.6665 236.6665


# Or perhaps we need the small-sample version of the AIC ?
# Hurvich & Tsai (1989), and Burnham and Anderson (2002)
# Small-sample correction of the AIC
getAICc <- function(AIC, k, n = 50){
  # Get AICc from value of AIC, number of parameters k and sample size n
  AIC + (2*k*(k+1))/(n-k-1)
}
aicc1 <- getAICc(aic1, 1)
aicc2 <- getAICc(aic2, 2)
aicc3 <- getAICc(aic3, 3)
aicc4 <- getAICc(aic4, 4)
aicc5 <- getAICc(aic5, 5)

# Compare regular AIC and AICc
tt1 <- rbind(aic1, aic2, aic3, aic4, aic5)   # Regular AIC
tt2 <- rbind(aicc1, aicc2, aicc3, aicc4, aicc5)   # Corrected AICc
comp <- data.frame('regular AIC' = tt1, 'corrected AIC' = tt2)
print(comp, 2)

     # regular.AIC corrected.AIC
# aic1    405.6669      405.7502
# aic2    369.3568      369.6122
# aic3    365.0232      365.5449
# aic4    333.7411      334.6299
# aic5    236.6665      238.0301

# Not much difference really !


# Quick check with package AICcmodavg (see also Section 18.7)
library(AICcmodavg)
cand.set <- list(fm1, fm2, fm3, fm4, fm5)
modnames <- c('Null', '+rock', '+oak', '+chip1', '+chip2')
aictab(cand.set = cand.set, modnames = modnames)
# Looks good.


# ---------------------------------------------------------
# (7.2) DIC: a "somewhat Bayesian estimate of ELPD" (Bayes)
# ---------------------------------------------------------

# Use mean of the deviance as measure of model fit and apply a data-based penalty

# We use JAGS to fit the five models. Each time, we draw a sample from the
# posterior predictive distribution of the data (replicate counts 'Crep' below),
# which we will later also use with the R package loo.

# The DIC is calculated with the samples of the deviance of the model.
# This is automatically computed by JAGS, but for fun, we also 
# define our own deviance to better understand it.

# Bundle and summarize data for each analysis
str(dataList1 <- list(n = n, C = dat$C))  
str(dataList2 <- list(n = n, C = dat$C, X = Xmat[,1:2]))  
str(dataList3 <- list(n = n, C = dat$C, X = Xmat[,1:3]))  
str(dataList4 <- list(n = n, C = dat$C, X = Xmat[,1:4]))  
str(dataList5 <- list(n = n, C = dat$C, X = Xmat[,1:5]))  

# Write JAGS model file for model 1
cat(file="model1.txt", "
  model {
  # Priors
  beta ~ dnorm(0, 0.0001)

  # Likelihood, and computation of log-density and replicate data 'Crep'
  for (i in 1:n) {
    C[i] ~ dpois(exp(beta))
    logdens[i] <- logdensity.pois(C[i], exp(beta))
    Crep[i] ~ dpois(exp(beta))
  }
  # Define our own deviance for fun
  myDeviance <- -2 * sum(logdens)
  }
")

# Write JAGS model file for model 2
cat(file="model2.txt", "
  model {
  # Priors
  for(k in 1:2){
    beta[k] ~ dnorm(0, 0.0001)
  }
  # Likelihood, and computation of log-density and replicate data 'Crep'
  for (i in 1:n) {
    lambda[i] <- exp(inprod(X[i,], beta))
    C[i] ~ dpois(lambda[i])
    logdens[i] <- logdensity.pois(C[i], lambda[i])
    Crep[i] ~ dpois(lambda[i])
  }
  # Define our own deviance for fun
  myDeviance <- -2 * sum(logdens)
  }
")

# Write JAGS model file for model 3
cat(file="model3.txt", "
  model {
  # Priors
  for(k in 1:3){
    beta[k] ~ dnorm(0, 0.0001) # All params incl. intercept
  }
  # Likelihood, and computation of log-density and replicate data 'Crep'
  for (i in 1:n) {
    lambda[i] <- exp(inprod(X[i,], beta))
    C[i] ~ dpois(lambda[i])
    logdens[i] <- logdensity.pois(C[i], lambda[i])
    Crep[i] ~ dpois(lambda[i])
  }
  # Define our own deviance for fun
  myDeviance <- -2 * sum(logdens)
  }
")

# Write JAGS model file for model 4
cat(file="model4.txt", "
  model {
  # Priors
  for(k in 1:4){
    beta[k] ~ dnorm(0, 0.0001) # All params incl. intercept
  }
  # Likelihood, and computation of log-density and replicate data 'Crep'
  for (i in 1:n) {
    lambda[i] <- exp(inprod(X[i,], beta))
    C[i] ~ dpois(lambda[i])
    logdens[i] <- logdensity.pois(C[i], lambda[i])
    Crep[i] ~ dpois(lambda[i])
  }
  # Define our own deviance for fun
  myDeviance <- -2 * sum(logdens)
  }
")

# Write JAGS model file for model 5
cat(file="model5.txt", "
  model {
  # Priors
  for(k in 1:5){
    beta[k] ~ dnorm(0, 0.0001) # All params incl. intercept
  }
  # Likelihood, and computation of log-density and replicate data 'Crep'
  for (i in 1:n) {
    lambda[i] <- exp(inprod(X[i,], beta))
    C[i] ~ dpois(lambda[i])
    logdens[i] <- logdensity.pois(C[i], lambda[i])
    Crep[i] ~ dpois(lambda[i])
  }
  # Define our own deviance for fun
  myDeviance <- -2 * sum(logdens)
  }
")

# Function to generate starting values: we don't give any for these simple models
inits <- NULL

# Parameters to estimate
params <- c("beta", "logdens", "Crep", "myDeviance")

# MCMC settings
# Need rather long settings to get stable estimates of DIC !
na <- 1000  ;  ni <- 32000  ;  nb <- 2000  ; nc <- 4  ; nt <- 30

# Call JAGS, check convergence, summarize posteriors
out1 <- jags(dataList1, inits, params, "model1.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
out2 <- jags(dataList2, inits, params, "model2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
out3 <- jags(dataList3, inits, params, "model3.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
out4 <- jags(dataList4, inits, params, "model4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
out5 <- jags(dataList5, inits, params, "model5.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)

# Check convergence and print posterior summaries
# traceplot(out1)      # not shown
# traceplot(out2)      # not shown
# traceplot(out3)      # not shown
# traceplot(out4)      # not shown
# traceplot(out5)      # not shown
print(out1, 3)
print(out2, 3)
print(out3, 3)
print(out4, 3)
print(out5, 3)


# DIC
# What jagsUI returns for the 5 models
tt1 <- rbind(out1$DIC, out2$DIC, out3$DIC, out4$DIC, out5$DIC)

# "number of parameters" for the 5 models
tt2 <- rbind(out1$pD, out2$pD, out3$pD, out4$pD, out5$pD)
comp <- data.frame('DIC' = tt1, 'pD' = tt2)
print(comp)

       # DIC        pD
# 1 405.5821 0.9250468
# 2 370.7272 2.8652191
# 3 366.6263 4.1006406
# 4 337.4705 6.7054118
# 5 236.2562 4.7220331


# Now repeat that 'by hand'
# -------------------------

# Calculate pD manually using Gelman's approximation
# Don't forget deviance = -2 * log p (y | theta_hat)
dev1 <- out1$samples[,"deviance"]
dev2 <- out2$samples[,"deviance"]
dev3 <- out3$samples[,"deviance"]
dev4 <- out4$samples[,"deviance"]
dev5 <- out5$samples[,"deviance"]

# btw, in case you wonder what this is
str(dev1)


# Calculate the variance of deviance samples in each chain, divide by 2
# Have 4 chains !
pD1 <- sapply(dev1, var)/2
pD2 <- sapply(dev2, var)/2
pD3 <- sapply(dev3, var)/2
pD4 <- sapply(dev4, var)/2
pD5 <- sapply(dev5, var)/2

# And average them over all four chains
mean(pD1)
mean(pD2)
mean(pD3)
mean(pD4)
mean(pD5)


# Calculate the mean of deviance samples in each chain and add pD
# This is for 4 chains again
dic1 <- sapply(dev1, mean) + pD1
dic2 <- sapply(dev2, mean) + pD2
dic3 <- sapply(dev3, mean) + pD3
dic4 <- sapply(dev4, mean) + pD4
dic5 <- sapply(dev5, mean) + pD5

# and average them over chains
dic1_jagsUI <- mean(dic1)
dic2_jagsUI <- mean(dic2)
dic3_jagsUI <- mean(dic3)
dic4_jagsUI <- mean(dic4)
dic5_jagsUI <- mean(dic5)

# Output DIC and p_DIC
tt1 <- c(dic1_jagsUI, dic2_jagsUI, dic3_jagsUI, dic4_jagsUI, dic5_jagsUI)
tt2 <- c(mean(pD1), mean(pD2), mean(pD3), mean(pD4), mean(pD5))
data.frame('DIC' = tt1, 'p_DIC' = tt2)

       # DIC     p_DIC
# 1 405.5821 0.9250468
# 2 370.9609 3.0577671
# 3 366.5821 4.1006406
# 4 337.9955 6.7054118
# 5 236.7365 4.7220331


# The above calculations are based on R2WinBUGS/R2jags methods
# apparently this approximation of pD is not very good (even Gelman agrees with this)
# https://statmodeling.stat.columbia.edu/2006/07/06/number_of_param/
# Here's DIC calculated directly by rjags, which is supposed to be more accurate

# This does not work with parallel output so have to re-run models
out1.np <- jags(dataList1, inits, params, "model1.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)
dic1_rjags <- rjags::dic.samples(out1.np$model, n.iter = 1000, type = "pD")

out2.np <- jags(dataList2, inits, params, "model2.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)
dic2_rjags <- rjags::dic.samples(out2.np$model, n.iter = 1000, type = "pD")

out3.np <- jags(dataList3, inits, params, "model3.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)
dic3_rjags <- rjags::dic.samples(out3.np$model, n.iter = 1000, type = "pD")

out4.np <- jags(dataList4, inits, params, "model4.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)
dic4_rjags <- rjags::dic.samples(out4.np$model, n.iter = 1000, type = "pD")

out5.np <- jags(dataList5, inits, params, "model5.txt", n.iter = ni, n.burnin = nb,
  n.chains = nc, n.thin = nt, n.adapt = na, parallel = FALSE)
dic5_rjags <- rjags::dic.samples(out5.np$model, n.iter = 1000, type = "pD")

# Output DIC and p_DIC by Plummer (2002)
(dic1_rjags)
(dic2_rjags)
(dic3_rjags)
(dic4_rjags)
(dic5_rjags)

# Extract values of DIC and p_DIC from rjags version
DIC1_rjags <- sum(dic1_rjags[[1]] + dic1_rjags[[2]])
DIC2_rjags <- sum(dic2_rjags[[1]] + dic2_rjags[[2]])
DIC3_rjags <- sum(dic3_rjags[[1]] + dic3_rjags[[2]])
DIC4_rjags <- sum(dic4_rjags[[1]] + dic4_rjags[[2]])
DIC5_rjags <- sum(dic5_rjags[[1]] + dic5_rjags[[2]])

pD1_rjags <- sum(dic1_rjags[[2]])
pD2_rjags <- sum(dic2_rjags[[2]])
pD3_rjags <- sum(dic3_rjags[[2]])
pD4_rjags <- sum(dic4_rjags[[2]])
pD5_rjags <- sum(dic5_rjags[[2]])


# Compare the jagsUI and the rjags version of DIC and p_DIC
data.frame('DIC_jagsUI' = c(dic1_jagsUI, dic2_jagsUI, dic3_jagsUI, dic4_jagsUI, dic5_jagsUI),
           'pD_jagsUI' = c(mean(pD1), mean(pD2), mean(pD3), mean(pD4), mean(pD5)), 
           'DIC_rjags' = c(DIC1_rjags, DIC2_rjags, DIC3_rjags, DIC4_rjags, DIC5_rjags), 
           'pD_rjags' = c(pD1_rjags, pD2_rjags, pD3_rjags, pD4_rjags, pD5_rjags))

  # DIC_jagsUI pD_jagsUI DIC_rjags pD_rjags
# 1   405.5821 0.9250468  405.7537 1.044954
# 2   370.9609 3.0577671  369.7980 1.907052
# 3   366.5821 4.1006406  365.6748 3.046952
# 4   337.9955 6.7054118  335.0307 4.292559
# 5   236.7365 4.7220331  236.3742 4.844784

# We should probably recommend use of rjags::dic.samples vs what jagsUI reports



# --------------------------------------------------
# (7.4) WAIC: a "fully Bayesian estimate of ELPPD"
# --------------------------------------------------

## Extract posterior draws of the log-density of C at each site
ld.samples1 <- out1$sims.list$logdens
ld.samples2 <- out2$sims.list$logdens
ld.samples3 <- out3$sims.list$logdens
ld.samples4 <- out4$sims.list$logdens
ld.samples5 <- out5$sims.list$logdens

# what's the format of this ?
str(ld.samples1)


## Compute log-pointwise-predictive-density over data set
lppd1 <- sum(log(colMeans(exp(ld.samples1))))
lppd2 <- sum(log(colMeans(exp(ld.samples2))))
lppd3 <- sum(log(colMeans(exp(ld.samples3))))
lppd4 <- sum(log(colMeans(exp(ld.samples4))))
lppd5 <- sum(log(colMeans(exp(ld.samples5))))

## Compute penalty
pW1 <- sum(apply(ld.samples1, 2, var))
pW2 <- sum(apply(ld.samples2, 2, var))
pW3 <- sum(apply(ld.samples3, 2, var))
pW4 <- sum(apply(ld.samples4, 2, var))
pW5 <- sum(apply(ld.samples5, 2, var))

## Return WAIC
waic1 <- -2*(lppd1-pW1)
waic2 <- -2*(lppd2-pW2)
waic3 <- -2*(lppd3-pW3)
waic4 <- -2*(lppd4-pW4)
waic5 <- -2*(lppd5-pW5)

cbind(waic1, waic2, waic3, waic4, waic5)

       # waic1    waic2    waic3    waic4    waic5
# [1,] 410.844 374.8369 374.0893 348.5436 244.5751



# --------------------------------------------------
# (7.5) loo-IC by Vehtari et al. (2017)
# --------------------------------------------------

# Try loo-ic from loo package
library(loo)
(loo_ic1 <- loo(out1$sims.list$logdens))
(loo_ic2 <- loo(out2$sims.list$logdens))
(loo_ic3 <- loo(out3$sims.list$logdens))
(loo_ic4 <- loo(out4$sims.list$logdens))
(loo_ic5 <- loo(out5$sims.list$logdens))
LOOIC <- c(loo_ic1$est[3,1], loo_ic2$est[3,1], loo_ic3$est[3,1], loo_ic4$est[3,1], loo_ic5$est[3,1])




# --------------------------------------------------
# (7.6) WAIC from NIMBLE
# --------------------------------------------------

library(nimble)
library(MCMCvis)

# Bundle and summarize data
str(dataList <- list(n = n, C = dat$C, X = Xmat))  

# Write NIMBLE model file
model5 <- nimbleCode( {
# Priors
for(k in 1:5){
  beta[k] ~ dnorm(0, 0.0001) # All params incl. intercept
}
# Likelihood, and computation of log-density and replicate data yrep
for (i in 1:n) {
  lambda[i] <- exp(inprod(X[i,1:5], beta[1:5]))
  C[i] ~ dpois(lambda[i])
  logdens[i] <- dpois(C[i], lambda[i], log = TRUE)
  Crep[i] ~ dpois(lambda[i])
}
# Derived quantities: deviance
myDeviance <- -2 * sum(logdens[1:n])
} )


# Inits
inits <- function(){ list(beta = rnorm(5, 0, 1))}

# Parameters to estimate
params <- c("beta", "logdens", "Crep", "myDeviance")

# MCMC settings
na <- 1000  ;  ni <- 12000  ;  nb <- 2000  ; nc <- 4  ; nt <- 10

# Call NIMBLE (ART 25 sec), check convergence, and summarize posteriors
system.time( out5N <- 
    nimbleMCMC(code = model5,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE, WAIC = TRUE) )
	
par(mfrow = c(2,2), ask = TRUE)	
coda::traceplot(out5N$samples)                   # not shown
(nsum <- nimble_summary(out5N$samples, params))  # not shown

MCMCsummary(out5N$samples)
MCMCplot(out5N$samples)

# Extract WAIC
out5N$WAIC

# nimbleList object of type waicNimbleList
# Field "WAIC":
# [1] 244.3892
# Field "lppd":
# [1] -111.0352
# Field "pWAIC":
# [1] 11.15936

# Note that the calc goes like this
-2 * out5N$WAIC$lppd + 2 * out5N$WAIC$pWAIC

# [1] 244.3892


# -------------------------------------------------------------------------------
# (8) The Grande Comparison
# -------------------------------------------------------------------------------

# Compare the 'in-sample' lpd, the true 'out-of-sample' lpd, the lpd_loocv, lpd_5foldCV, 
#      and the AIC, DIC, WAIC and the LOO-IC
comp <- rbind(c(lpd1is, lpd2is, lpd3is, lpd4is, lpd5is),
              c(lpd1oos, lpd2oos, lpd3oos, lpd4oos, lpd5oos),
			  lpd_loocv,
			  lpd_5fold,
			  c(AIC(fm1), AIC(fm2), AIC(fm3), AIC(fm4), AIC(fm5)),
              c(out1$DIC, out2$DIC, out3$DIC, out4$DIC, out5$DIC),
              c(waic1, waic2, waic3, waic4, waic5),
              c(loo_ic1$est[3,1], loo_ic2$est[3,1], loo_ic3$est[3,1], loo_ic4$est[3,1], loo_ic5$est[3,1]))
dimnames(comp) <- list(c("lpd_is", "lpd_oos", "lpd_loocv", "lpd_5fold", "AIC", "DIC", "WAIC", "LOO-IC"),
                    c("Null", "+rock", "+oak", "+chip1", "+chip2"))
print(comp, 5)

             # Null   +rock    +oak  +chip1  +chip2
# lpd_is    -201.83 -182.68 -179.51 -162.87 -113.33
# lpd_oos   -251.49 -262.67 -234.05 -210.13 -181.88
# lpd_loocv -208.00 -190.25 -191.96 -179.59 -128.81
# lpd_5fold -206.54 -192.61 -195.78 -181.06 -130.23
# AIC        405.67  369.36  365.02  333.74  236.67
# DIC        405.61  370.94  366.60  337.18  236.86
# WAIC       410.84  374.84  374.09  348.54  244.58
# LOO-IC     410.51  374.95  373.97  349.02  244.99


# Scale the first 4 metrics in the same way as a deviance (i.e., multiply by -2)
comp2 <- comp
comp2[1:4,] <- -2 * comp2[1:4,]

print(comp2, 1)

          # Null +rock +oak +chip1 +chip2
# lpd_is     404   365  359    326    227
# lpd_oos    503   525  468    420    364
# lpd_loocv  416   381  384    359    258
# lpd_5fold  413   385  392    362    260
# AIC        406   369  365    334    237
# DIC        406   371  367    337    237
# WAIC       411   375  374    349    245
# LOO-IC     411   375  374    349    245



# 18.7 Model averaging
# --------------------

set.seed(18)
dat <- simDat18(nSites = 100, beta1.vec = c(0.5, 0.3, 0.2, 0.3, -0.4),
  ncov2 = 50, beta2.vec = rnorm(50, 0, 0.1), show.plot = TRUE)

# Fit 8 models: all combos of rock, oak and chip1
data <- data.frame(C = dat$C, rock = dat$rock, oak = dat$oak, chip1 = dat$chip1, chip2 = dat$chip2)
fm1 <- glm(C ~ 1, family = poisson, data = data)
fm2 <- glm(C ~ rock, family = poisson, data = data)
fm3 <- glm(C ~ oak, family = poisson, data = data)
fm4 <- glm(C ~ chip1, family = poisson, data = data)
fm5 <- glm(C ~ rock + oak, family = poisson, data = data)
fm6 <- glm(C ~ rock + chip1, family = poisson, data = data)
fm7 <- glm(C ~ oak + chip1, family = poisson, data = data)
fm8 <- glm(C ~ rock + oak + chip1, family = poisson, data = data)

summary(fm8) # Look at most complex model (not shown)

library(AICcmodavg)
cand.set <- list(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8)
modnames <- c('Null', 'rock', 'oak', 'chip1',
  'rock+ oak', 'rock+chip1', 'oak+ chip1', 'rock+oak+ chip1')
aictab(cand.set = cand.set, modnames = modnames)

newdata <- data.frame(rock = seq(-3, 3, length.out = 100), oak = 0, chip1 = 0)

# Predictions where we accommodate overdispersion
c.hat <- fm8$devi/fm8$df.resi # Our estimate of overdispersion
ma_pred <- modavgPred(cand.set = cand.set, modnames = modnames, newdata = newdata, c.hat = c.hat)
ma_pred # print out predictions (not shown)
plot(newdata$rock, ma_pred$mod.avg.pred, xlab = 'Rock', ylab = 'Expected count',
  type = 'l', lwd = 3, ylim = c(0, 9), frame = FALSE) # not shown
polygon(c(newdata$rock, rev(newdata$rock)), c(ma_pred$lower, rev(ma_pred$upper)),
  border = NA, col = rgb(0,0,0,0.1)) # not shown


# 18.8 Regularization: penalization, shrinkage, ridge, and lasso 
# ---------------------------------------------------------------

# -------------------------------------------------------
# Ridge-regression regularization with ticky rattlesnakes
# -------------------------------------------------------

# We regulate ALL coefficients
# We try out with 64 coefficients estimated from 100 data points


# Create one data set with 100 sites with 64 non-zero coefficients
# -----------------------------------------------------------------
set.seed(18)
str(dat <- simDat18(nSites = 100, beta1.vec = c(0.5, 0.4, 0.8, 1, -1), ncov2 = 60, 
  beta2.vec = rnorm(60, 0, 0.1)))
summary(fm <- glm(C ~ rock + oak + chip1 + chip2, family = poisson, data = as.data.frame(dat))) # works NOT
summary(fm <- glm(dat$C ~ dat$rock + dat$oak + dat$chip1 + dat$chip2, family = poisson))


#### Loop over multiple values of the regulator lam.pen
# --------------------------------------------------------------------------------------------

# Choose a range of values for the ridge-regression regulator lam_pen
n.regs <- 40             # Number of values to try out for regulator
round(lam.pen.vals <- exp(seq(0.001, 12, length.out = n.regs)),1)

 # [1]      1.0      1.4      1.9      2.5      3.4
 # [6]      4.7      6.3      8.6     11.7     16.0
# [11]     21.7     29.5     40.2     54.6     74.3
# [16]    101.1    137.5    187.0    254.4    346.1
# [21]    470.8    640.3    871.0   1184.8   1611.6
# [26]   2192.2   2982.0   4056.2   5517.4   7505.1
# [31]  10208.7  13886.4  18889.0  25693.7  34949.8
# [36]  47540.3  64666.6  87962.6 119650.9 162754.8


# Array to hold log lambda for left-out point 
# total number of points * number mcmc samples (see below) * number of values for regulator
lamCV0 <- array(NA, dim = c(dat$nSites, 2000, n.regs))   
lamCV100 <- array(NA, dim = c(dat$nSites, 2000, n.regs))   

# Array to hold sampled values of log predictive density of the left-out point 
# total number of points * number mcmc samples (see below) * number of values for regulator
LPPD0 <- array(NA, dim = c(dat$nSites, 2000, n.regs))   # No stabilization  
LPPD100 <- array(NA, dim = c(dat$nSites, 2000, n.regs)) # Stabilization at exp(+/- 100) 


# Launch LOO-CV and loop over all data points in the data set, repeat for every value for regulator
for(l in 1:n.regs){      # Loop over all regulator values
  # Select the value of the regulator 'lam.pen'
  lam.pen <- lam.pen.vals[l]   
  
  for(k in 1:dat$nSites){       # Loop over all data points in the entire data set
    cat(paste("\n *** Regulator value", round(lam.pen,2), "and leave-out point", k, "***\n"))

    # Extract response and explanatory variables for the left-out value
    lo.C <- dat$C[k]
    pred.covs <- cbind(dat$rock, dat$oak, dat$chip1, dat$chip2)[k,]
    pred.Xrest <- dat$Xrest[k, ]

    # Bundle and summarize data (in-sample data always with 1 fewer site)
	# and then we supply the left-out data for computation of log predictive density
    dataList <- list(C = dat$C[-k], rock = dat$rock[-k], oak = dat$oak[-k], 
      chip1 = dat$chip1[-k], chip2 = dat$chip2[-k], n = length(dat$C)-1,
      Xrest = dat$Xrest[-k,], lam.pen = lam.pen, ncov2 = dat$ncov2, 
      lo.C = lo.C, pred.covs = pred.covs, pred.Xrest = pred.Xrest)

    # Write JAGS model file
    cat(file="model18ridge.txt", "
    model {
    # Priors
    # No regularization for intercept
	alpha ~ dnorm(0, 0.001)

    # Regularization for ALL coefficients
    # lam.pen is the value of the regulator (called lambda in the literature)
    # Variance is 1/lam.pen, so precision = 1/variance = lam.pen
    # Primary coefficients
    for(k in 1:4){
      beta1[k] ~ dnorm(0, lam.pen)
    }
    # More 'nuisance params'
    for(k in 1:ncov2){
      beta2[k] ~ dnorm(0, lam.pen)
    }

    # Likelihood
    for (i in 1:n) {                 # Loop over all data points
      C[i] ~ dpois(lambda[i])        # The response variable (C above)
      lambda[i] <- exp(alpha + beta1[1]* rock[i] + beta1[2]* oak[i] + 
                     beta1[3] * chip1[i] + beta1[4] * chip2[i] + inprod(Xrest[i,], beta2))
    } 
    # Compute log predictive density of left-out count under the model and based on its covariates
 
    # Direct, without stabilization
    log.lam.pred0 <- alpha + beta1[1]* pred.covs[1] + beta1[2]* pred.covs[2] + beta1[3] * pred.covs[3] +
          beta1[4] * pred.covs[4] + inprod(pred.Xrest, beta2)
    lam.pred0 <- exp(log.lam.pred0)
	logdens.lo.C.0 <- logdensity.pois(lo.C, lam.pred0) # Without the stabilization
  
    # Stabilization at exp(+/- 100)
    log.lam.pred100 <- min(max(log.lam.pred0, -100), 100)
    lam.pred100 <- exp(log.lam.pred100)
	logdens.lo.C.100 <- logdensity.pois(lo.C, lam.pred100)
	} 
    ")

    # Function to generate starting values
    inits <- function(){list(alpha = rnorm(1), beta1 = rnorm(4),  beta2 = rnorm(dat$ncov2))}

    # Parameters to estimate
    params <- c("lam.pred0", "lam.pred100", "logdens.lo.C.0", "logdens.lo.C.100")

    # MCMC settings
    na <- 1000  ;  ni <- 5000  ;  nb <- 1000  ; nc <- 4  ; nt <- 8
    # na <- 3  ;  ni <- 3  ;  nb <- 1  ; nc <- 2  ; nt <- 1

    # Call JAGS (ART <1 min), check convergence and summarize posteriors
    out18ridge <- jags(dataList, inits, params, "model18ridge.txt", 
      n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
    print(out18ridge, 3)

    # Save samples of log-lambda and log-density for left-out datum
	# with different degrees of numerical stabilization
    lamCV0[k,,l] <- out18ridge$sims.list$lam.pred0
    lamCV100[k,,l] <- out18ridge$sims.list$lam.pred100

    LPPD0[k,,l] <- out18ridge$sims.list$logdens.lo.C.0
	LPPD100[k,,l] <- out18ridge$sims.list$logdens.lo.C.100
  }
}


#  Compare JAGS and R solutions
lamCV <- lamCV0 
LPPD <- LPPD0  

lamCV <- lamCV100
LPPD <- LPPD100

str(lamCV)
str(dat$C)

str(lppd_R <- array(NA, dim = dim(lamCV0)))

for(i in 1:100){
  for(k in 1:2000){
    lppd_R[i,k, ] <- dpois(dat$C[i], lamCV[i,k, ], log = TRUE)   
  }
}

# Sum test to compare JAGS and the R solution
sum(lppd_R)  ;  sum(LPPD)
### Gets heap of NAs in lppd_R
prod(dim(lppd_R))
sum(is.na(lppd_R))


# Plot and inspect LPPD profile to detect maximum
# Compute lppd for entire data set under the model
# Get posterior mean only
score <- numeric(n.regs)
for(l in 1:n.regs){
  score[l] <- sum(log(rowMeans(exp(LPPD100[,,l]), na.rm = TRUE)))
}
lam.pen.vals
lam.pen.vals[which(score == max(score, na.rm = TRUE))]
# [1] 29.52783

# Plot the E[lppd] vs the log(lam). The lambda with the highest E[lppd]
# is the optimal amount of regularization for prediction and thus the most parsimonious model.
plot(log(lam.pen.vals), score, pch = 16, col = rgb(0,0,0,0.3), cex = 2,
  type = 'b', frame = FALSE, main = 'Trace of ELPPD as a function of value of regulator')
abline(v = log(lam.pen.vals[which(score == max(score, na.rm = TRUE))]), lwd = 1, col = 'blue')

log(lam.pen.vals[which(score == max(score))])
lam.pen.vals[which(score == max(score))]

# > log(lam.pen.vals[which(score == max(score))])
# [1] 3.385333
# > lam.pen.vals[which(score == max(score))]
# [1] 29.52783

round(score,2)
# > round(score,2)
 # [1] -200.55 -194.73 -188.28 -182.70 -177.99 -172.70 -168.58 -164.55
 # [9] -162.04 -160.05 -159.08 -158.82 -159.65 -161.59 -164.38 -168.91
# [17] -174.92 -182.35 -190.68 -201.53 -213.86 -226.77 -241.15 -256.01
# [25] -270.02 -283.20 -293.75 -303.97 -311.76 -318.86 -324.55 -328.38
# [33] -332.45 -335.08 -337.36 -338.89 -340.14 -340.88 -341.53 -342.19


round(score-max(score),2)
# > round(score-max(score),2)
 # [1]  -41.72  -35.91  -29.46  -23.87  -19.16  -13.87   -9.75   -5.73
 # [9]   -3.22   -1.22   -0.25    0.00   -0.83   -2.77   -5.56  -10.09
# [17]  -16.10  -23.53  -31.86  -42.71  -55.04  -67.95  -82.32  -97.18
# [25] -111.19 -124.38 -134.93 -145.14 -152.93 -160.03 -165.72 -169.55
# [33] -173.62 -176.26 -178.54 -180.06 -181.31 -182.05 -182.71 -183.37


# Get full posterior
# Not sure what is the difference of this one to the one above
SCORE <- array(NA, dim = c(n.regs, dim(LPPD0)[2]))
for(l in 1:n.regs){
  for(s in 1:dim(LPPD0)[2]){
    SCORE[l, s] <- sum(LPPD0[,s,l])
  }
}
boxplot(t(SCORE), outline = FALSE, axes = FALSE)
axis(2)
axis(1, labels = round(log(lam.pen.vals), 2), at = 1:n.regs)
apply(SCORE, 1, mean)
pm <- apply(SCORE, 1, mean)
log(lam.pen.vals[which(pm == max(pm))])
lam.pen.vals[which(pm == max(pm))]
round(pm - max(pm))

# > log(lam.pen.vals[which(pm == max(pm))])
# [1] 4.923667
# > lam.pen.vals[which(pm == max(pm))]
# [1] 137.5059

 # [1] -29204 -12198  -4607  -1838  -1181   -583   -399   -226   -150
# [10]    -97    -61    -37    -23    -12     -4     -2      0     -2
# [19]     -6    -11    -18    -26    -34    -41    -50    -58    -66
# [28]    -72    -78    -83    -87    -90    -92    -94    -96    -97
# [37]    -98    -98    -99    -99




##### Now refit the models with all values of the regulators to all 100 data points and 
# monitor all parameter estimates


# Create array to hold posterior summaries for every fit
estimates <- array(NA, dim = c(66, 11, n.regs))
dimnames(estimates) <- list(rownames(out18ridge$summary), 
  colnames(out18ridge$summary), NULL)  # requires one fit of below

# Launch LOO-CV and loop over all data points in the data set, repeat for every value for regulator
for(l in 1:n.regs){      # Loop over all regulator values
  # Select the value of the regulator 'lam.pen'
  lam.pen <- lam.pen.vals[l]   
  cat(paste("\n *** Regulator value", round(lam.pen,2), "\n"))

  # Bundle and summarize data: this is ALL data now
  dataList <- list(C = dat$C, rock = dat$rock, oak = dat$oak, 
    chip1 = dat$chip1, chip2 = dat$chip2, n = length(dat$C),
    Xrest = dat$Xrest, lam.pen = lam.pen, ncov2 = dat$ncov2)

  # Write JAGS model file
  cat(file="model18ridgeX.txt", "
  model {
  # Priors
  alpha ~ dnorm(0, 0.001)

  # Regularization for ALL coefficients
  # lam.pen is the value of the regulator (called lambda in the literature)
  # Variance is 1/lam.pen, so precision = 1/variance = lam.pen
  # Primary parameters
  for(k in 1:4){
    beta1[k] ~ dnorm(0, lam.pen)
  }
  # More 'nuisance params'
  for(k in 1:ncov2){
    beta2[k] ~ dnorm(0, lam.pen)
  }

  # Likelihood
  for (i in 1:n) {                 # Loop over all data points
    C[i] ~ dpois(lambda[i])        # The response variable (C above)
    lambda[i] <- exp(alpha + beta1[1]* rock[i] + beta1[2]* oak[i] + 
                 beta1[3] * chip1[i] + beta1[4] * chip2[i] + inprod(Xrest[i,], beta2))
  } 
  } 
  ")

  # Function to generate starting values
  inits <- function(){list(alpha = rnorm(1), beta1 = rnorm(4),  beta2 = rnorm(dat$ncov2))}

  # Parameters to estimate
  params <- c("alpha", "beta1", "beta2")

  # MCMC settings
  na <- 1000  ;  ni <- 10000  ;  nb <- 2000  ; nc <- 4  ; nt <- 8
  #na <- 10  ;  ni <- 6  ;  nb <- 1  ; nc <- 2  ; nt <- 1

  # Call JAGS (ART <1 min), check convergence and summarize posteriors
  out18ridge <- jags(dataList, inits, params, "model18ridgeX.txt", 
    n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)
  #jagsUI::traceplot(out18ridge)
  print(out18ridge, 3)

  # Save posterior summary of fit
  estimates[,,l] <- out18ridge$summary
}



# Plot of parameter estimates (posterior means) vs. value of regulator
# 'Fan plot'

# Extract posterior means only
tmp <- estimates[,1,]
head(tmp)
matplot(log(lam.pen.vals), t(tmp[c(2:65),]), type = 'l', lty = 1, lwd = 3, col = 1:64,
  xlab = 'Value of regulator (log)', ylab = 'Coefficient', frame = FALSE)
abline(v = log(lam.pen.vals[which(score == max(score))]), lwd = 2, col = 'blue')
abline(h = 0, lwd = 2, lty = 2)

# Same with colors only for the 4 main parameters, but grey trajectories for the others
rrr <- range(tmp[c(2:5),])
matplot(log(lam.pen.vals), t(tmp[c(6:65),]), type = 'l', lty = 1, lwd = 2, col = rgb(0,0,0,0.3),
 xlab = 'Value of regulator (log)', ylab = 'Coefficient', ylim = rrr, frame = FALSE)
matpoints(log(lam.pen.vals), t(tmp[c(2:5),]), type = 'l', lty = 1, lwd = 2, col = 1:4)
abline(h = 0, lwd = 2, lty = 2)
abline(v = log(lam.pen.vals[which(score == max(score))]), lwd = 1, col = 'blue')


# ELPPD trajectory and fan plot
# -----------------------------
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# Plot the E[lppd] vs the log(lam). The lambda with the highest E[lppd]
# is the optimal amount of regularization for prediction and thus the most parsimonious model.
plot(log(lam.pen.vals), score, pch = 16, col = rgb(0,0,0,0.3), cex = 2,
  type = 'b', frame = FALSE, xlab = 'Value of regulator (log)', ylab = "log score", main = 'Predictive accuracy')
abline(v = log(lam.pen.vals[which(score == max(score))]), lwd = 2, col = 'blue')

matplot(log(lam.pen.vals), t(tmp[c(2:25),]), type = 'l', lty = 1, lwd = 3, col = 1:24,
  xlab = 'Value of regulator (log)', ylab = 'Coefficient', main = "Parameter estimates", frame = FALSE)
abline(v = log(lam.pen.vals[which(score == max(score))]), lwd = 2, col = 'blue')
abline(h = 0, lwd = 2, lty = 2)


# 18.9 Summary and outlook
# ------------------------
# (no code)
