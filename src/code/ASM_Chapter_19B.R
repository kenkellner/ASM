# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -------------------------------------------------------
# Chapter 19B  --  Binomial N-mixture models
# -------------------------------------------------------

# Last changes: 6 July 2024


# 19B.1 Introduction
# -----------------
# (no code)


# 19B.2 Data generation 
# --------------------

# Load Swiss MHB elevation data for 267 quadrats
library(unmarked)
data(crossbill)
str(crossbill)           # not shown

# Pull out elevation data, jitter some and sort
set.seed(191)
mhbElev <- sort(jitter(crossbill$ele, factor = 2))
summary(mhbElev)

# Histogram of elevation (not shown)
par(mar = c(5,5,5,3), cex.lab = 1.5, cex.axis = 1.5)
hist(mhbElev, breaks = 50, xlim = c(0, 3000), 
  main = 'Mean elevation of 267 MHB quadrats in Switzerland')

# Load the full Swiss landscape data and grab the elevation data
data(Switzerland)
ch <- Switzerland
str(ch)
chElev <- ch$elevation # Swiss elevation data in km units
summary(chElev)

# Make a plot of elevation (Fig. 19B.2)
library(raster)
par(mfrow = c(1, 1), mar = c(3,4,4,8), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = chElev))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Swiss elevation map (in m)", zlim = c(0, 4500))

# Simulate MHB "look-alike data" for the bullfinch
set.seed(191)

# Pick design constants
nSites <- 267           # Number of quadrats, or sites
nVisits <- 3            # Number of occasions, or visits per site

# Create scaled version of MHB elevation covariate with units of 1 km
mhbElevScaled <- mhbElev/1000

# Pick values for the regression parameters in expected abundance
alpha.lam <- -3           # Intercept
beta1.lam <- 8.5          # Linear effect of elevation
beta2.lam <- -3.5         # Quadratic effect of elevation

# Compute expected abundance (lambda) for MHB quadrats
lambda <- exp(alpha.lam + beta1.lam * mhbElevScaled + beta2.lam * mhbElevScaled^2) 

# Compute realized abundance (N) for MHB quadrats
N <- rpois(n = nSites, lambda = lambda)
table(N)              # Distribution of abundances across sites
sum(N > 0) / nSites   # Empirical occupancy probability

# Figure 19B.3 left
# State process (true states)
par(mfrow = c(1, 3), mar = c(6,6,5,3), cex.lab = 2, cex.axis = 2, cex.main = 2)
plot(mhbElev, N, main = "Bullfinch abundance", pch = 16, cex = 2,
  col = rgb(0,0,0,0.4), frame = FALSE, xlim = c(0, 3000),
  ylim = c(0, 15), xlab = 'Elevation (m)', ylab = 'Abundance per 1km2') 
lines(mhbElev, lambda, lwd = 7, col = rgb(1,0,0,0.4))

# Compute optimal elevation for lambda (in kilometric units !)
(opt.elev.true <- - beta1.lam / ( 2 * beta2.lam) )

# Pick values for the 2 regression parameters in detection and get p
alpha.p <- 2
beta.p <- -2
p <- plogis(alpha.p + beta.p * mhbElevScaled)

# Save true parameters into vector
truth <- c(alpha.lam=alpha.lam, beta1.lam=beta1.lam, 
  beta2.lam=beta2.lam, alpha.p=alpha.p, beta.p=beta.p)

# Figure 19B.3 middle
# Observation process (characterized by detection probability)
plot(mhbElev, p, xlab = "Elevation (m)", 
  ylab = " Detection probability (p)", type = 'l', frame = FALSE, 
  xlim = c(0, 3000), lwd = 7, col = rgb(1,0,0,0.4), 
  main = "Bullfinch detectability per survey", ylim = c(0, 1))

# Simulate the observation process
C <- array(dim = c(nSites, nVisits))
for(j in 1:nVisits){
  C[,j] <- rbinom(n = nSites, size = N, prob = p)
}

# Figure 19B.3 right
# Observed counts = 'relative abundance' = 'abundance index'
matplot(mhbElev, C, ylim = c(0, 15), xlab = "Elevation (m)", las = 1,
  ylab = "Counts (C)", main = "Observed bullfinch counts", pch = 16, 
  cex = 2, col = rgb(0,0,0,0.4), frame = FALSE, xlim = c(0, 3000))
lines(mhbElev, lambda, type = "l", col = rgb(1,0,0,0.5), lwd = 7)
lines(mhbElev, lambda * p, type = "l", col = "black", lwd = 5)

# Compare true abundance with the replicated counts
cbind('True state' = N, 'Obs Visit 1' = C[,1], 'Obs Visit 2' = C[,2], 'Obs Visit 3' = C[,3])# Look true state and at three observed states at each site

sum(apply(C, 1, sum) > 0) # Apparent distribution (prop. occupied sites)
sum(N > 0)                # True occupancy

library(raster)

# Compute true expected abundance (lambda) for all Swiss quadrats
chElevScaled <- chElev / 1000
lambdaCH <- exp(alpha.lam + beta1.lam * chElevScaled + beta2.lam * chElevScaled^2) 

# Inspect quadrat-level expected abundance
summary(lambdaCH)
hist(lambdaCH)      # not shown

# Plot true Swiss bullfinch abundance map (Fig. 19B.4)
par(mfrow = c(1, 1), mar = c(3,4,4,6), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = lambdaCH))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "True Swiss SDM for the bullfinch (expected abundance)", zlim = c(0, 9))

# Add up to get Swiss total population size
( Ntot_true <- sum(lambdaCH) )


# 19B.3 A detection-naive analysis of the maximum count per site 
# --------------------------------------------------------------

maxC <- apply(C, 1, max)
fm <- glm(maxC ~ mhbElevScaled + I(mhbElevScaled^2), family = poisson)
summary(fm) # not shown
pred.elev <- seq(0.2, 2.5, length.out = 1000)  # Prediction covariate
lpred <- predict(fm, type = 'link', se = TRUE, 
    newdata = data.frame(mhbElevScaled = pred.elev))
pred <- exp(lpred$fit)
LCL <- exp(lpred$fit-2*lpred$se)
UCL <- exp(lpred$fit+2*lpred$se)

sum(maxC)  # Detection-naïve estimate of total N in 267 sample quadrats
sum(N)     # True total N in 267 sample quadrats

# Fig. 19B.5
par(mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(1000*pred.elev, pred, type = 'l', col = rgb(0,0,1, 0.4), lwd = 3,
  xlab = "Elevation (m)", ylab = "Number", xlim = c(0, 3000), ylim = c(0, 10), main = "Confounding of state and observation processes", frame = FALSE)
polygon(c(1000*pred.elev, rev(1000*pred.elev)), c(LCL, rev(UCL)), col = rgb(0,0,1, 0.2), border = NA)
points(mhbElev, maxC, col = rgb(0,0,0, 0.4), pch = 16, cex = 1.3)
lines(mhbElev, lambda, type = "l", col = rgb(1,0,0,0.5), lwd = 5)

(opt.elev2 <- -fm$coef[2] / ( 2 * fm$coef[3]))
(opt.elev2 - opt.elev.true) # Difference to true opt.elev in simulation
abline(v = 1000*opt.elev2, col = rgb(0,0,1,0.4), lwd = 3) # Observed
abline(v = 1000*opt.elev.true, col = rgb(1,0,0,0.4), lwd = 3) # True

# Make predictions for all Switzerland: prototypical SDM (Fig. 19B.6)
pred.lam_obs <- predict(fm, type = "response", 
  newdata = data.frame(mhbElevScaled = chElev/1000))
summary(pred.lam_obs)

par(mfrow = c(1, 1), mar = c(3,4,4,6), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pred.lam_obs))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Detection-naive Swiss SDM for the bullfinch (expected max counts)", zlim = c(0, 6))

# Estimated national "relative population size"
( SwissTotalMaxC <- sum(pred.lam_obs) )

# load required packages
library(ASMbook); library(jagsUI); library(nimble); library(rstan); library(TMB)


# 19B.4 Likelihood analysis with canned functions in the R package unmarked 
# -------------------------------------------------------------------------
# (no code)


#   19B.4.1 Fitting the model to the simulated bullfinch data set
# -------------------------------------------------------------

# Load unmarked, format data and summarize
library(unmarked)
summary( umf <- unmarkedFramePCount(y = C, 
  siteCovs = data.frame(elev=mhbElevScaled, elev2=mhbElevScaled ^2)) )

out19B.4 <- pcount(~ elev ~ elev + elev2, data = umf)

# Check for sensitivity of solutions to choice of K
fhm11 <- pcount(~elev ~ elev+elev2, data=umf, K = 11)
fhm13 <- pcount(~elev ~ elev+elev2, data=umf, K = 13)
fhm20 <- pcount(~elev ~ elev+elev2, data=umf, K = 20)
fhm50 <- pcount(~elev ~ elev+elev2, data=umf, K = 50)
fhDef <- pcount(~elev ~ elev+elev2, data=umf)  # With default K = 108
fhm200 <- pcount(~elev ~ elev+elev2, data=umf, K = 200)

# Show only the part where something happens ...
aic <- c(fhm20@AIC, fhm50@AIC, fhDef@AIC, fhm200@AIC)
k <- c(20, 50, 108, 200)
loglik <- sapply(list(fhm20, fhm50, fhDef, fhm200), logLik)

# ... or the whole range of values tried out for K
aic <- c(fhm11@AIC, fhm13@AIC, fhm20@AIC, fhm50@AIC, fhDef@AIC, fhm200@AIC)
k <- c(11, 13, 20, 50, 108, 200)
loglik <- sapply(list(fhm11, fhm13, fhm20, fhm50, fhDef, fhm200), logLik)

# Plot both AIC and log-likelihood (Fig. 19B.7)
par(mfrow=c(2,1), mar = c(6,6,4,2), cex.lab = 1.5, cex.axis = 1.5)
plot(k, aic, type='l', ylab="AIC", xlab="Value of K", xlim=c(9, 200), frame = FALSE, ylim = c(1785, 1800))
points(k, aic, col='red', pch=16, cex = 2)

plot(k, loglik, type='l', ylab="log-likelihood", xlab="Value of K",
     xlim=c(9, 200), frame = FALSE, ylim = c(-895, -888))
points(k, loglik, col='blue', pch=16, cex = 2)

print(out19B.4)
unm_est <- coef(out19B.4)      # Save estimates

# Make predictions of abundance and detection(Fig. 19B.8)
state.pred <- predict(out19B.4, type = 'state')
det.pred <- predict(out19B.4, type = 'det')
p.pred <- matrix(det.pred[,1], nrow = nSites, byrow = TRUE) # reformat
p.LCL <- matrix(det.pred[,3], nrow = nSites, byrow = TRUE)  # reformat
p.UCL <- matrix(det.pred[,4], nrow = nSites, byrow = TRUE)  # reformat

# Plot predicted elevation profiles for state and observation processes
par(mfrow = c(1, 2), mar = c(5,5,4,2), cex.lab = 1.5, cex.axis = 1.5)
# Expected abundance (lambda)
plot(mhbElev, state.pred[,1], xlab = 'Elevation (m)', ylab = 'lambda', frame = FALSE, col = 'blue', lwd = 3, main = 'State process', type = 'l', ylim = c(0, 15))
lines(mhbElev, lambda, lwd = 3, col = 'red')
polygon(c(mhbElev, rev(mhbElev)), c(state.pred[,3], rev(state.pred[,4])), col = rgb(0,0,1, 0.2), border = NA)
legend('topleft', legend = c('Estimate', 'Truth'), lwd = 3, col = c('blue', 'red'), bty = 'n', cex = 1.3)

# Detection
plot(mhbElev, p.pred[,1], xlab = 'Elevation (m)', ylab = 'Detection prob.', frame = FALSE, col = 'blue', lwd = 3, main = 'Observation process', type = 'l', ylim = c(0, 1))
lines(mhbElev, p, lwd = 3, col = 'red')
polygon(c(mhbElev, rev(mhbElev)), c(p.LCL[,1], rev(p.UCL[,1])), col = rgb(0,0,1, 0.2), border = NA)

# Compare parameter estimates with truth
comp <- cbind(truth=truth, unmarked=unm_est)
print(comp, 4)

# Compare estimated total population size in surveyed sample to truth
true_Ntotal <- sum(N)
naive_Ntotal <- sum(apply(C, 1, max))
unm_Ntotal <- sum(bup(ranef(out19B.4)))
comp <- c(truth=true_Ntotal, naive=naive_Ntotal, unmarked=unm_Ntotal)
print(comp, 1)


#   19B.4.2 Spatial prediction of expected abundance 
# --------------------------------------------------------------

# Make predictions for all Switzerland
pred.lam <- predict(out19B.4, type = "state", 
  newdata = data.frame(elev = chElevScaled, elev2 = chElevScaled^2))
str(pred.lam)
head(pred.lam)
summary(pred.lam)

# The Holy Grail (Fig. 19B.9)
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))

par(mfrow = c(2, 2), mar = c(3,4,4,6), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pred.lam[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Swiss SDM for the bullfinch\n(estimated true abundance, lambda)", zlim = c(0, 10))
r2 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pred.lam[,2]))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Prediction SE", zlim = c(0, 2))
r3 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pred.lam[,3]))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Uncertainty (Lower limit of 95% CI)", zlim = c(0, 14))
r4 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pred.lam[,4]))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Uncertainty (Upper of 95% CI)", zlim = c(0, 14))

# Estimated national population size
( SwissTotalN <- sum(pred.lam[,1]) )

# Compare with true national population size
Ntot_true


#   19B.4.3 Computing SEs for a derived quantity using the detal method
#           and the bootstrap
# ---------------------------------------------------------------------

# Point estimate of optimal elevation (in metres; unmarked solution)
opt.elev.um <- -coef(out19B.4)[2] / (2 * coef(out19B.4)[3])
opt.elev.um

# Get the ingredients for the beast
out19B.4                   # Output from pcount()
str(out19B.4@opt)          # Internal output from optim()
( VC <- solve(out19B.4@opt$hessian) ) # Get VC matrix from Hessian

( b1_hat <- coef(out19B.4)[2]) # Point estimate of elev linear
( b2_hat <- coef(out19B.4)[3]) # ... of elev squared
( var_b1 <- VC[2,2])       # Variance of elev linear
( var_b2 <- VC[3,3])       # ... of elev squared
( cov_b1_b2 <- VC[2,3])    # Covariance between elev linear and squared

# Get delta-method approximation of SE of optimum elevation
var_elev.opt <- (-b1_hat / (2*b2_hat))^2 * 
                (((2^2 * var_b2) / (2^2 * b2_hat^2)) +
                (var_b1 / b1_hat^2 ) - 
                ((2 * cov_b1_b2) / (b1_hat * b2_hat) ) )
sqrt(var_elev.opt)                # SE is square root of variance

# Get parametric bootstrap approximation of SE of optimum elevation
library(MASS)              # Load MASS for mvrnorm()
?mvrnorm                   # Check syntax

# Get mean vector and vc matrix for beta1 and beta2
( mu <- c(b1_hat, b2_hat)) # Collect means
( vc <- VC[2:3, 2:3])      # Get relevant part of VC matrix from above

# Draw large number of values from this multivariate normal distribution
boot_beta <- mvrnorm(n = 100000, mu = mu, Sigma = vc)

# Compute the function that gives optimum elevation
opt.elev_boot <- -boot_beta[,1] / (2 * boot_beta[,2])

# Plot bivariate sampling distribution of beta1 and beta2 and
# sampling distribution of -b1 / (2 x b2)     # Fig. 19B.10
par(mfrow = c(1, 2), mar = c(5,5,5,6), cex.main = 1.5, cex.lab = 1.5,
  cex.axis = 1.5)
plot(boot_beta, asp = 1, xlab = 'beta1', ylab = 'beta2', frame = FALSE,
  pch = 16, col = rgb(0, 0, 0, 0.2), main = 'Bivariate sampling distribution\n of beta1 and beta2')
hist(1000*opt.elev_boot, breaks = 50, main = 'Sampling distribution\nof optimum elevation', xlab = 'Elevation (m)', xlim = c(1000, 1600))

# Get parametric-bootstrapped SE and CI of optimum elevation
# SE and CIs are SD and 0.025/0.975th percentiles
SE_opt.elev <- sd(opt.elev_boot)
( CI_opt.elev <- quantile(opt.elev_boot, c(0.025, 0.975)))

# Compare two non-Bayesian solutions to the SE of optimum elevation
comp <- c(sqrt(var_elev.opt), SE_opt.elev)
names(comp) <- c("Delta method", "Parametric Bootstrap")
print(comp, 4)

# 19B.5 Bayesian analysis with JAGS 
# ---------------------------------

# Bundle and summarize data
elev2 <- mhbElevScaled^2         # squared elevation
str(dataList <- list(C=C, elev = mhbElevScaled, elev2 = elev2, 
  nSites = nSites, nVisits = nVisits) )

# Write JAGS model file
cat(file="model19B.5.txt", "
model {
# Priors
mean.lam ~ dunif(0, 1)
alpha.lam <- log(mean.lam)
beta1.lam ~ dnorm(0, 0.001)
beta2.lam ~ dnorm(0, 0.001)
mean.p ~ dunif(0, 1)
alpha.p <- logit(mean.p)
beta.p ~ dnorm(0, 0.001)

# Likelihood
# Biological model for true abundance
for (i in 1:nSites) {             # Loop over sites
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha.lam + beta1.lam * elev[i] + beta2.lam * elev2[i]
}

# Observation model for replicated counts
for (i in 1:nSites) {             # Loop over sites
  for (t in 1:nVisits) {          # Loop over all n observations
    C[i,t] ~ dbin(p[i,t], N[i]) 
    logit(p[i,t]) <- alpha.p + beta.p * elev[i]
  }
}

# Derived quantities
totalN <- sum(N[])   # Estimate total population size across all sites
opt.elev <- -beta1.lam / (2 * beta2.lam)  # Optimum elevation

# Assess model fit using Chi-squared discrepancy
for (i in 1:nSites) {             # Loop over sites
  for (t in 1:nVisits) {          # Loop over all n observations
    # Compute fit statistic for observed data
    eval[i,t] <- p[i,t]*N[i]
    E[i,t] <- pow((C[i,t] - eval[i,t]),2) / (eval[i,t] + 0.01)# Avoid div. by 0 !
    # Generate replicate data and compute fit stats for them
    C.new[i,t] ~ dbin(p[i,t], N[i])
    E.new[i,t] <- pow((C.new[i,t] - eval[i,t]),2) / (eval[i,t] + 0.01)
  }
}
fit <- sum(E[,])
fit.new <- sum(E.new[,])
}
")

# Function to generate starting values
Nst <- apply(C, 1, max) + 1
inits <- function(){list(N = Nst, mean.lam=runif(1, 0, 1), beta1.lam=rnorm(1, 0, 1), 
                         beta2.lam=rnorm(1, 0, 1), mean.p=runif(1), beta.p=rnorm(1, 0, 1))}

# Parameters to estimate
params <- c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p", 
            "mean.lam", "mean.p", "totalN", "fit", "fit.new", "opt.elev", "N")

# MCMC settings
na <- 10000;  ni <- 150000;  nb <- 50000; nc <- 4; nt <- 100 
# na <- 50;  ni <- 1100;  nb <- 100; nc <- 4; nt <- 2        # test

# Call JAGS (ART 11 min)
out19B.5 <- jags(dataList, inits, params, "model19B.5.txt", n.iter = ni,
  n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na,
  parallel = TRUE)

# Traceplots (Fig. 19B.11)
jagsUI::traceplot(out19B.5)      # All parameters
jagsUI::traceplot(out19B.5, c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p", "opt.elev", "totalN"), layout = c(3, 3)) # Key parameters only

# Check out the Rhat values
which(out19B.5$summary[,8] > 1.1)          # anybody not converged?
out19B.5$summary[which(out19B.5$summary[,8] > 1.1) ,8] # how bad?

# Look at traceplot for parameter that has not converged
# jagsUI::traceplot(out19B.5, c("N[2]"), layout = c(2, 2))

# Compute Bayesian p-value and print
( bpv <- mean(out19B.5$sims.list$fit.new > out19B.5$sims.list$fit) )

# Fig. 19B.12
par(mar = c(6,6,5,4), cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3)
plot(out19B.5$sims.list$fit, out19B.5$sims.list$fit.new, 
  main = paste("Bayesian p-value:", round(bpv, 4)),
  xlab = "Discrepancy for actual data set", 
  ylab = "Discrepancy for perfect data sets", 
  pch= 16, cex = 1.2, col = rgb(0,0,0,0.2), frame = FALSE)
abline(0,1, lwd = 2, col = "black")

print(out19B.5, 3)

par(mfrow = c(3,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
hist(out19B.5$sims.list$alpha.lam, col = "grey", main = "alpha.lam", xlab = "", breaks = 50)
abline(v = alpha.lam, lwd = 3, col = "red")
hist(out19B.5$sims.list$beta1.lam, col = "grey", main = "beta1.lam", xlab = "", breaks = 50)
abline(v = beta1.lam, lwd = 3, col = "red")
hist(out19B.5$sims.list$beta2.lam, col = "grey", main = "beta2.lam", xlab = "", breaks = 50)
abline(v = beta2.lam, lwd = 3, col = "red")
hist(out19B.5$sims.list$alpha.p, col = "grey", main = "alpha.p", xlab = "", breaks = 50)
abline(v = alpha.p, lwd = 3, col = "red")
hist(out19B.5$sims.list$beta.p, col = "grey", main = "beta.p", xlab = "", breaks = 50)
abline(v = beta.p, lwd = 3, col = "red")
hist(out19B.5$sims.list$totalN, col = "grey", , main = "Total N", xlab = "", breaks = 50)
abline(v = sum(N), lwd = 3, col = "red")

par(mfrow = c(2,3), mar = c(5,5,5,3), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, ask = TRUE)
  # for(i in 1:nSites){      # Do this for plots of all sites
  for(i in 151:156){         # Do this to produce Fig. 19B.14
 
  hist(out19B.5$sims.list$N[,i], col = "grey", 
    xlim = range(c(max(C[i,]), out19B.5$sims.list$N[,i])), 
    main = paste("Site", i), xlab = "Population size N", freq = FALSE)
  abline(v = Nst[i]-1, lwd = 3, col = "black")
  abline(v = N[i], lwd = 3, lty = 2, col = "red")
}

# Fig. 19B.15
par(mar = c(5,5,5,3), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(mhbElev, N, main = "", xlab = "Elevation (m)", ylab = "Abundance", las = 1, ylim = c(0, max(N)), pch = 16, col = rgb(1,0,0,0.5), cex = 1.3, frame = FALSE)
lines(mhbElev, lambda, col = rgb(1,0,0,0.5), lwd = 5)

Nmix.pred <- exp(out19B.5$mean$alpha.lam + out19B.5$mean$beta1.lam * mhbElevScaled + out19B.5$mean$beta2.lam * mhbElevScaled^2)
lines(mhbElev, Nmix.pred, type = "l", col = rgb(0,0,1,0.5), lwd = 5)
lines(1000*pred.elev, pred, col = rgb(0,0,0,0.5), lwd = 5)
legend("topright", lwd = 5, col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), rgb(0,0,0,0.5)), legend = c("Nmix", "True", "Detection-naive"))

# Compare likelihood with Bayesian estimates and with truth
# Parameter estimates
jags_est <- unlist(out19B.5$mean)[1:5]
comp <- cbind(truth = truth, unmarked = unm_est, JAGS = jags_est)
print(comp, 3)

# Point estimates of optimal elevation (in metres)
comp <- cbind('truth' = 1000*opt.elev.true, 'p-naive' = 1000*opt.elev2, 
  'unmarked' = 1000*opt.elev.um, 'JAGS' = 1000*out19B.5$mean$opt.elev)
print(comp, 3)

# Compare three solutions to the SE of optimum elevation
comp <- c(sqrt(var_elev.opt), SE_opt.elev, out19B.5$sd$opt.elev)
names(comp) <- c("Delta method", "Parametric Bootstrap", "Posterior SD")
print(comp, 4)

# Compare total population size in sampled quadrats to truth
jags_Ntotal <- out19B.5$q50$totalN
comp <- c(truth = true_Ntotal, naive = naive_Ntotal, unmarked = unm_Ntotal, JAGS = jags_Ntotal)
print(comp, 2)

# Make predictions for all Switzerland for a prototypical SDM
# Prelims
samps <- out19B.5$sims.list     # Grab all posterior samples
str(samps)                      # Look at overview
nsamp <- length(samps$alpha.lam)# Sample size

# Create array to hold predictions for all of Switzerland
pred.lam2 <- array(NA, dim = c(nsamp, length(chElev)))

# Fill the array using MCMC draws and Swiss elevation data
for(i in 1:nsamp){
  if(i %% 50 == 0) cat(paste("Predicting for MCMC draw", i, "\n"))
  pred.lam2[i,] <- exp(samps$alpha.lam[i] + samps$beta1.lam[i] * chElevScaled + samps$beta2.lam[i] * chElevScaled^2)
}

# Compute posterior medians, sds and 95%CRIs
# (NOTE: we use the median again instead of the mean, since distribution skewed)
pm <- apply(pred.lam2, 2, median)    # Posterior medians
psd <- apply(pred.lam2, 2, sd)       # Posterior sd
lcl <- apply(pred.lam2, 2, function(x) quantile(x, 0.025)) # Lower CRL
ucl <- apply(pred.lam2, 2, function(x) quantile(x, 0.975)) # Upper CRL

# The Bayesian version of the Holy Grail (see Chapter 8 in AHM1 book)
#   (Fig. 19B.15)
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
par(mfrow = c(2, 2), mar = c(3,4,4,6), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pm))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Swiss SDM for the bullfinch\n(estimated true abundance, lambda)", zlim = c(0, 15))
r2 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = psd))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Uncertainty (posterior sd of lambda)", zlim = c(0, 2.5))
r3 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = lcl))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Uncertainty (Lower credible limit)", zlim = c(0, 15))
r4 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = ucl))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, 
  main = "Uncertainty (Upper credible limit)", zlim = c(0, 15))

# Compute estimated Swiss population size, with posterior sd and CRI
post_Ntot <- apply(pred.lam2, 1, sum)

# Plot full posterior distribution (Fig. 19B.17)
hist(post_Ntot, breaks = 60, col = "gold", main = "Swiss national population size of simulated bullfinches (Posterior distribution)", xlab = "Number of territories", ylab = "Density", freq = FALSE)
abline(v = Ntot_true, col = "red", lwd = 4)
abline(v = median(post_Ntot), col = "blue", lwd = 4, lty = 2)

# Compute median, sd and CRI
(pm_Ntot <- median(post_Ntot))
(psd_Ntot <- sd(post_Ntot))
(CRI_Ntot <- quantile(post_Ntot, c(0.025, 0.975)))


# 19B.6 Bayesian analysis with NIMBLE 
# -----------------------------------

# Bundle and summarize data
elev2 <- mhbElevScaled^2         # squared elevation
str(dataList <- list(C=C, elev = mhbElevScaled, elev2 = elev2, 
  nSites = nSites, nVisits = nVisits) )

# Write NIMBLE model file
model19B.6 <- nimbleCode( {

# Priors
mean.lam ~ dunif(0, 1)
alpha.lam <- log(mean.lam)
beta1.lam ~ dnorm(0, 0.001)
beta2.lam ~ dnorm(0, 0.001)
mean.p ~ dunif(0, 1)
alpha.p <- logit(mean.p)
beta.p ~ dnorm(0, 0.001)

# Likelihood
# Biological model for true abundance
for (i in 1:nSites) {             # Loop over sites
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha.lam + beta1.lam * elev[i] + beta2.lam * elev2[i]
}

# Observation model for replicated counts
for (i in 1:nSites) {             # Loop over sites
  for (t in 1:nVisits) {          # Loop over all n observations
    C[i,t] ~ dbin(p[i,t], N[i]) 
    logit(p[i,t]) <- alpha.p + beta.p * elev[i]
  }
}

# Derived quantities
totalN <- sum(N[1:nSites]) # Estimate total pop. size across all sites
opt.elev <- -beta1.lam / (2 * beta2.lam)  # Optimum elevation

} )

# Inits
Nst <- apply(C, 1, max) + 1
inits <- function(){list(N = Nst, mean.lam=runif(1, 0, 1), beta1.lam=rnorm(1, 0, 1), beta2.lam=rnorm(1, 0, 1), mean.p=runif(1), beta.p=rnorm(1, 0, 1))}

# Parameters monitored: same as before
params <- c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p", "mean.lam", "mean.p", "totalN", "opt.elev", "N")
params <- c("alpha.lam", "beta1.lam", "beta2.lam", "alpha.p", "beta.p", "mean.lam", "mean.p", "totalN", "opt.elev")

# MCMC settings
# Number of samples returned is floor((niter-nburnin)/thin)
# ni <- 150000  ;  nb <- 50000  ; nc <- 4  ; nt <- 100  # Like JAGS
ni <- 300000  ;  nb <- 200000  ; nc <- 4  ; nt <- 200   # safer

# Call NIMBLE (ART 38 min), check convergence and summarize posteriors
system.time( out19B.6 <- 
    nimbleMCMC(code = model19B.6,
    constants = dataList,
    inits = inits, monitors = params,
    nburnin = nb, niter = ni, thin = nt, nchains = nc,
    samplesAsCodaMCMC = TRUE)    )
par(mfrow=c(2,2), ask = TRUE); coda::traceplot(out19B.6)
nsum <- nimble_summary(out19B.6, params)
print(as.matrix(nsum), 4)          # summary, not shown
nimble_est <- nsum[1:5,1]          # Save estimates
nimble_Ntotal <- nsum[8,4]         # Posterior median


# 19B.7 Bayesian analysis with Stan 
# ---------------------------------

# Bundle and summarize data
# Build data list (a bit different from JAGS/Nimble)
# Minimum and maximum possible abundance at each site
Kmin <- maxC                # Maximum count at each site
K = max(C) + 20             # Over max count plus some 'buffer'
dataList <- list(C=C, nSites=nSites, nVisits=nVisits,
  elev=mhbElevScaled, elev2= mhbElevScaled ^2, Kmin=Kmin, K=K)
str(dataList)

# Write Stan model
cat(file="model19B_7.stan", "

data{
  int nSites;
  int nVisits;
  int C[nSites, nVisits];
  vector[nSites] elev;
  vector[nSites] elev2;
  int K;
  int Kmin[nSites];
}

parameters{
  real alpha_lam;
  real beta1_lam;
  real beta2_lam;
  real alpha_p;
  real beta_p;
}

transformed parameters{
  vector[nSites] log_lambda;
  vector[nSites] logit_p;
  for (i in 1:nSites){
    log_lambda[i] = alpha_lam + beta1_lam * elev[i] + beta2_lam * elev2[i];
    logit_p[i] = alpha_p + beta_p * elev[i];
  }
}

model{
  vector[nSites] lik; //likelihood for each site

  //Priors
  alpha_lam ~ normal(0, 2);
  beta1_lam ~ normal(0, 100);
  beta2_lam ~ normal(0, 100);
  alpha_p ~ normal(0, 10);
  beta_p ~ normal(0, 100);

  for (i in 1:nSites){
    lik[i] = 0;
    for (k in Kmin[i]:K){
      lik[i] += exp(poisson_log_lpmf(k | log_lambda[i]) +
                    binomial_logit_lpmf(C[i,1:nVisits] | k, logit_p[i]));
    }
    target += log(lik[i]);
  }
}

generated quantities {
  real opt_elev;
  opt_elev = -beta1_lam / (2 * beta2_lam);
}
")

# Parameters to estimate
params <- c("alpha_lam", "beta1_lam", "beta2_lam", "alpha_p", 
  "beta_p", "opt_elev")

# HMC settings
ni <- 1500   ;   nb <- 500   ;  nc <- 4   ;  nt <- 1

# Call STAN (ART 18 min), assess convergence and print results table
options(mc.cores = parallel::detectCores()-1) # Run Stan in parallel
system.time(
out19B.7 <- stan(file = "model19B_7.stan", data = dataList, 
  pars = params, warmup = nb, iter = ni, chains = nc, thin = nt) )
rstan::traceplot(out19B.7)        # not shown, but see next paragraph
print(out19B.7, dig = 3)          # not shown
stan_est <- summary(out19B.7)$summary[1:5,1]   # Save estimates

# Get posterior of sumN

# Get posterior of parameters (beta vector)
stan_beta_post <- as.matrix(out19B.7)[,1:5]

# Function that takes posterior of beta and generates posterior of sumN
get_post_sumN <- function(Beta_samples, K){
  nsims <- nrow(Beta_samples)
  sumN_sims <- rep(NA, nsims)
  for (i in 1:nsims){
    beta <- Beta_samples[i,]
    lambda <- exp(beta[1] + beta[2] * mhbElevScaled + beta[3] * mhbElevScaled^2)
    p <- plogis(beta[4] + beta[5] * mhbElevScaled)

    N <- rep(NA, nSites)
    for (n in 1:nSites){
      Kvals <- max(C[n,]):K
      Kprob <- rep(NA, length(Kvals))

      for (k in 1:length(Kvals)){
        pLam <- dpois(Kvals[k], lambda[n])
        pP <- 1
        for (j in 1:nVisits){
          pP <- pP * dbinom(C[n,j], Kvals[k], p[n])
        }
        Kprob[k] <- pLam * pP
      }
      Kprob <- Kprob / sum(Kprob)
      N[n] <- sample(Kvals, 1, prob=Kprob)
    }
    sumN_sims[i] <- sum(N)
  }
  sumN_sims
}

# Bootstrap that function to get CIs as well and produce plot
system.time(     # ART 210 sec
   post_sumN <- get_post_sumN(stan_beta_post, K=30)
   )
stan_Ntotal <- median(post_sumN)      # Note use of posterior median again

# Draw figure (Fig. 19B.18)
hist(post_sumN, main="Posterior of sumN from Stan", breaks = 40)
abline(v=stan_Ntotal, col='blue', lwd = 5)
abline(v=true_Ntotal, col='red', lwd = 5)
legend("topright", lty=1, col=c('red','blue'), legend=c('truth','estimate'), lwd = 5, cex = 1.2, bty = 'n')


# 19B.8 Bayesian analysis with canned functions in the R package ubms 
# -------------------------------------------------------------------

library(ubms)

# Load unmarked, format data and summarize
# Produces same as for unmarked in Section 19B.4 of course
library(unmarked)
summary( umf <- unmarkedFramePCount(y = C, 
  siteCovs = data.frame(elev=mhbElevScaled, elev2=mhbElevScaled^2) ) )

# Fit the model in ubms in parallel (ART 7 min)
options(mc.cores=4) # number of parallel cores to use
system.time(
   out19B.8 <- stan_pcount(~elev ~ elev + elev2, umf, chains=4))
out19B.8
ubms_est <- coef(out19B.8)    # Save estimates

# Compare estimates with truth
comp <- cbind(truth = truth, unmarked = unm_est, JAGS = jags_est, 
  NIMBLE = nimble_est, Stan = stan_est, ubms = ubms_est)
print(comp, 3)

# Compare total population size to truth
postN_ubms <- posterior_predict(out19B.8, "z") # Get posterior of latent abundance 
ubms_NtotalX <- apply(postN_ubms, 1, sum) # Calculate sumN for each draw
ubms_Ntotal <- median(ubms_NtotalX)       # Calculate posterior median
comp <- c(truth = true_Ntotal, naïve = naive_Ntotal,
  unmarked = unm_Ntotal, JAGS = jags_Ntotal, NIMBLE = nimble_Ntotal,
  Stan = stan_Ntotal, ubms = ubms_Ntotal)
print(comp, 2)


# 19B.9 Do-it-yourself maximum likelihood estimates 
# -------------------------------------------------

# Definition of NLL for Binomial N-mixture model for lizard counts
# Variant 1 (from Royle & Dorazio, 2008), faster than Variant 2 below
NLL <- function(param, C, x, K, nSites, nVisits){
  alpha.lam <- param[1]    # Abundance intercept (log scale)
  beta1.lam  <- param[2]   # Abundance slope on x
  beta2.lam  <- param[3]   # Abundance slope on x^2
  alpha.p   <- param[4]    # Detection intercept (logit scale)
  beta.p    <- param[5]    # Detection slope on x

  # Covariate model for expected abundance lambda ('GLM 1')
  lambda <- exp(alpha.lam + beta1.lam * x + beta2.lam * x^2)
  # Covariate model for detection probability p ('GLM 2')
  p <- plogis(alpha.p + beta.p*x)

  L <- rep(NA, 3) # Vector for likelihood contribution of each site
  tmp.like <- matrix(NA, nrow = K+1, ncol = nVisits)
  for(i in 1:nSites){
    gN <- dpois(0:K, lambda[i])
    gN <- gN / sum(gN)
    dvec <- C[i,]    # Extract counts at site i
    tmp.like[,1] <- dbinom(dvec[1], 0:K, p[i])
    tmp.like[,2] <- dbinom(dvec[2], 0:K, p[i])
    tmp.like[,3] <- dbinom(dvec[3], 0:K, p[i])
    likvec <- apply(tmp.like, 1, prod)
    L[i] <- sum(likvec * gN)  # Contribution to L from 1 site
  }
  NLL <- -1 * sum(log(L))
  return(NLL)
}

# Variant 2: much slower, but easier to understand than Variant 1

NLL <- function(param, C, x, K, nSites, nVisits){
  alpha.lam <- param[1]    # Abundance intercept (log scale)
  beta1.lam  <- param[2]   # Abundance slope on x
  beta2.lam  <- param[3]   # Abundance slope on x^2
  alpha.p   <- param[4]    # Detection intercept (logit scale)
  beta.p    <- param[5]    # Detection slope on x

  # Covariate model for expected abundance lambda ('GLM 1')
  lambda <- exp(alpha.lam + beta1.lam * x + beta2.lam * x^2)
  # Covariate model for detection probability p ('GLM 2')
  p <- plogis(alpha.p + beta.p*x)

  L <- rep(NA, nSites) # Vector for likelihood contribution of each site

  for(i in 1:nSites){
    pN <- dpois(0:K, lambda[i])
    pY <- rep(1, K+1)
    for (k in 0:K){
      for (j in 1:nVisits){
        pY[k+1] <- pY[k+1] * dbinom(C[i,j], k, p[i])
      }
    }
    L[i] <- sum(pN * pY)
  }
  NLL <- -1 * sum(log(L))
  return(NLL)
}

# Minimize that NLL to find MLEs and get SEs (ART 23 or 55 sec)
inits <- c('alpha.lam' = 0, 'beta1.lam' = 0, 'beta2.lam' = 0,
  'alpha.p' = 0, 'beta.lp' = 0)
inits <- rep(0, 5)
names(inits) <- names(truth)
system.time( out19B.9 <- optim(inits, NLL, C = C, x = mhbElevScaled, 
  K = max(C) + 100, nSites = nSites, nVisits = 3, hessian=TRUE, 
  method = "BFGS", control=list(trace=1, REPORT=2)) )
get_MLE(out19B.9, 4)
diy_est <- out19B.9$par      # Save estimates

# Get AIC: 2 * deviance + 2 * number of parameters
(AIC <- 2 * out19B.9$value + 2 * length(out19B.9$par))

# Compare with the solutions and AIC from unmarked
summary(out19B.4)               # not shown

# Get estimate of sum N (ART 50 sec)
library(MASS)
Beta <- out19B.9$par
Sigma <- solve(out19B.9$hessian)
param_samples <- mvrnorm(1000, Beta, Sigma)
system.time(
  post_sumN <- get_post_sumN(param_samples, 30)
)
diy_Ntotal <- median(post_sumN)

# Plot posterior of totalN from DIY-MLEs -- not shown
hist(post_sumN, main="Posterior of sumN from DIY MLEs", breaks = 40)
abline(v=median(post_sumN), col='blue', lwd = 5)
abline(v=true_Ntotal, col='red', lwd = 5)
legend("topright", lty=1, col=c('red','blue'), legend=c('truth','estimate'), lwd = 4, cex = 1, bty = 'n')


# 19B.10 Likelihood analysis with TMB 
# -----------------------------------

# Bundle and summarize data
Kmin <- apply(C, 1, max)
K = max(C) + 20
dataList <- list(C=C, nSites=nSites, nVisits=nVisits,
  elev=mhbElevScaled, elev2= mhbElevScaled^2, Kmin=Kmin, K=K)
str(dataList)              # not shown

# Write TMB model file
cat(file="model19B_10.cpp",
"#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Describe input data
  DATA_MATRIX(C);
  DATA_INTEGER(nSites);
  DATA_INTEGER(nVisits);
  DATA_VECTOR(elev);
  DATA_VECTOR(elev2);
  DATA_INTEGER(K);
  DATA_IVECTOR(Kmin);

  //Describe parameters
  PARAMETER(alpha_lam);
  PARAMETER(beta1_lam);
  PARAMETER(beta2_lam);
  PARAMETER(alpha_p);
  PARAMETER(beta_p);

  Type loglik = 0.0;    //Initialize log-likelihood at 0

  vector<Type> lambda(nSites);
  vector<Type> p(nSites);
  vector<Type> lik(nSites);

  for (int i=0; i<nSites; i++){
    lambda(i) = exp(alpha_lam + beta1_lam * elev(i) + beta2_lam * elev2(i));
    p(i) = invlogit(alpha_p + beta_p * elev(i));

    lik(i) = 0;
    for (int k=Kmin(i); k<(K+1); k++){
      Type pN = dpois(Type(k), lambda(i), false);
      Type pY = 1;
      for (int j=0; j<nVisits; j++){
        pY *= dbinom(C(i,j), Type(k), p(i), false);
      }
      lik(i) += pN * pY;
    }
    loglik += log(lik(i));
  }

  Type opt_elev = -beta1_lam / (2 * beta2_lam);
  ADREPORT(opt_elev);

  return -loglik;
}
")

# Compile and load TMB function
compile("model19B_10.cpp")
dyn.load(dynlib("model19B_10"))

# Provide dimensions and starting values for parameters
params <- list(alpha_lam = 0, beta1_lam = 0, beta2_lam = 0, 
  alpha_p = 0, beta_p = 0)

# Create TMB object
out19B.10 <- MakeADFun(data = dataList, parameters = params,
                     DLL = "model19B_10", silent=TRUE)

# Optimize TMB object and print results
starts <- rep(0, 5)
opt <- optim(starts, fn = out19B.10$fn, gr = out19B.10$gr, 
  method = "BFGS", hessian = TRUE)
(tsum <- tmb_summary(out19B.10))
tmb_est <- opt$par      # Save estimates

# Get estimate of sum N from TMB solutions and then plot
bootstrap.n <- 10000       # Choose how many samples from MVN to draw
Beta <- tmb_est
Sigma <- solve(opt$hessian)
param_samples <- mvrnorm(bootstrap.n, Beta, Sigma)
system.time(
  post_sumN <- get_post_sumN(param_samples, 30)  # ART 9 min
)
tmb_Ntotal <- median(post_sumN)       # Compute posterior median

# Plot (not shown)
hist(post_sumN, main="Posterior of sum(N) from TMB", breaks=50)
abline(v=true_Ntotal, col='red', lwd = 5)
abline(v=post_sumN, col='blue', lwd = 5)
legend('topright', col=c('red','blue'), lty=1, lwd = 4, legend=c('truth','estimate (median)'), bty = 'n')


# 19B.11 Comparison of the estimates 
# ----------------------------------

# Compare results with truth
comp <- cbind(cbind(truth = truth, unmarked = unm_est, JAGS = jags_est, 
   NIMBLE = nimble_est, Stan = stan_est, ubms = ubms_est,
   DIY = diy_est, TMB = tmb_est))
print(comp, 4)

# Compare total population size to truth
comp <- c(truth = true_Ntotal, naive = naive_Ntotal, 
  unmarked = unm_Ntotal, JAGS = jags_Ntotal, NIMBLE = nimble_Ntotal,
  Stan = stan_Ntotal, ubms = ubms_Ntotal, DIY = diy_Ntotal,
  TMB = tmb_Ntotal)
print(comp, 2)


# 19B.12 Summary and outlook 
# --------------------------
# (no code)
