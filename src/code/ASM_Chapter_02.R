## Code from file Chap 02.v.02.docx
## Datetime: 2023-07-07 13:13:35
## Loading packages in case they are not loaded in text
suppressMessages(library(jagsUI))
suppressMessages(library(nimble))
suppressMessages(library(rstan))
suppressMessages(library(TMB))
suppressMessages(library(ISMbook))
suppressMessages(library(DHARMa))

# Look up help text for the Poisson
?dpois                 # Also '?Poisson'

# Poisson probability mass function (PMF) evaluated for y = 0
dpois(0, lambda = 2)   # Probability of getting a value 0
dpois(0, lambda = 2, log = TRUE)     # Same on log scale

# Get same density 'by hand' to emphasize dpois() is just shorthand !
lam <- 2; y <- 0
(lam^y)*exp(-lam) / factorial(y) # Probability of getting a value 0

# Cumulative distribution function (CDF), or just 'distr. function'
ppois(3, lambda = 2)   # Probability of getting a value of 0, 1, 2 or 3
sum(dpois(0:3, 2))     # Same, summing up probs 'by hand'
plot(ppois(0:9, 2), type = 'h', lend = 'butt', lwd = 20, ylim = c(0, 1), xlab = 'Value of Y', ylab = 'P(Y <= y)', main = "CDF of Y") # not shown

# Quantile function: value of y for which CDF(y) has a certain value
qpois(0.85, lambda = 2)

# Random number generator (RNG) function
set.seed(2016)           # Set seed if want same numbers as we have
rpois(n = 10, lambda = 2)# 10 Poisson(2) random numbers


# Plot PMF of Poisson with lambda = 2 (Fig. 2–2)
par(mar = c(6,8,5,4), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
plot(0:10, dpois(0:10, lambda = 2), xlab = 'Realized value of r.v.', ylab = 'Probability', main = 'Poisson(2) PMF', type = 'h', lend = 'butt', lwd = 20, col = 'gray30', frame = FALSE, las = 1)


?dnorm  # Also '?Normal': Look up the help text for the Normal

# Gaussian, or normal, probability density function (PDF) for y = 650
dnorm(650, mean = 600, sd = 30)
dnorm(650, mean = 600, sd = 30, log = TRUE)     # Same on log scale

# Get the density 'by hand' to emphasize dnorm() is just a shortcut
mu <- 600; sig <- 30; y <- 650
1/(sig*sqrt(2*pi)) * exp(-(y - mu)^2 / (2 * sig^2))

# Cumulative distribution function (CDF), or just 'distr. function'
pnorm(600, mean = 600, sd = 30) # Prob. of value between -Inf and 600
# Next is same, but integrating the PDF 'by hand'
f <- function(x, mean, sd) dnorm(x, mean = mean, sd = sd)
integrate(f, lower = -Inf, upper = 600, mean = 600, sd = 30)
plot(pnorm(seq(400, 800, by = 1), 600, 30), type = 'l', lwd = 5, ylim = c(0, 1), xlab = 'Value of Y', ylab = 'Prob. density', main = "CDF of Y", frame = FALSE) # not shown

# Quantile function: value of y for which CDF(y) has a certain value
qnorm(0.95, mean = 600, sd = 30)

# Random number generator (RNG) function
set.seed(2016)
rnorm(n = 10, mean = 600, sd = 30)# 10 Normal(600, 30) random numbers


# Compute probabilities for classes of binned Normal rv
limits <- seq(500, 700, by = 5)
midpts <- seq(500 + 5/2, 700 - 5/2, by = 5)
cumProb <- pnorm(limits, mean = 600, sd = 30)
probs <- diff(cumProb)

# Plot probability density function (PDF) of Normal(600, 30)
par(mfrow = c(1, 2), mar = c(6,6,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)               # Fig. 2–3
curve(dnorm(x, mean = 600, sd = 30), 500, 700, xlab = 'Body mass (g)', ylab = 'Probability density', main = 'PDF of Normal(600, 30) r.v.', type = 'l', lwd = 3, col = 'gray30', frame = FALSE, las = 1)

# Plot probabilities for 5g-binned mass random variable
plot(midpts, probs, xlab = 'Binned body mass (g)', ylab = 'Probability', main = 'PMF of binned Normal r.v.', type = 'h', lend = 'butt', lwd = 10, col = 'gray30', frame = FALSE, las = 1)


# Look how different a beta can look depending on its parameters !
par(mfrow = c(2,2))       # not shown
curve(dbeta(x,10,10), 0, 1, lwd = 3, col = 1, main = "beta(10, 10)")
curve(dbeta(x,0.1,0.1), 0, 1, lwd = 3, col = 1, main = "beta(0.1, 0.1)")
curve(dbeta(x,1,100), 0, 1, lwd = 3, col = 1, main = "beta(1, 100)")
curve(dbeta(x,100,1), 0, 1, lwd = 3, col = 1, main = "beta(100, 1)")


set.seed(2)

# Two data sets representing counts with expected count lam
lam <- 2
y1 <- rpois(10, lambda = lam)  # Small data set (n = 10)
y2 <- rpois(100, lambda = lam) # Large data set (n = 100)
table(y1)  ;  table(y2)


# Define likelihood function in R
L <- function(lambda, y) {
  Li <- dpois(y, lambda) # likelihood contribution of each data point i
  L <- prod(Li)     # Likelihood for entire data set is a product
  return(L)
}

# Define log-likelihood function in R
LL <- function(lambda, y) {
  LLi <- dpois(y, lambda, log = TRUE) # log-likelihood contribution of i
  LL <- sum(LLi)    # Log-likelihood for entire data set is a sum
  return(LL)
}

# Define negative log-likelihood function in R
NLL <- function(lambda, y) {
  LL <- dpois(y, lambda, log=TRUE) # log-likelihood contribution of i
  NLL <- -sum(LL)   # *neg* log-likelihood for entire data set is a sum
  return(NLL)
}


# Evaluate all three functions for large range of possible values for lambda
possible.lam <- seq(0.001, 4, by = 0.001)
like1 <- loglike1 <- negloglike1 <- numeric(length(possible.lam))
for(i in 1:length(like1)){
  like1[i] <- L(possible.lam[i], y1)         # Likelihood
  loglike1[i] <- LL(possible.lam[i], y1)     # log-Likelihood
  negloglike1[i] <- NLL(possible.lam[i], y1) # negative log-likelihood
}

# Plot profiles (this is Fig. 2-5)
par(mfrow=c(1, 3), mar=c(6,6,5,3), cex.lab=1.5, cex.axis=1.5, cex.main=2)
plot(possible.lam, like1, xlab = 'lambda', ylab = 'L', main = 'Likelihood', frame = FALSE)
abline(v = possible.lam[which(like1 == max(like1))], col = 'grey', lwd = 3, lty = 3)
plot(possible.lam, loglike1, xlab = 'lambda', ylab = 'LL', main = 'log-Likelihood', frame = FALSE)
abline(v = possible.lam[which(loglike1 == max(loglike1))], col = 'grey', lwd = 3, lty = 3)
plot(possible.lam, negloglike1, xlab = 'lambda', ylab = '-LL', main = 'Negative log-Likelihood', frame = FALSE)
abline(v = possible.lam[which(negloglike1 == min(negloglike1))], col = 'grey', lwd = 3, lty = 3)

# Determine MLE 'by hand'
possible.lam[which(like1 == max(like1))]           # Maximum likelihood
possible.lam[which(loglike1 == max(loglike1))]     # Maximum log-likelihood
possible.lam[which(negloglike1 == min(negloglike1))] # Minimum negative log-likelihood


# Determine MLE using function minimisation in R with optim()
inits <- c('lambda' = 1)
out <- optim(inits, NLL, y = y1) # Optimize function for y1 over lambda
out
(MLE <- out$par)                 # Grab the MLEs


# Evaluate L, LL and NLL for large sample
like2 <- loglike2 <- negloglike2 <- numeric(length(possible.lam))
for(i in 1:length(like2)){
  like2[i] <- L(possible.lam[i], y2)         # Likelihood
  loglike2[i] <- LL(possible.lam[i], y2)     # log-Likelihood
  negloglike2[i] <- NLL(possible.lam[i], y2) # negative log-likelihood
}

# Compute value of the maximized log-likelihood for both data sets
(lmax1 <- max(loglike1, na.rm = TRUE))
(lmax2 <- max(loglike2, na.rm = TRUE))

# Plot LL in the vicinity of the MLE for both samples (Fig. 2-6)
ylim <- c(-5, 0)      # Choose scaling of y axis
plot(possible.lam, loglike1-max(loglike1), xlab = 'lambda', 
ylab = 'Scaled Log-likelihood', type = 'l', frame = FALSE, col = 'gray', 
lwd = 5, ylim =ylim, xlim = c(1, 3))
lines(possible.lam, loglike2-max(loglike2), col = 'gray40', lwd = 5)
abline(v = possible.lam[which(loglike1==max(loglike1))], lwd = 3, lty = 3, col = 'grey')      # MLE for small data set
abline(v = possible.lam[which(loglike2==max(loglike2))], lwd = 3, lty = 3, col = 'grey40')    # MLE for large data set
legend('bottomright', col = c('gray', 'gray40'), lwd = 4,
  bty = 'n', legend = c("n = 10", "n = 100"), cex = 1.5)


# Method 1 for CIs: profile-likelihood- or LRT-based 95% CI
# ---------------------------------------------------------
(lmax <- max(loglike1, na.rm = TRUE))   # Maximized ll
(llim <- lmax - 1.92)   # value of that is 3.84/2 units below

# Plot profile focused on some range around MLE (Fig. 2-7)
par(mar = c(6, 6, 6, 3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
ooo <- which(abs(loglike1-lmax) < 3.5)
xlim <- c(0.9 * possible.lam[min(ooo)], 1.1*possible.lam[max(ooo)]) 
ylim <- c(-20, -16)
plot(possible.lam[ooo], loglike1[ooo], type = 'l', lwd = 3, 
  xlab = 'lambda', ylab = 'Log-likelihood', frame = FALSE, xlim = xlim, 
  ylim = ylim, main = 'Log-likelihood curve in vicinity of the MLE\n with LRT-based 95% CI')

# Add to graph vertical line that is chisq(1) / 2 + max(ll)
idx <- which(loglike1 == max(loglike1)) # order which contains max
arrows(possible.lam[idx], llim, possible.lam[idx], lmax, 
  length = 0.25, angle = 30, code = 3, lwd = 2)
abline(h = llim, lty = 2, lwd = 3)
text(2.4, mean(c(llim, lmax)), "1.92", cex = 2, srt = 90)


# Compute two CI limits separate for left and right branch of profile
l.left <- loglike1[1:idx]                      # Left branch
l.right <- loglike1[(idx+1):length(loglike1)]  # Right branch

# Compare lower limit (llim) with values in l.left and then with those in l.right
idx.left <- which(abs(l.left - llim) == min(abs(l.left - llim), na.rm = TRUE))
idx.right <- which(abs(l.right - llim) == min(abs(l.right - llim), na.rm = TRUE))

# Compute profile confidence limits ....
(LCL <- possible.lam[1:idx][idx.left])
(UCL <- possible.lam[(idx+1):length(loglike1)][idx.right])

# add the confidence limits into the plot
segments(possible.lam[idx.left], -21, possible.lam[idx.left], llim, 
  lwd = 2, lty = 3)
segments(possible.lam[idx+idx.right], -21, possible.lam[idx+idx.right],
  llim, lwd = 2, lty = 3)
arrows(possible.lam[idx.left], -20, possible.lam[idx+idx.right],
  -20, length = 0.25, angle = 30, code = 3, lwd = 2)
text(2.4, -19.5, "95% CI", cex = 2, srt = 0)

# Compare with canned profile CI method in R for our Poisson GLM
exp(confint(glm(y1 ~ 1, family = 'poisson')))


# Method 2 for SEs and CIs: 'inverting the Hessian'
# -------------------------------------------------
inits <- c(lambda = 1)

# Small data set (n = 10)
(out1 <- optim(inits, NLL, y = y1, hessian=TRUE))
(MLE1 <- out1$par)                 # Grab MLEs
(VC1 <- solve(out1$hessian))       # Get variance-covariance matrix
(VC1 <- 1/out1$hessian)            # Same for a one-parameter model !
(ASE1 <- sqrt(diag(VC1)))          # Extract asymptotic SEs

# Large data set (n = 100)
(out2 <- optim(inits, NLL, y = y2, hessian=TRUE))
(MLE2 <- out2$par)                 # Grab MLEs
(VC2 <- solve(out2$hessian))       # Get variance-covariance matrix
(VC2 <- 1/out2$hessian)            # Same for a one-parameter model !
(ASE2 <- sqrt(diag(VC2)))          # Extract asymptotic SEs

# Compare with the truth MLEs and ASE's for small and large data set
print(cbind('truth' = lam, MLE1, ASE1, MLE2, ASE2), 5)


# Do-it-yourself Wald-type CI for smaller data set
print(CI1 <- c(MLE1 - 1.96 * ASE1, MLE1 + 1.96 * ASE1), 4)

# LRT-based CI for smaller data set (from R)
print(exp(confint(glm(y1 ~ 1, family = 'poisson'))), 4)

# DIY Wald-type CI for larger data set
print(CI2 <- c(MLE2 - 1.96 * ASE2, MLE2 + 1.96 * ASE2), 4)

# LRT-based CI for larger data set (from R)
print(exp(confint(glm(y2 ~ 1, family = 'poisson'))), 4)


# Method 3.1 for SEs and CIs: Parametric bootstrap
# ------------------------------------------------
simrep <- 10000        # Number of bootstrap samples: a large number
estiPB <- array(NA, dim = c(simrep, 2)) # Array to hold estimates
colnames(estiPB) <- c('Small data set', 'Large data set')

for(i in 1:simrep){
  if(i %% 500 == 0) cat(paste("Iter", i, "\n")) # Counter
  yb1 <- rpois(10, lambda = MLE1)    # Draw another small data set
  fm1 <- summary(glm(yb1~1, family = 'poisson')) # re-fit the same model
  estiPB[i,1] <- exp(fm1$coef[1])    # Save estimate on natural scale
  yb2 <- rpois(100, lambda = MLE2)   # Draw another large data set
  fm2 <- summary(glm(yb2~1, family = 'poisson')) # re-fit the same model
  estiPB[i,2] <- exp(fm2$coef[1])    # Save estimates
}

# Parametrically bootstrapped sampling distr. of the MLE (Fig. 2–8)
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)
xlim <- range(estiPB[,1])
hist(estiPB[,1], col = 'grey', main = 'Bootstrapped sampling distribution of MLE\nof lambda (small data set; parametric bootstrap)', xlim = xlim)
hist(estiPB[,2], col = 'grey', main = 'Bootstrapped sampling distribution of MLE\nof lambda (large data set; parametric bootstrap)', xlim = xlim)

# Parametric-bootstrapped standard errors
se1.pb <- sd(estiPB[,1])            # Small data set
se2.pb <- sd(estiPB[,2])            # Large data set

# Parametric-bootstrapped CIs
CI1.pb <- quantile(estiPB[,1], c(0.025, 0.975)) # Small data set
CI2.pb <- quantile(estiPB[,2], c(0.025, 0.975)) # Large data set



# Method 3.2 for SEs and CIs: Nonparametric bootstrap
# ---------------------------------------------------
simrep <- 10000       # Number of bootstrap samples: a large number
estiNPB <- array(NA, dim = c(simrep, 2)) # Array to hold estimates
colnames(estiNPB) <- c('Small data set', 'Large data set')

for(i in 1:simrep){
  if(i %% 500 == 0) cat(paste("Iter", i, "\n")) # Counter
  yb1 <- sample(y1, 10, replace = TRUE)  # Re-sample a small data set
  fm1 <- summary(glm(yb1~1, family = 'poisson')) # Re-fit model
  estiNPB[i,1] <- exp(fm1$coef[1])       # Save estimate
  yb2 <- sample(y2, 100, replace = TRUE) # Re-sample a large data set
  fm2 <- summary(glm(yb2~1, family = 'poisson')) # Re-fit model
  estiNPB[i,2] <- exp(fm2$coef[1])       # Save estimate
}

# Non-parametrically bootstrapped sampling distributions of the MLE
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.6)          # not shown
xlim <- range(estiNPB[,1])
hist(estiNPB[,1], col = 'grey', main = 'Bootstrapped sampling distribution of MLE\nof lambda (small data set; non-parametric bootstrap)', xlim = xlim)
hist(estiNPB[,2], col = 'grey', main = 'Bootstrapped sampling distribution of MLE\nof lambda (large data set; non-parametric bootstrap)', xlim = xlim)

# Nonparametric-bootstrapped standard errors
se1.npb <- sd(estiNPB[,1])               # Small data set
se2.npb <- sd(estiNPB[,2])               # Large data set

# Nonparametric-bootstrapped CIs
CI1.npb <- quantile(estiNPB[,1], c(0.025, 0.975))# Small data set
CI2.npb <- quantile(estiNPB[,2], c(0.025, 0.975))# Large data set



# Standard errors
# ---------------
# Small data set: asymptotic normality and the two bootstraps
print(cbind(ASE1, 'PB-SE' = se1.pb, 'NPB-SE' = se1.npb), 4)

# Large data set: asymptotic normality and the two bootstraps
print(cbind(ASE2, 'PB-SE' = se2.pb, 'NPB-SE' = se2.npb), 4)

# 95% Confidence intervals
# ------------------------
# Quickly re-calculate LRT-based CI
CIp1 <- exp(confint(glm(y1 ~ 1, family = 'poisson')))
CIp2 <- exp(confint(glm(y2 ~ 1, family = 'poisson')))

# CIs for small data set: LRT-based, normal approximation, bootstraps
print(cbind(CIp1, CI1, CI1.pb, CI1.npb), 4)

# CIs for large data set: LRT-based, normal approximation, bootstraps
print(cbind(CIp2, CI2, CI2.pb, CI2.npb), 4)


# Swiss Bee-eater data set (national counts of pairs 1990-2020)
# -------------------------------------------------------------
# Counts of known pairs in the country 1990-2020
y <- c(0, 2, 7, 5, 1, 2, 4, 8, 10, 11, 16, 11, 10, 13, 19, 31, 20, 26, 19, 21, 34, 35, 43, 53, 66, 61, 72, 120, 102, 159, 199)
year <- 1990:2020     # Define year range
x <- (year-1989)      # Scaled, but not centered, year covariate
x <- x-16             # Now year covariate is also centered

# Plot bee-eater data (Fig. 2–10)
par(mfrow = c(1, 2), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
hist(y, xlab = 'Count (y)', ylab = 'Frequency', breaks = 20,
  col = 'gray', main = 'Frequency distribution of counts')
plot(year, y, xlab = 'Year (x)', ylab = 'Count (y)', frame = FALSE, cex = 1.5,
  pch = 16,  col = 'gray20', main = 'Relationship y ~ x')
fm <- glm(y ~ x, family = 'poisson') # Add Poisson GLM line of best fit
lines(year, predict(fm, type = 'response'), lwd = 3, col = 'red', lty = 3)


# Define an R function for the NLL of a Poisson log-linear model
NLL <- function(param, y, x) {
  alpha <- param[1]        # Intercept
  beta <- param[2]         # Slope
  lambda <- exp(alpha + beta * x)
  LLi <- dpois(y, lambda, log=TRUE) # log-likelihood for an observation
# LLi <- y*log(lambda)-lambda-lfactorial(y) # Same 'by hand' !
  LL <- -sum(LLi)          # NLL for all observations in vector y
  return(LL)
}

# Try out, evaluate a few candidates for alpha, beta (lower = better)
NLL(c(0, 0), y, x)
NLL(c(0, 1), y, x)
NLL(c(1, 0), y, x)


# Brute-force search for the MLEs in grid
nside <- 1000   # Grid resolution governs quality of approximation
try.alpha <- seq(2.8, 2.96, length.out = nside)  # intercept
try.beta <- seq(0.135, 0.16, length.out = nside) # slope
nll.grid <- array(NA, dim = c(nside, nside))

# Evaluate NLL over the grid: try 1 Million of (alpha, beta) pairs
for(i in 1:nside){
  for(j in 1:nside){  
     nll.grid[i,j] <- NLL(c(try.alpha[i], try.beta[j]), y, x)
  }
}

# Get pair of values associated with minimum function value
(best <- which(nll.grid == min(nll.grid)))
(best.point <- cbind(rep(try.alpha, nside), rep(try.beta, each = nside))[best,])


# Plot the likelihood surface, or likelihood landscape (Fig. 2–11)
par(mar = c(6,6,5,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=try.alpha, y=try.beta, z=nll.grid, col = mapPalette(100),
axes = FALSE, xlab = "Intercept (alpha)", ylab = "Slope (beta)")
contour(x=try.alpha, y=try.beta, z=nll.grid, add = TRUE, lwd = 1.5, labcex = 1.5)
axis(1, at = seq(min(try.alpha), max(try.alpha), by = 0.05))
axis(2, at = seq(min(try.beta), max(try.beta), by = 0.005))
box()
points(best.point[1], best.point[2], col = 'black', pch = 'X', cex = 1.5)  # Best point from our grid search


# Minimize NLL and get VC matrix by inverting the Hessian
inits <- c(alpha = 0, beta = 0)      # Inits
(out <- optim(inits, NLL, y = y, x = x, hessian = TRUE, method = 'BFGS'))
MLE <- out$par                 # Grab MLEs
(VC <- solve(out$hessian))     # Get variance-covariance matrix
ASE <- sqrt(diag(VC))          # Extract asymptotic SEs
LCL <- MLE - 1.96 * ASE        # Lower limit of Wald-type CI
UCL <- MLE + 1.96 * ASE        # Upper limit of Wald-type CI
print(cbind(MLE, ASE, LCL, UCL), 4) # Print MLEs, SEs and CI

# Compare with function glm() and get 95% profile CIs
fm <- glm(y~x, family = 'poisson')
summary(fm)
confint(fm)



# Illustration of the concept of Monte Carlo integration (Fig. 2–13)
set.seed(2)
par(mfrow = c(2,2), mar = c(3,5,4,2))

# 10 draws from distribution
out1 <- rnorm(10, mean = 5, sd = 2)
hist(out1, breaks = 50, col="grey", main="Can't believe it's normal !")
mean(out1)  ;  sd(out1)

# 100 draws from distribution
out2 <- rnorm(100, mean = 5, sd = 2)
hist(out2, breaks = 50, col="grey", main = 'Well, perhaps ...')
mean(out2)  ;  sd(out2)

# 1,000 draws from distribution
out3 <- rnorm(1000, mean = 5, sd = 2)
hist(out3, breaks = 50, col="grey", main = 'Oh !')
mean(out3)  ;  sd(out3)

# 1,000,000 draws from distribution
out4 <- rnorm(10^6, mean = 5, sd = 2)
hist(out4, breaks = 100, col="grey", main = 'Oh, wow!')
mean(out4)  ;  sd(out4)


# Define R function for log-likelihood of Poisson log-linear regression
LL <- function(param, y, x) {
  alpha <- param[1]        # Intercept
  beta <- param[2]         # Slope
  lambda <- exp(alpha + beta * x)   # Log-linear model
  LLi <- dpois(y, lambda, log=TRUE) # LL contribution of each datum i
  LL <- sum(LLi)      # LL for all observations in the data vector y
  return(LL)
}


# Define R function for un-normalized log-posterior of the model
log.posterior <- function(alpha, beta, y, x){
  loglike <- LL(c(alpha, beta), y, x)
  logprior.alpha <- dnorm(alpha, mean = 0, sd = 100, log = TRUE)
  logprior.beta <- dnorm(beta, mean = 0, sd = 100, log = TRUE)
  return(loglike + logprior.alpha + logprior.beta)
}


# Choose number of iterations for which to run the algorithm
niter <- 10000

# Create an R object to hold the posterior draws produced
out.mcmc <- matrix(NA, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# Initialize chains for both parameters: these are the initial values
alpha <- rnorm(1)             # Init for alpha
beta <- rnorm(1)              # Init for beta

# R object to hold acceptance indicator for both params
acc <- matrix(0, niter, 2, dimnames = list(NULL, c("alpha", "beta")))

# Evaluate current value of the log(posterior density) (for the inits)
logpost.curr <- log.posterior(alpha, beta, y, x)

# Choose values for the tuning parameters of the algorithm
sd.tune <- c(0.1, 0.01)         # first is for alpha, second for beta


# Run MCMC algorithm
for(t in 1:niter){
  if(t %% 1000 == 0)     # counter
    cat("iter", t, "\n")
  # First, update log-linear intercept (alpha)
  # ------------------------------------------
  # Propose candidate value of alpha
  alpha.cand <- rnorm(1, alpha, sd.tune[1]) # note tuning 'parameter'

  # Evaluate log(posterior) for proposed new alpha and 
  #   for current beta for the data y and covs x
  logpost.cand <- log.posterior(alpha.cand, beta, y, x)
  # Compute Metropolis acceptance ratio r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate alpha if it meets criterion (u < r)
  if(runif(1) < r){
    alpha <- alpha.cand
    logpost.curr <- logpost.cand
    acc[t,1] <- 1    # Indicator for whether candidate alpha accepted
  }


  # Second, update log-linear slope (beta)
  # -------------------------------------
  beta.cand <- rnorm(1, beta, sd.tune[2]) # note again tuning parameter

  # Evaluate the log(posterior) for proposed new beta and 
  #   for the current alpha for the data y with covs x 
  logpost.cand <- log.posterior(alpha, beta.cand, y, x)
  # Compute Metropolis acceptance ratio
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta <- beta.cand
    logpost.curr <- logpost.cand
    acc[t,2] <- 1    # Indicator for whether candidate beta accepted
  }
  out.mcmc[t,] <- c(alpha, beta) # Save samples for iteration i
  # NOTE: if proposed new values not accepted, then in this step 
  # we will copy their 'old' values and insert them into object 'out'
}


# Compute overall acceptance ratio separately for alpha and beta
(acc.ratio <- apply(acc, 2, mean))  # Should be in range 25-40%

# Check the traceplots for convergence of chains
# Plot all MC draws right from the start (Fig. 2–14)
par(mfrow = c(2, 1))
plot(1:niter, out.mcmc[,1], main = paste('alpha (acc. ratio =', round(acc.ratio[1], 2),')'), xlab = "MC iteration", 
  ylab = 'Posterior draw', type = 'l')
plot(1:niter, out.mcmc[,2], main = paste('beta (acc. ratio =', round(acc.ratio[2], 2),')'), xlab = "MC iteration", 
  ylab = 'Posterior draw', type = 'l')


# Based on that plot choose the amount of burnin
nb <- 500              # Discard this number of draws at the start

# Compute post-burnin acceptance ratio
(acc.ratio.pb <- apply(acc[nb:niter,], 2, mean))

# Repeat traceplots only for post-burnin draws (Fig. 2–15)
par(mfrow = c(2, 1))
plot(nb:niter, out.mcmc[nb:niter,1], main =
  paste('alpha (post-burnin acc. ratio=', round(acc.ratio.pb[1], 2),')'), xlab="Iteration", ylab='Posterior draw', type='l')
plot(nb:niter, out.mcmc[nb:niter,2], main =
  paste('beta (post-burnin acc. ratio=', round(acc.ratio.pb[2], 2),')'), xlab="Iteration", ylab='Posterior draw', type='l')


# Get posterior mean, sd and 95% CRI (all post-burnin)
post.mn <- apply(out.mcmc[nb:niter,], 2, mean)  # posterior mean
post.sd <- apply(out.mcmc[nb:niter,], 2, sd)    # posterior sd
CRI <- apply(out.mcmc[nb:niter,], 2, 
  function(x) quantile(x, c(0.025, 0.975)))     # 95% CRI

# Produce marginal posterior summary for post-burnin part of chains
tmp <- cbind(post.mn, post.sd, t(CRI))
print(tmp, 4)

# Compare Bayesian with DIY-MLE inferences
get_MLE(out, 4)           # function from the ASMbook package


# Plot the joint posterior distribution of alpha and beta (Fig. 2–16)
par(mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
plot(out.mcmc[nb:niter,1], out.mcmc[nb:niter,2], xlab='alpha', 
  ylab='beta', frame=FALSE, cex=1.2, pch=16, col=rgb(0,0,0,0.1))


# Plot the marginal posterior distributions of alpha, beta (Fig. 2–17)
par(mfrow = c(1, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
hist(out.mcmc[nb:niter,1], xlab = 'alpha', col = 'grey', main = '')
hist(out.mcmc[nb:niter,2], xlab = 'beta', col = 'grey', main = '')


# Shorter run than default
# (can use for true.vals the MLEs from above for graphics)
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(0, 0),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.1, 0.1),
  niter=2000, nburn=1000)) # Bad mixing, strong autocorrelation

# Smaller tuning param to increase acceptance
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits = c(0,0),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.1, 0.01),
  niter=2000, nburn=1000))

# Larger step length (=tuning.params)
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(0, 0),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.5, 0.5),
  niter=2000, nburn=1000)) # Ugly Manhattan traceplots, not efficient, no convergence after 1000 steps

# Try to break by picking far-out inits (pretty amazing !)
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(40, 40),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.1, 0.1),
  niter=5000, nburn=2000) )

# Try to break some more ... (perhaps even more amazing)
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(-40, 40),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.1, 0.1),
  niter=5000, nburn=2000) )  # Note need longer burnin

# ...and more 
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(40, 0),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.1, 0.1),
  niter=5000, nburn=2000) )   # Remarkable search trajectory !

# Same with bigger step length (and shorter chain)
str(tmp <- demoMCMC(Cdat=y, Xdat=x, true.vals=MLE, inits=c(40, 0),
  prior.sd.alpha=100, prior.sd.beta=100, tuning.params=c(0.5, 0.2),
  niter=3000, nburn=2000) )   # Remarkable search trajectory !


