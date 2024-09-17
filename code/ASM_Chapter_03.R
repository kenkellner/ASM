
# ---------------------------------------------------------------
# Code for the book 'Applied Statistical Modeling for Ecologists' 
# by Kéry & Kellner, Elsevier, 2024
# ---------------------------------------------------------------

# -----------------------------------------------------------------------
# Chapter 3  --  Linear regression models and their extensions to 
#                generalized linear, hierarchical and integrated models
# -----------------------------------------------------------------------

# Last changes: 11 June 2024


# 3.1 Introduction
# ----------------
# (no code)


# 3.2 Statistical distributions for the random variables in our model
# -------------------------------------------------------------------
# (no code)


#   3.2.1 Normal distribution
# ---------------------------

n <- 100000                                 # Sample size
mu <- 600                                   # Value of population mean
sigma <- 30                                 # Value of population SD
sample <- rnorm(n = n, mean = mu, sd = sigma)     # Draw random numbers
print(sample, dig = 4)                      # ... print them
hist(sample, col = "grey", breaks = 60, freq = F) # ... plot them


#   3.2.2 Uniform distribution
# ----------------------------

n <- 100000                                 # Sample size
a <- lower.limit <- 0                       # Lower limit of RV
b <- upper.limit <- 10                      # Upper limit of RV
sample <- runif(n = n, min = a, max = b)    # Draw random numbers
print(sample, dig = 4)                      # ... print them
hist(sample, col = "grey", breaks = 30, freq = F) # ... plot them


#   3.2.3 Binomial distribution: the "coin-flip distribution"
# -----------------------------------------------------------

n <- 100000                                 # Sample size
N <- 16                                     # Number of individuals that flip a coin
p <- 0.8                                    # Probability of being counted (seen), dead or a male
sample <- rbinom(n = n, size = N, prob = p) # Draw random numbers
print(sample, dig = 4)                      # ... print them
plot(table(sample)/n, ylab = 'Probability', type = 'h', lwd = 20, col = rgb(0, 0, 0, 0.5),
xlim = c(0, 16), lend = 'butt', frame = F)  # ... plot them


#   3.2.4 Poisson distribution
# ----------------------------

n <- 100000                                 # Sample size
lambda <- 5                                 # Average count per window (density or intensity)
sample <- rpois(n = n, lambda = lambda)     # Draw random numbers
print(sample, dig = 4)                      # ... print them
plot(table(sample)/n, ylab = 'Probability', type = 'h', lwd = 20,col = rgb(0, 0, 0, 0.5),
lend = 'butt', frame = F)                   # ... plot them


#   3.2.5 Summary remarks on statistical distributions
# ----------------------------------------------------
# (no code)


# 3.3 Link functions to model parameters on a transformed scale
# -------------------------------------------------------------
# (no code)


# 3.4 Linear modeling of covariate effects
# ----------------------------------------

mass <- c(60, 80, 50, 70, 90, 110)   # Response variable (continuous)
pop <- factor(c(1, 1, 2, 2, 3, 3))   # A categorical covariate (factor)
reg <- factor(c(1, 1, 1, 1, 2, 2))   # Another factor
hab <- factor(c(1, 2, 3, 1, 2, 3))   # Yet another factor
length <- c(40, 45, 39, 50, 52, 56)  # A continuous covariate

# A continuous covariate with quantitative values
mass
# A categorical covariate or factor
pop

#   3.4.1 The "model of the mean"
# -------------------------------

lm(mass ~ 1)

model.matrix(mass ~ 1)


#   3.4.2 Comparison of two groups
# --------------------------------

model.matrix( ~ reg)

lm(mass ~ reg)

lm(mass ~ reg - 1) # not shown
model.matrix( ~ reg - 1)


#   3.4.3 Simple linear regression: modeling the effect of one continuous covariate
# ---------------------------------------------------------------------------------

model.matrix(~ length)

lm(mass ~ length)

lm(mass ~ I(length - mean(length)))

model.matrix( ~ length - 1)


#   3.4.4 Modeling the effects of one categorical covariate, or factor
# --------------------------------------------------------------------

model.matrix(~pop)

lm(mass ~ pop) # Fit model in effects parameterization (R default)

model.matrix(~ pop - 1)

lm(mass ~ pop - 1) # Fit model in means parameterization


#   3.4.5 Modeling the effects of two factors
# -------------------------------------------

lm(mass ~ reg + hab)

model.matrix( ~ reg + hab)

lm(mass ~ reg * hab)

model.matrix(~ reg * hab)

table(reg, hab)

alias(lm(mass ~ reg*hab))

lm(mass ~ reg * hab - 1 - reg - hab)

model.matrix(mass ~ reg * hab - 1 - reg - hab)


#   3.4.6 Modeling the effects of a factor and continuous covariate
# -----------------------------------------------------------------

model.matrix( ~ pop + length)

lm(mass ~ pop + length)              # Default model fitting
DM <- model.matrix(~pop + length)    # create design matrix DM
lm(mass ~ DM-1)                      # Fit same model via explicit DM

par(mfrow = c(1, 2)) # Fig. 3–5 left
fm <- lm(mass ~ pop + length)        # Refit model
plot(length, mass, col = c(rep("red", 2), rep("blue", 2), rep("green", 2)),
pch = 16, cex = 1.5, ylim = c(40, 120), frame = FALSE)
abline(fm$coef[1], fm$coef[4], col = rgb(1, 0, 0, 0.5), lwd = 2)
abline(fm$coef[1] + fm$coef[2], fm$coef[4], col = rgb(0, 0, 1, 0.5), lwd = 2)
abline(fm$coef[1] + fm$coef[3], fm$coef[4], col = rgb(0, 1, 0, 0.5), lwd = 2)

model.matrix( ~ pop * length)

lm(mass ~ pop * length)

fm <- lm(mass ~ pop * length)        # Refit model, Fig. 3–5 right
plot(length, mass, col = c(rep("red", 2), rep("blue", 2), rep("green", 2)),
pch = 16, cex = 1.5, ylim = c(40, 120), frame = FALSE)
abline(fm$coef[1], fm$coef[4], col = rgb(1, 0, 0, 0.5), lwd = 2)
abline(fm$coef[1] + fm$coef[2], fm$coef[4] + fm$coef[5], col = rgb(0, 0, 1, 0.5), lwd = 2)
abline(fm$coef[1] + fm$coef[3], fm$coef[4] + fm$coef[6], col = rgb(0, 1, 0, 0.5), lwd = 2)

model.matrix( ~ pop + length - 1)

lm(mass ~ pop + length - 1)

model.matrix( ~ pop * length - 1 - length)

lm(mass ~ pop * length - 1 - length)
summary(lm(mass ~ pop * length - 1 - length)) # not shown


# 3.5 Brief overview of linear, generalized linear, (generalized) linear
#     mixed, hierarchical, and integrated models
# ----------------------------------------------------------------------
# (no code)


# 3.6 Summary and outlook
# -----------------------
# (no code)
