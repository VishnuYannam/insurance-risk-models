library(MASS)
library(actuar)
library(fitdistrplus)
library(stats)
library(sads)
library(nortest)
library(tibble)
library(kSamples)
library("pbapply")

###TASK 1###
#Read the loss data

##Edit CSV Column name to loss##
loss_data <- read.csv("Loss.csv")
losses <- loss_data$Loss

#Fit the Distributions using MLE and estimate model parameters

#Log-normal Distribution
fit_lognormal <- fitdistr(losses, "lognormal")
mle_mu <- fit_lognormal$estimate["meanlog"]
mle_sigma <- fit_lognormal$estimate["sdlog"]


#Gamma Distribution
fit_gamma <- fitdistr(losses, "gamma")
mle_shape <- fit_gamma$estimate["shape"]
mle_rate <- fit_gamma$estimate["rate"]

#Pareto Distribution
pareto_log_likelihood <- function(params) {
  #Define alpha and xmin from parameters
  alpha <- params[1]
  xmin <- params[2]
  #Ensure parameters are positive
  if (alpha <= 0 || xmin <= 0) return(Inf)
  n <- length(losses)
  - (n * log(alpha) + n * alpha * log(xmin) - (alpha + 1) * sum(log(pmax(losses - xmin, 1e-10))))
}

initial_params <- c(1, min(losses))
mle_result <- optim(par = initial_params, fn = pareto_log_likelihood, method = "L-BFGS-B")

cat("Estimated Shape Parameter (alpha):", mle_result$par[1], "\n")
cat("Estimated Scale Parameter (xmin):", mle_result$par[2], "\n")

#Generate summary table for parameters
mle_mu_lognormal <- fit_lognormal$estimate["meanlog"]
mle_sigma_lognormal <- fit_lognormal$estimate["sdlog"]

mle_shape_gamma <- fit_gamma$estimate["shape"]
mle_rate_gamma <- fit_gamma$estimate["rate"]

alpha_pareto <- mle_result$par[1]
xmin_pareto <- mle_result$par[2]

# Create summary table
summary_table_parameters <- tibble(
  Distribution = c("Log-normal", "Gamma", "Pareto"),
  `Parameter 1` = c(paste("Meanlog =", round(mle_mu_lognormal, 3)),
                    paste("Shape =", round(mle_shape_gamma, 3)),
                    paste("Alpha =", round(alpha_pareto, 3))),
  `Parameter 2` = c(paste("Sdlog =", round(mle_sigma_lognormal, 3)),
                    paste("Rate =", round(mle_rate_gamma, 3)),
                    paste("Xmin =", round(xmin_pareto, 3)))
)

print(summary_table_parameters)

##QQ Plots##
## Generating random numbers based on estimated parameters
random.lnorm <- rlnorm(100000, meanlog = mle_mu, sdlog = mle_sigma)
random.gamma <- rgamma(100000, shape=mle_shape_gamma, rate=mle_rate_gamma)
random.pareto <- rpareto2(100000, min = 0, shape = alpha_pareto, scale = xmin_pareto)

qqplot(losses, random.lnorm)
abline(a=0,b=1, col="red")
qqplot(losses, random.gamma)
abline(a=0,b=1,col="red")
qqplot(losses, random.pareto)
abline(a = 0, b = 1, col = "red")
######################

##KS Tests##
#Log-normal
#KS Test
ks_lognorm <- ks.test(losses, "plnorm", meanlog = mle_mu, sdlog = mle_sigma)
print(ks_lognorm)


#Gamma
#KS Test
ks_gamma <- ks.test(losses, "pgamma", shape = mle_shape, rate = mle_rate)
print(ks_gamma)

#Pareto
#KS Test
pareto_cdf <- function(x, alpha, xmin) {
  ifelse(x < xmin, 0, 1 - (xmin/x)^alpha)
}

ks_pareto <- ks.test(losses, pareto_cdf, alpha = mle_result$par[1], xmin = mle_result$par[2])
print(ks_pareto)

##AIC and BIC for Distributions##
n <- length(losses)

# Log-normal Distribution
k_lognorm <- 2
log_likelihood_lognorm <- sum(dlnorm(losses, meanlog = mle_mu, sdlog = mle_sigma, log = TRUE))
aic_lognorm <- 2 * k_lognorm - 2 * log_likelihood_lognorm
bic_lognorm <- k_lognorm * log(n) - 2 * log_likelihood_lognorm
cat("AIC for Log-normal Distribution:", aic_lognorm, "\n")
cat("BIC for Log-normal Distribution:", bic_lognorm, "\n")

# Gamma Distribution
k_gamma <- 2
log_likelihood_gamma <- sum(dgamma(losses, shape = mle_shape, rate = mle_rate, log = TRUE))
aic_gamma <- 2 * k_gamma - 2 * log_likelihood_gamma
bic_gamma <- k_gamma * log(n) - 2 * log_likelihood_gamma
cat("AIC for Gamma Distribution:", aic_gamma, "\n")
cat("BIC for Gamma Distribution:", bic_gamma, "\n")

# Pareto Distribution
pareto_log_likelihood_val <- -pareto_log_likelihood(c(mle_result$par[1], mle_result$par[2]))
k_pareto <- 2
aic_pareto <- 2 * k_pareto - 2 * pareto_log_likelihood_val
bic_pareto <- k_pareto * log(n) - 2 * pareto_log_likelihood_val
cat("AIC for Pareto Distribution:", aic_pareto, "\n")
cat("BIC for Pareto Distribution:", bic_pareto, "\n")

##SUMMARY TABLE##
# Create a summary table
summary_table <- data.frame(
  Distribution = c("Log-normal", "Gamma", "Pareto"),
  AIC = c(aic_lognorm, aic_gamma, aic_pareto),
  BIC = c(bic_lognorm, bic_gamma, bic_pareto)
)

# Print the summary table
print(summary_table)


###TASK 2###
#### Task 2 #####

set.seed(1)

# Parameters
c0 <- 90
t <- 60
lambda <- 5
theta <- 0.075
a <- 5
b <- 1.5
size <- 1e6
alpha <- 0.7

# Calculate premium and reinsurance premium
premium <- (1 + theta) * lambda * (a / b)
reinsurancepremium <- (1 + loading) * (1 - alpha) * (a / b) * lambda

# Function to check for ruin with proportional reinsurance
ruincheck <- function() {
  claimsnumber <- rpois(1, lambda * t)
  arrtime <- sort(runif(claimsnumber, 0, t))
  claimsize <- alpha * rgamma(claimsnumber, shape = a, rate = b)
  surplus <- c0 + premium * arrtime - cumsum(claimsize)
  min(surplus) < 0
}

# Simulated ruin probability for proportional reinsurance
prob_ruin_proportional <- mean(pbreplicate(size, ruincheck()))

# Excess of loss parameters
d <- 2.57359
loadingB <- 0.125

# Define a function to calculate reinsurance premium for EoL
calculateReinsurancePremium <- function(d, loadingB) {
  integrand <- function(y) {
    return((y - d) * (dgamma(y, shape = a, rate = b)))
  }
  premium <- lambda * (1 + loadingB) * integrate(integrand, lower = d, upper = Inf)$value
  return(premium)
}

premiumReinsuranceB <- calculateReinsurancePremium(d, loadingB)

premiumReinsuranceB



#Function to check for ruin with Excess of Loss (EoL) reinsurance
ruincheckEoL <- function() {
  N <- rpois(1, lambda * t)
  arrtime <- sort(runif(N, 0, t))
  adjclmsize <- pmin(rgamma(N, shape = a, rate = b), d)
  sp <- c0 + (premium - premiumReinsuranceB) * arrtime - cumsum(adjclmsize)
  min(sp) < 0
}

# Simulated ruin probability for EoL reinsurance
prob_ruin_eol <- mean(replicate(size, ruincheckEoL()))

###Q1 (b)###
# Define the equation for root finding
eqn <- function(R) {
  return(-R * (1 + theta) * a / b + (1 - R / b) ^ (-a) - 1)
}
eqn <- Vectorize(eqn)

rt = uniroot(eqn, lower = 0.01, upper = 0.10)
upper_bound_Proportional <- exp(-rt$root * c0)

objective_fn_EoL <- function(R) {
  return((-R * (1 + theta) * a / b + (1 - R / b) ^ (-a) - 1)^2)
}
root_optim_EoL <- optim(par = 0.05, fn = objective_fn_EoL, method = "Brent", lower = 0, upper = 1)
root_value_EoL <- root_optim_EoL$par
upper_bound_EoL <- exp(-root_value_EoL * c0)

results <- data.frame(
  Option = c("Proportional", "EoL"),
  Simulated_Ruin_Probability = c(prob_ruin_proportional, prob_ruin_eol),
  Upper_Bound = c(upper_bound_Proportional, upper_bound_EoL)
)


# Retrieving the optimal values #
optimal_alpha_val <- nlm(function(x) -r_val_optimal_alpha(x), 0.5)$estimate
optimal_d_val <- nlm(r_val_optimal_d, 2)$estimate

optimal_alpha_val
optimal_d_val

# Print the results table
print(results)

