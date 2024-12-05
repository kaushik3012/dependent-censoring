# # Load necessary libraries
# library(copula)
# library(MASS)
# 
# # Set parameters for lognormal margins
# mu_T <- c(2.2, 2.5)  # Mean parameters for T
# sigma_T <- c(1.0, 1.0)  # Standard deviation parameters for T
# mu_C <- c(2.0, 2.0)  # Mean parameters for C
# sigma_C <- c(0.25, 0.50)  # Standard deviation parameters for C
# 
# # Define copula families: Frank, Clayton, Gumbel, Gauss
# copulas <- list(frankCopula(), claytonCopula(), gumbelCopula(), normalCopula())
# 
# # Set Kendall's tau values
# taus <- c(0.2, 0.5, 0.7)
# 
# # Function to generate lognormal margins
# generate_lognormal <- function(n, mu, sigma) {
#   rlnorm(n, meanlog = mu, sdlog = sigma)
# }
# 
# # Function to simulate data for a given copula and sample size
# simulate_data <- function(n, copula, tau, mu_T, sigma_T, mu_C, sigma_C) {
#   # Set copula parameter using Kendall's tau
#   copula <- setTheta(copula, iTau(copula, tau))
#   
#   # Simulate from the copula
#   u <- rCopula(n, copula)
#   
#   # Apply lognormal margins
#   T <- generate_lognormal(n, mu_T, sigma_T)
#   C <- generate_lognormal(n, mu_C, sigma_C)
#   
#   # Get observed time Y = min(T, C) and censoring indicator delta
#   Y <- pmin(T, C)
#   delta <- as.numeric(T <= C)
#   
#   return(data.frame(Y, delta))
# }
# 
# # Simulation settings
# n_values <- c(200, 500)  # Sample sizes
# n_replications <- 200  # Number of replications
# results <- list()
# 
# # Perform simulation for each scenario
# for (scenario in 1:2) {
#   for (copula in copulas) {
#     for (n in n_values) {
#       for (tau in taus) {
#         estimates <- replicate(n_replications, {
#           data <- simulate_data(n, copula, tau, mu_T[scenario], sigma_T[scenario], mu_C[scenario], sigma_C[scenario])
#           # Apply estimation here (placeholder, as actual MLE code is needed)
#           # For simplicity, the mean of Y and delta as a proxy
#           mean(data$Y)
#         })
#         
#         # Store results
#         results[[paste(scenario, deparse(substitute(copula)), n, tau, sep = "_")]] <- estimates
#       }
#     }
#   }
# }
# 
# # Summary of results (placeholder)
# results




library(copula)
library(stats4)

# Log-normal parameters for T and C
mu_T <- 2.2
sigma_T <- 1.0
mu_C <- 2.0
sigma_C <- 0.5

# Function to generate lognormal data
gen_lognormal <- function(n, mu, sigma) {
  exp(rnorm(n, mean = mu, sd = sigma))
}

# Simulation settings
n <- 200  # Sample size
tau_vals <- c(0.2, 0.5, 0.7)  # Kendall's tau values for dependence
num_simulations <- 200  # Number of simulations

# Function to generate copula-based dependent data
gen_copula_data <- function(n, copula_family, tau, mu_T, sigma_T, mu_C, sigma_C) {
  # Set copula parameters
  if (copula_family == "Frank") {
    cop <- frankCopula(param = iTau(frankCopula(), tau))
  } else if (copula_family == "Clayton") {
    cop <- claytonCopula(param = iTau(claytonCopula(), tau))
  } else if (copula_family == "Gumbel") {
    cop <- gumbelCopula(param = iTau(gumbelCopula(), tau))
  } else if (copula_family == "Gauss") {
    cop <- normalCopula(param = iTau(normalCopula(), tau))
  }
  
  # Generate copula data
  u <- rCopula(n, cop)
  
  # Generate T and C with lognormal marginals
  T <- qlnorm(u[, 1], meanlog = mu_T, sdlog = sigma_T)
  C <- qlnorm(u[, 2], meanlog = mu_C, sdlog = sigma_C)
  
  return(data.frame(T = T, C = C))
}

# Function to perform MLE
log_likelihood <- function(params, copula_family, data) {
  mu_T <- params[1]
  sigma_T <- params[2]
  mu_C <- params[3]
  sigma_C <- params[4]
  theta <- params[5]
  
  if (copula_family == "Frank") {
    cop <- frankCopula(param = theta)
  } else if (copula_family == "Clayton") {
    cop <- claytonCopula(param = theta)
  } else if (copula_family == "Gumbel") {
    cop <- gumbelCopula(param = theta)
  } else if (copula_family == "Gauss") {
    cop <- normalCopula(param = theta)
  }
  
  u_T <- plnorm(data$T, meanlog = mu_T, sdlog = sigma_T)
  u_C <- plnorm(data$C, meanlog = mu_C, sdlog = sigma_C)
  
  cop_density <- dCopula(cbind(u_T, u_C), cop)
  
  # Calculate log likelihood
  logL <- sum(log(cop_density)) +
    sum(dlnorm(data$T, meanlog = mu_T, sdlog = sigma_T, log = TRUE)) +
    sum(dlnorm(data$C, meanlog = mu_C, sdlog = sigma_C, log = TRUE))
  
  return(-logL)  # Negative log-likelihood for minimization
}

# Running the simulation for different copula families and tau values
results <- list()

for (copula_family in c("Frank", "Clayton", "Gumbel", "Gauss")) {
  for (tau in tau_vals) {
    sim_data <- gen_copula_data(n, copula_family, tau, mu_T, sigma_T, mu_C, sigma_C)
    
    # Initial parameter guesses
    init_params <- c(mu_T, sigma_T, mu_C, sigma_C, 1)
    
    # MLE using optim
    fit <- optim(init_params, log_likelihood, copula_family = copula_family, data = sim_data, 
                 method = "BFGS", hessian = TRUE)
    
    # Store the results
    results[[paste(copula_family, tau, sep = "_")]] <- list(fit = fit, hessian = fit$hessian)
  }
}

# Printing results
results