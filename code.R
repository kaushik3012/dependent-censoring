library(copula)

mu_T <- c(2.2, 2.5)  # Mean parameters for T
sigma_T <- c(1.0, 1.0)  # Standard deviation parameters for T
mu_C <- c(2.0, 2.0)  # Mean parameters for C
sigma_C <- c(0.25, 0.50)  # Standard deviation parameters for C
tau <- 0.2
params <- c(mu_T[2], sigma_T[2], mu_C[2], sigma_C[2], tau)

# Function to simulate data for a given copula and sample size
simulate_data <- function(n, params, copula_family) {
  mu_T <- params[1]
  sigma_T <- params[2]
  mu_C <- params[3]
  sigma_C <- params[4]
  tau <- params[5]
  if (copula_family == "Frank") {
    cop <- frankCopula()
  } else if (copula_family == "Clayton") {
    cop <- claytonCopula()
  } else if (copula_family == "Gumbel") {
    cop <- gumbelCopula()
  } else if (copula_family == "Gauss") {
    cop <- normalCopula()
  }
  cop <- setTheta(cop, iTau(cop, tau))
  u <- rCopula(n, cop)
  T <- qlnorm(u[,1], meanlog = mu_T, sdlog = sigma_T)
  C <- qlnorm(u[,2], meanlog = mu_C, sdlog = sigma_C)
  Y <- pmin(T, C)
  delta <- as.numeric(T <= C)
  return(data.frame(Y, delta))
}


# ... [Previous code for mu_T, sigma_T, mu_C, sigma_C, tau, params, and simulate_data function remains unchanged]

likelihood_of_y <- function(x, params, copula_family) {
  mu_T <- params[1]
  sigma_T <- params[2]
  mu_C <- params[3]
  sigma_C <- params[4]
  tau <- params[5]
  
  if (copula_family == "Frank") {
    cop <- frankCopula(iTau(frankCopula(), tau))
  } else if (copula_family == "Clayton") {
    cop <- claytonCopula(iTau(claytonCopula(), tau))
  } else if (copula_family == "Gumbel") {
    cop <- gumbelCopula(iTau(gumbelCopula(), tau))
  } else if (copula_family == "Gauss") {
    cop <- normalCopula(iTau(normalCopula(), tau))
  }
  
  F_T <- plnorm(x[1], meanlog = mu_T, sdlog = sigma_T)
  F_C <- plnorm(x[1], meanlog = mu_C, sdlog = sigma_C)
  
  if (x[2] == 1) {
    f_T <- dlnorm(x[1], meanlog = mu_T, sdlog = sigma_T)
    cond_cop <- cCopula(cbind(F_T, F_C), indices=2, cop)[1]
    ans <- log(f_T) + log(1-cond_cop)
  } else {
    f_C <- dlnorm(x[1], meanlog = mu_C, sdlog = sigma_C)
    cond_cop <- cCopula(cbind(F_C, F_T), indices=2, cop)[1]
    ans <- log(f_C) + log(1-cond_cop)
  }
  
  return(ifelse(is.finite(ans), ans, -1e10))  # Return a large negative value if ans is not finite
}

# Function to perform MLE
log_likelihood <- function(params, copula_family, data) {
  if (any(params <= 0) || params[5] <= -1 || params[5] >= 1) {
    return(1e10)  # Return a large positive value for invalid parameters
  }
  
  log_likelihoods <- apply(data, 1, likelihood_of_y, params=params, copula_family=copula_family)
  logL <- sum(log_likelihoods)
  
  return(ifelse(is.finite(logL), -logL, 1e10))  # Return a large positive value if logL is not finite
}

# Simulate and estimate parameters
n_simulations <- 2  # Number of simulations
n <- 200  # Sample size for each simulation
init_params <- c(1, 1.0, 1, 1, 0.5)  # Initial parameters

# Store results
estimates <- matrix(NA, nrow=n_simulations, ncol=length(init_params))
asderr <- matrix(NA, nrow=n_simulations, ncol=length(init_params))

for (i in 1:n_simulations) {
  data <- simulate_data(n, params, copula_family = "Clayton")
  
  # Estimate parameters using MLE
  fit <- try(optim(par=init_params, fn=log_likelihood, copula_family = "Clayton", data = data,
                   method = "L-BFGS-B", hessian = TRUE, 
                   lower=c(0.01,0.01,0.01,0.01,-0.8), upper=c(10,10,10,10,0.8), 
                   control=list(maxit=100)), silent=TRUE)
  
  if (!inherits(fit, "try-error") && fit$convergence == 0) {
    estimates[i, ] <- fit$par
    
    # Calculate asymptotic standard errors from the Hessian matrix
    if (!is.null(fit$hessian) && !any(is.na(fit$hessian)) && det(fit$hessian) != 0) {
      asderr[i, ] <- sqrt(diag(solve(fit$hessian)))
    } else {
      asderr[i, ] <- rep(NA, length(init_params))
    }
  } else {
    estimates[i, ] <- rep(NA, length(init_params))
    asderr[i, ] <- rep(NA, length(init_params))
  }
}

# Calculate metrics
aver.est <- colMeans(estimates, na.rm=TRUE)
sd.aver.est <- apply(estimates, 2, sd, na.rm=TRUE)
aver.asderr <- colMeans(asderr, na.rm=TRUE)
RMSE <- sqrt(colMeans((estimates - matrix(rep(params, n_simulations), nrow=n_simulations, byrow=TRUE))^2, na.rm=TRUE))

# Output the results
results <- list(
  aver.est = aver.est,
  sd.aver.est = sd.aver.est,
  aver.asderr = aver.asderr,
  RMSE = RMSE
)

print(results)