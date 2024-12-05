library(copula)
library(knitr)

# c(0.2, 0.5, 0.7)
# c(200, 500)

# Define parameters
mu_T <- 2.2
sigma_T <- 1.0
mu_C <- 2.0
sigma_C <- 0.25
tau_values <- c(0.2)
n_values <- c(500)
n_simulations <- 20  # Number of simulations for each combination
available_copulas <- c("Frank", "Clayton", "Gumbel", "Gauss")
copula_family <- available_copulas[1]

# Function to simulate data (remains unchanged)
simulate_data <- function(n, params, copula_family) {
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
  
  u <- rCopula(n, cop)
  
  T <- qlnorm(u[,1], meanlog = mu_T, sdlog = sigma_T)
  C <- qlnorm(u[,2], meanlog = mu_C, sdlog = sigma_C)
  
  Y <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  return(data.frame(Y, delta))
}

# Likelihood function (remains unchanged)
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
    ans <- log(f_T) + log(1 - cond_cop)
  } else {
    f_C <- dlnorm(x[1], meanlog = mu_C, sdlog = sigma_C)
    cond_cop <- cCopula(cbind(F_C, F_T), indices=2, cop)[1]
    ans <- log(f_C) + log(1 - cond_cop)
  }
  
  return(ifelse(is.finite(ans), ans, -1e10))
}

# Log-likelihood function (remains unchanged)
log_likelihood <- function(params, copula_family, data) {
  if (any(params <= 0) || abs(params[5]) >= 1) {
    return(1e10)
  }
  
  log_likelihoods <- apply(data, 1, likelihood_of_y, params=params, copula_family=copula_family)
  logL <- sum(log_likelihoods)
  
  return(ifelse(is.finite(logL), -logL, 1e10))
}

# Function to run simulations and compute metrics
run_simulations <- function(n, tau, n_simulations) {
  params <- c(mu_T, sigma_T, mu_C, sigma_C, tau)
  # init_params <- c(0.5, 0.5, 0.5, 0.5, 0.5)
  init_params <- params
  estimates <- matrix(NA, nrow=n_simulations, ncol=length(init_params))
  asderr <- matrix(NA, nrow=n_simulations, ncol=length(init_params))
  
  for (i in 1:n_simulations) {
    data <- simulate_data(n, params, copula_family = copula_family)
    fit <- try(optim(par=init_params, fn=log_likelihood, copula_family = copula_family, data = data,
                     method = "L-BFGS-B", hessian = TRUE, 
                     lower=c(0.01,0.01,0.01,0.01,-0.8), upper=c(10,10,10,10,0.8), 
                     control=list(maxit=100, reltol=1e-8)))
    
    if (!inherits(fit, "try-error") && fit$convergence == 0) {
      estimates[i, ] <- fit$par
      
      if (!is.null(fit$hessian) && !any(is.na(fit$hessian)) && det(fit$hessian) != 0) {
        asderr[i, ] <- sqrt(diag(solve(fit$hessian)))
      }
    }
    
    cat("Completed", i, "simulations for n =", n, "and tau =", tau, "\n")
  }
  
  aver.est <- colMeans(estimates, na.rm=TRUE)
  sd.aver.est <- apply(estimates, 2, sd, na.rm=TRUE)
  aver.asderr <- colMeans(asderr, na.rm=TRUE)
  RMSE <- sqrt(colMeans((estimates - matrix(rep(params, n_simulations), nrow=n_simulations, byrow=TRUE))^2, na.rm=TRUE))
  
  return(list(aver.est=aver.est, sd.aver.est=sd.aver.est, aver.asderr=aver.asderr, RMSE=RMSE))
}

# Run simulations and store results
results <- list()
for (n in n_values) {
  for (tau in tau_values) {
    key <- paste(n, tau, sep="_")
    results[[key]] <- run_simulations(n, tau, n_simulations)
  }
}

# Function to create a formatted table row
create_row <- function(tau, n, results, metric) {
  key <- paste(n, tau, sep="_")
  res <- results[[key]][[metric]]
  c(tau, formatC(res[1:4], digits=2, format="f"), formatC(res[5], digits=2, format="f"), formatC(res[5], digits=2, format="f"))
}

# Create table
table_data <- do.call(rbind, lapply(n_values, function(n) {
  do.call(rbind, lapply(tau_values, function(tau) {
    rbind(
      create_row(tau, n, results, "aver.est"),
      create_row(tau, n, results, "sd.aver.est"),
      create_row(tau, n, results, "aver.asderr"),
      create_row(tau, n, results, "RMSE")
    )
  }))
}))

# Add row names
rownames(table_data) <- rep(c("aver.est", "sd.aver.est", "aver.asderr", "RMSE"), length(n_values) * length(tau_values))

# Create and print the table
table <- kable(table_data, 
               col.names = c("τ", "μ_T", "σ_T", "μ_C", "σ_C", "θ", "τ"),
               caption = paste("Simulation results for the",copula_family,"copula"),
               align = 'c')

print(table)

### Example Output
"
Table: Simulation results for the Clayton copula

|            |  τ  | μ_T  | σ_T  | μ_C  | σ_C  |  θ   |  τ   |
  |:-----------|:---:|:----:|:----:|:----:|:----:|:----:|:----:|
  |aver.est    | 0.2 | 2.19 | 1.00 | 1.99 | 0.24 | 0.19 | 0.19 |
  |sd.aver.est | 0.2 | 0.01 | 0.02 | 0.06 | 0.00 | 0.24 | 0.24 |
  |aver.asderr | 0.2 | 0.11 | 0.09 | 0.05 | 0.02 | 0.22 | 0.22 |
  |RMSE        | 0.2 | 0.01 | 0.01 | 0.04 | 0.01 | 0.17 | 0.17 |
  |aver.est    | 0.5 | 2.27 | 1.13 | 1.99 | 0.23 | 0.47 | 0.47 |
  |sd.aver.est | 0.5 | 0.06 | 0.08 | 0.03 | 0.00 | 0.11 | 0.11 |
  |aver.asderr | 0.5 | 0.12 | 0.11 | 0.04 | 0.02 | 0.14 | 0.14 |
  |RMSE        | 0.5 | 0.09 | 0.14 | 0.02 | 0.02 | 0.09 | 0.09 |
  |aver.est    | 0.7 | 2.10 | 0.97 | 1.98 | 0.25 | 0.74 | 0.74 |
  |sd.aver.est | 0.7 | 0.14 | 0.14 | 0.05 | 0.02 | 0.09 | 0.09 |
  |aver.asderr | 0.7 | 0.09 | 0.08 | 0.03 | 0.02 | 0.06 | 0.06 |
  |RMSE        | 0.7 | 0.14 | 0.11 | 0.04 | 0.02 | 0.07 | 0.07 |
  |aver.est    | 0.2 | 2.28 | 1.04 | 2.03 | 0.25 | 0.12 | 0.12 |
  |sd.aver.est | 0.2 | 0.06 | 0.04 | 0.02 | 0.02 | 0.16 | 0.16 |
  |aver.asderr | 0.2 | 0.06 | 0.05 | 0.04 | 0.02 | 0.15 | 0.15 |
  |RMSE        | 0.2 | 0.09 | 0.05 | 0.04 | 0.02 | 0.14 | 0.14 |
  |aver.est    | 0.5 | 2.20 | 1.02 | 2.03 | 0.24 | 0.37 | 0.37 |
  |sd.aver.est | 0.5 | 0.02 | 0.01 | 0.02 | 0.02 | 0.11 | 0.11 |
  |aver.asderr | 0.5 | 0.07 | 0.06 | 0.03 | 0.02 | 0.13 | 0.13 |
  |RMSE        | 0.5 | 0.01 | 0.03 | 0.03 | 0.02 | 0.15 | 0.15 |
  |aver.est    | 0.7 | 2.17 | 0.94 | 1.99 | 0.27 | 0.74 | 0.74 |
  |sd.aver.est | 0.7 | 0.03 | 0.03 | 0.01 | 0.01 | 0.02 | 0.02 |
  |aver.asderr | 0.7 | 0.06 | 0.05 | 0.02 | 0.01 | 0.04 | 0.04 |
  |RMSE        | 0.7 | 0.04 | 0.07 | 0.01 | 0.02 | 0.04 | 0.04 |"

