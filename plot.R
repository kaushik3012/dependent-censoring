library(copula)

# Define parameters for Scenario 1
mu_T <- 2.2
sigma_T <- 1.0
mu_C <- 2.0
sigma_C <- 0.25

f_Y <- function(x, copula, tau, mu_T, sigma_T, mu_C, sigma_C) {
  
  cop <- setTheta(copula, iTau(copula, tau))
  
  F_T <- plnorm(x, meanlog = mu_T, sdlog = sigma_T)
  F_C <- plnorm(x, meanlog = mu_C, sdlog = sigma_C)
  
  # if (x[2] == 1) {
    f_T <- dlnorm(x, meanlog = mu_T, sdlog = sigma_T)
    cond_cop <- cCopula(cbind(F_T, F_C), indices=2, cop)[1]
    # if (is.nan(cond_cop) || 1-cond_cop <= 0) return(0)
    ans <- f_T*(1 - cond_cop)
  # } else {
    f_C <- dlnorm(x[1], meanlog = mu_C, sdlog = sigma_C)
    cond_cop <- cCopula(cbind(F_C, F_T), indices=2, cop)[1]
    # if (is.nan(cond_cop) || 1-cond_cop <= 0) return(0)
    ans <- ans + f_C *(1 - cond_cop)
  # }
  
  return(ans)
}


# Function to calculate survival function of Y
S_Y <- function(y, copula, tau, mu_T, sigma_T, mu_C, sigma_C) {
  cop <- setTheta(copula, iTau(copula, tau))
  F_T <- plnorm(y, meanlog = mu_T, sdlog = sigma_T)
  F_C <- plnorm(y, meanlog = mu_C, sdlog = sigma_C)
  
  F_Y <- F_T+F_C - pCopula(cbind(F_T,F_C), cop)
  return(1-F_Y)
}

# Function to calculate hazard function of Y
h_Y <- function(y, copula, tau, mu_T, sigma_T, mu_C, sigma_C) {
  f_Y(y, copula, tau, mu_T, sigma_T, mu_C, sigma_C)/S_Y(y, copula, tau, mu_T, sigma_T, mu_C, sigma_C)
}

# Create plots
plot_copulas <- function(tau_values) {
  copulas <- list(
    Clayton = claytonCopula(),
    Frank = frankCopula(),
    Gumbel = gumbelCopula(),
    Gauss = normalCopula()
  )
  
  colors <- c("blue", "red", "green", "cyan")
  line_types <- c(2, 3, 4, 5)
  
  y_values <- seq(0, 30, length.out = 100)
  
  for (tau in tau_values) {
    # Density plot
    plot(NULL, xlim = c(0, 30), ylim = c(0, 0.2), 
         main = paste("tau =", tau), xlab = "y", ylab = "f_Y")
    for (i in seq_along(copulas)) {
      lines(y_values, sapply(y_values, f_Y, copula = copulas[[i]], tau = tau, 
                             mu_T = mu_T, sigma_T = sigma_T, mu_C = mu_C, sigma_C = sigma_C),
            col = colors[i], lty = line_types[i])
    }
    legend("topright", names(copulas), col = colors, lty = line_types, cex = 0.8)
    
    # Survival function plot
    plot(NULL, xlim = c(0, 30), ylim = c(0, 1), 
         main = paste("tau =", tau), xlab = "y", ylab = "S_Y")
    for (i in seq_along(copulas)) {
      lines(y_values, sapply(y_values, S_Y, copula = copulas[[i]], tau = tau, 
                             mu_T = mu_T, sigma_T = sigma_T, mu_C = mu_C, sigma_C = sigma_C),
            col = colors[i], lty = line_types[i])
    }
    legend("topright", names(copulas), col = colors, lty = line_types, cex = 0.8)
    
    # Hazard function plot
    plot(NULL, xlim = c(0, 30), ylim = c(0, 1), 
         main = paste("tau =", tau), xlab = "y", ylab = "h_Y")
    for (i in seq_along(copulas)) {
      lines(y_values, sapply(y_values, h_Y, copula = copulas[[i]], tau = tau, 
                             mu_T = mu_T, sigma_T = sigma_T, mu_C = mu_C, sigma_C = sigma_C),
            col = colors[i], lty = line_types[i])
    }
    legend("topright", names(copulas), col = colors, lty = line_types, cex = 0.8)
  }
}

# Set up the plotting area
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

# Create plots for three τ values
tau_values <- c(0.2, 0.5, 0.7)
plot_copulas(tau_values)