# Load necessary libraries
library(MASS)       # For basic linear regression diagnostics
library(robustbase) # For lmrob (robust regression)
library(dplyr)      # For data manipulation
library(ggplot2)    # For plotting
library(quantreg)
library(sfsmisc)

# Simulation parameters for fixed betas
beta_0 <- 5                    # Intercept
beta_1 <- 2                    # Slope
sigma <- 1                     # Standard deviation of noise

set.seed(42)
n <- 100
contamination_level <- 0.2
ctam_sigma_level <- 2
ctam_mu_level <- 1

# Function to generate contaminated data
generate_data <- function(n, contamination_level, ctam_sigma, ctam_mu){
  x <- rnorm(n)                # Generate real predictor variable
  x_test <- rnorm(n)
  y <- beta_0 + beta_1 * x + rnorm(n, sd = sigma)  # Response without contamination
  y_test <- beta_0 + beta_1 * x_test + rnorm(n, sd = sigma)
  x_real <- x
  contaminated <- rep(FALSE, n)
  # Apply contamination
  n_contaminated <- round(n * contamination_level)
  if (n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    x[contaminated_indices] <- x[contaminated_indices] + rnorm(n_contaminated, mean = ctam_mu, sd = ctam_sigma)
    contaminated[contaminated_indices] <- TRUE
  }
  return(data.frame(x = x, y = y, contaminated=contaminated ))
}

# function to plot the fitted line and data
plot_data <- function(data, ols_coef, ols_diag_coef, robust_coef_m, robust_coef_mm, robust_coef_s) {   
  ggplot(data, aes(x = x, y = y)) + 
    geom_point(aes(color = ifelse(contaminated, "Contaminated", "Clean")), size = 2) + 
    geom_abline(aes(intercept = ols_coef[1], slope = ols_coef[2], color = "OLS")) + 
    geom_abline(aes(intercept = ols_diag_coef[1], slope = ols_diag_coef[2], color = "OLS-Diag")) + 
    geom_abline(aes(intercept = robust_coef_m[1], slope = robust_coef_m[2], color = "M-Estimator")) + 
    geom_abline(aes(intercept = robust_coef_mm[1], slope = robust_coef_mm[2], color = "MM-Estimator")) + 
    geom_abline(aes(intercept = robust_coef_s[1], slope = robust_coef_s[2], color = "S-Estimator")) + 
    theme_minimal() + 
    labs(x = "x",  
         y = "y") + 
    theme(legend.position = "right") +
    scale_color_manual(name = "Legend", 
                       values = c("Contaminated" = "brown", 
                                  "Clean" = "black", 
                                  "OLS" = "blue", 
                                  "OLS-Diag" = "green", 
                                  "M-Estimator" = "red", 
                                  "MM-Estimator" = "purple", 
                                  "S-Estimator" = "orange"))
}


data <- generate_data(n, contamination_level, ctam_sigma_level, ctam_mu_level)

# OLS
ols_fit <- lm(y ~ x, data = data)
ols_coef <- coef(ols_fit)

# OLS diag
cooks_distances <- cooks.distance(ols_fit)
cut_off <- 4 / nrow(data)
filtered_data <- data %>% filter(cooks_distances < cut_off)
ols_diag_fit <- lm(y ~ x, data = filtered_data)
ols_diag_coef <- coef(ols_diag_fit)

# M-estimator (Robust regression)
robust_fit_m <- rlm(y ~ x, data = data, k2 = 1.345) 
robust_coef_m <- coef(robust_fit_m)

# MM-estimator (Robust regression)
robust_fit_mm <- lmrob(y ~ x, data = data, method = "MM")
robust_coef_mm <- coef(robust_fit_mm)

# S-estimator (Robust regression)
robust_fit_s <- lqs(y ~ x, data = data, method = "S")
robust_coef_s <- coef(robust_fit_s)

# Plot
plot_data(data, ols_coef, ols_diag_coef, robust_coef_m, robust_coef_mm, robust_coef_s)

