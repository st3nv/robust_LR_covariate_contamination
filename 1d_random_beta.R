### Simple linear regression with covairate contamination
## Given (X_cta, y)
## y = X b + e, X is standard normal
## X_cta = X + t
## Assume t independent of X 

# Parameters
# Var(t), E(t)

# Estimators
# OLS: b = (X^T X)^{-1} X^T y
# OLS: with Cook's distance diagnostic
# Robust: M-estimator with huber loss 
# Example fitH <- rlm(income ~ education, data = Duncan, k2 = 1.345) 
# Robust: MM-estimator
# robust_fit_mm <- lmrob(income ~ education, data = Duncan, method = "MM")
# Robust: S-estimator
# fitS <- lqs(income ~ education, data = Duncan, method = "S")



# Metrics
# \hat{b}: bias/variance
# MSPE
# False positive rate in HT b = 0

# install.packages("MASS")
# install.packages("robustbase")
# install.packages("dplyr")

# Load necessary libraries
library(MASS)       # For basic linear regression diagnostics
library(robustbase) # For lmrob (robust regression)
library(dplyr)      # For data manipulation

# Simulation parameters for fixed betas
beta_0 <- 5                    # Intercept
beta_1 <- 2                    # Slope
sigma <- 1                     # Standard deviation of noise

n_list <- c(100)
contamination_levels <- c(0.1)
ctam_sigma_levels <- c(0.2)
ctam_mu_levels <- c(0)
cook_cutoff_levels <- c(4)
num_runs <- 100

# n_list <- c(50, 100, 200, 500, 1000)
# contamination_levels <- seq(from = 0, to = 0.6, by = 0.1)
# ctam_sigma_levels <- c(0.1, 0.2, 0.5, 1, 2)
# ctam_mu_levels <- c(-1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1)
# cook_cutoff_levels <- c(1, 2, 3, 4)
# num_runs <- 250                # Number of runs for each paramter combination

# Function to generate contaminated data
generate_data <- function(n, contamination_level, ctam_sigma, ctam_mu){
  x <- rnorm(n)                # Generate real predictor variable
  x_test <- rnorm(n)
  y <- beta_0 + beta_1 * x + rnorm(n, sd = sigma)  # Response without contamination
  y_test <- beta_0 + beta_1 * x_test + rnorm(n, sd = sigma)
  x_real <- x
  
  # Apply contamination
  n_contaminated <- round(n * contamination_level)
  if (n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    x[contaminated_indices] <- x[contaminated_indices] + rnorm(n_contaminated, mean = ctam_mu, sd = ctam_sigma)
  }
  return(data.frame(x = x, y = y, x_real = x_real, y_test = y_test, x_test = x_test))
}

# Function to perform OLS, robust regression, and OLS with diagnostics (Cook Distance)
run_simulation <- function(data, cook_cutoff) {
  # Ordinary Least Squares (OLS)
  ols_fit <- lm(y ~ x, data = data)
  ols_coef <- coef(ols_fit)

  # M-estimator (Robust regression)
  robust_fit_m <-  rlm(y ~ x, data = data, k2 = 1.345) 
  # robust_coef_m <- coef(robust_fit_m)

  # MM-estimator (Robust regression)
  robust_fit_mm <- lmrob(y ~ x, data = data, method = "MM")
  # robust_coef_mm <- coef(robust_fit_mm)

  # S-estimator (Robust regression)
  robust_fit_s <- lmrob(y ~ x, data = data, method = "S")
  # robust_coef_s <- coef(robust_fit_s)

  # OLS with Diagnostics (Cook's distance)
  cooks_distances <- cooks.distance(ols_fit)
  cut_off <- cook_cutoff / nrow(data)
  filtered_data <- data %>% filter(cooks_distances < cut_off)
  ols_diag_fit <- lm(y ~ x, data = filtered_data)
  ols_diag_coef <- coef(ols_diag_fit)

  return(list(ols_coef = ols_coef, ols_diag_coef = ols_diag_coef, robust_fit_m = robust_fit_m, robust_fit_mm = robust_fit_mm, robust_fit_s = robust_fit_s))
}

# Simulation loop over contamination levels and multiple runs
run_log <- data.frame(
  contamination_level = numeric(), 
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  n = numeric(),
  run_id = numeric(),
  cook_cutoff = numeric(),
  x_real = numeric(),
  y = numeric(),
  x = numeric(),
  y_test = numeric(),
  x_test = numeric(),
  ols_coef = numeric(),
  robust_m_coef = numeric(),
  robust_mm_coef = numeric(),
  robust_s_coef = numeric(),
  ols_diag_coef = numeric()
)


results <- data.frame(
  contamination_level = numeric(), 
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  cook_cutoff = numeric(),
  robust_m_est_bias = numeric(),
  robust_m_est_variance = numeric(),
  robust_mm_est_bias = numeric(),
  robust_mm_est_variance = numeric(),
  robust_s_est_bias = numeric(),
  robust_s_est_variance = numeric(),
  ols_est_bias = numeric(),
  ols_est_variance = numeric(),
  ols_diag_est_bias = numeric(),
  ols_diag_est_variance = numeric(),
  ols_mspe = numeric(),
  robust_m_mspe = numeric(),
  robust_mm_mspe = numeric(),
  robust_s_mspe = numeric(),
  ols_diag_mspe = numeric()
)

for (n in n_list) {
  for (contamination_level in contamination_levels) {
    for (ctam_sigma in ctam_sigma_levels) {
      for (ctam_mu in ctam_mu_levels) {
        for (cook_cutoff in cook_cutoff_levels) {
          # show a progress message
          print(paste("n:", n, "contamination_level:", contamination_level, "ctam_sigma:", ctam_sigma, "ctam_mu:", ctam_mu, "cook_cutoff:", cook_cutoff))
          # Define metrics
          robust_m_est <- c()
          robust_mm_est <- c()
          robust_s_est <- c()
          ols_est <- c()
          ols_diag_est <- c()
          robust_m_pred_squared_error <- c()
          robust_mm_pred_squared_error <- c()
          robust_s_pred_squared_error <- c()
          ols_pred_squared_error <- c()
          ols_diag_pred_squared_error <- c()
          for (run in 1:num_runs) {
            data <- generate_data(n, contamination_level, ctam_sigma, ctam_mu)
            results_list <- run_simulation(data, cook_cutoff)
            run_log_row <- data.frame(
              contamination_level = contamination_level, 
              ctam_sigma = ctam_sigma,
              ctam_mu = ctam_mu,
              n = n,
              run_id = run,
              cook_cutoff = cook_cutoff,
              # save x_real, y, x, y_test, x_test as string of list like [1,2,3,...]
              x_real = paste(data$x_real, collapse = ","),
              y = paste(data$y, collapse = ","),
              x = paste(data$x, collapse = ","),
              y_test = paste(data$y_test, collapse = ","),
              x_test = paste(data$x_test, collapse = ","),
              ols_coef = paste(results_list$ols_coef, collapse = ","),
              robust_m_coef = paste(results_list$robust_fit_m$coefficients, collapse = ","),
              robust_mm_coef = paste(results_list$robust_fit_mm$coefficients, collapse = ","),
              robust_s_coef = paste(results_list$robust_fit_s$coefficients, collapse = ","),
              ols_diag_coef = paste(results_list$ols_diag_coef, collapse = ","),
              ols_pred_squared_error = mean((data$y_test - data$x_test * results_list$ols_coef[2] - results_list$ols_coef[1])^ 2),
              robust_m_pred_squared_error = mean((data$y_test - data$x_test * results_list$robust_fit_m$coefficients[2] - results_list$robust_fit_m$coefficients[1])^ 2),
              robust_mm_pred_squared_error = mean((data$y_test - data$x_test * results_list$robust_fit_mm$coefficients[2] - results_list$robust_fit_mm$coefficients[1])^ 2),
              robust_s_pred_squared_error = mean((data$y_test - data$x_test * results_list$robust_fit_s$coefficients[2] - results_list$robust_fit_s$coefficients[1])^ 2),
              ols_diag_pred_squared_error = mean((data$y_test - data$x_test * results_list$ols_diag_coef[2] - results_list$ols_diag_coef[1])^ 2)
            )
            run_log <- rbind(run_log, run_log_row)

            # Calculate metrics
            robust_m_est <- c(robust_m_est, results_list$robust_fit_m$coefficients[2])
            robust_mm_est <- c(robust_mm_est, results_list$robust_fit_mm$coefficients[2])
            robust_s_est <- c(robust_s_est, results_list$robust_fit_s$coefficients[2])
            ols_est <- c(ols_est, results_list$ols_coef[2])
            ols_diag_est <- c(ols_diag_est, results_list$ols_diag_coef[2])
            robust_m_pred_squared_error <- c(robust_m_pred_squared_error, run_log_row$robust_m_pred_squared_error)
            robust_mm_pred_squared_error <- c(robust_mm_pred_squared_error, run_log_row$robust_mm_pred_squared_error)
            robust_s_pred_squared_error <- c(robust_s_pred_squared_error, run_log_row$robust_s_pred_squared_error)
            ols_pred_squared_error <- c(ols_pred_squared_error, run_log_row$ols_pred_squared_error)
            ols_diag_pred_squared_error <- c(ols_diag_pred_squared_error, run_log_row$ols_diag_pred_squared_error)
          }
          # Calculate metrics
          robust_m_est_bias <- mean(robust_m_est) - beta_1
          robust_m_est_variance <- var(robust_m_est)
          robust_mm_est_bias <- mean(robust_mm_est) - beta_1
          robust_mm_est_variance <- var(robust_mm_est)
          robust_s_est_bias <- mean(robust_s_est) - beta_1
          robust_s_est_variance <- var(robust_s_est)
          ols_est_bias <- mean(ols_est) - beta_1
          ols_est_variance <- var(ols_est)
          ols_diag_est_bias <- mean(ols_diag_est) - beta_1
          ols_diag_est_variance <- var(ols_diag_est)
          ols_mspe <- mean(ols_pred_squared_error)
          robust_m_mspe <- mean(robust_m_pred_squared_error)
          robust_mm_mspe <- mean(robust_mm_pred_squared_error)
          robust_s_mspe <- mean(robust_s_pred_squared_error)
          ols_diag_mspe <- mean(ols_diag_pred_squared_error)

          # Save results
          results_row <- data.frame(
            contamination_level = contamination_level, 
            ctam_sigma = ctam_sigma,
            ctam_mu = ctam_mu,
            cook_cutoff = cook_cutoff,
            robust_m_est_bias = robust_m_est_bias,
            robust_m_est_variance = robust_m_est_variance,
            robust_mm_est_bias = robust_mm_est_bias,
            robust_mm_est_variance = robust_mm_est_variance,
            robust_s_est_bias = robust_s_est_bias,
            robust_s_est_variance = robust_s_est_variance,
            ols_est_bias = ols_est_bias,
            ols_est_variance = ols_est_variance,
            ols_diag_est_bias = ols_diag_est_bias,
            ols_diag_est_variance = ols_diag_est_variance,
            ols_mspe = ols_mspe,
            robust_m_mspe = robust_m_mspe,
            robust_mm_mspe = robust_mm_mspe,
            robust_s_mspe = robust_s_mspe,
            ols_diag_mspe = ols_diag_mspe
          )
          results <- rbind(results, results_row)
        }
      }
    }
  }
}


# # Save results to RDS
# saveRDS(results, "1d_fixed_beta/results.RDS")
# saveRDS(run_log, "1d_fixed_beta/run_log.RDS")
# 
# # Save results to CSV, discarding the index column
# write.csv(results, "1d_fixed_beta/results.csv", row.names = FALSE)
# write.csv(run_log, "1d_fixed_beta/run_log.csv", row.names = FALSE)


# Show results
print(results)
