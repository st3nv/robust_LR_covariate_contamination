### 2D covariate contamination simulation
## Now there are two covariates, x1 and x2, and the contamination is on both covariates

# Parameters
# mu: vector of means for the contamination, Sigma: covariance matrix for the contamination

# Estimators
# OLS, OLS with diagnostics (Cook's distance), M-estimator, MM-estimator, S-estimator

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
sigma <- 1                     # Standard deviation of noise

# Parameters
p_n_list <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
n_list <- c(200)
contamination_levels <- c(0, 0.1, 0.2, 0.4, 0.5, 0.6)
ctam_sigma_levels <- c(1)
ctam_mu_levels <- c(0)
cook_cutoff_levels <- c(4)
num_runs <- 500                # Number of runs for each paramter combination

# p_n_list <- c(0.1)
# n_list <- c(50)
# contamination_levels <- c(0.4)
# ctam_sigma_levels <- c(1)
# ctam_mu_levels <- c(0)
# cook_cutoff_levels <- c(4)
# num_runs <- 250

# Function to generate contaminated data (p dimensional, noise uncorrelated)
generate_data <- function(n, p, beta, contamination_level, ctam_sigma, ctam_mu){
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X %*% beta + rnorm(n, sd = sigma)  # Response without contamination
  # Apply contamination
  n_contaminated <- round(n * contamination_level)
  if(n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    X[contaminated_indices, ] <- X[contaminated_indices, ] + matrix(rnorm(n_contaminated * p, mean = ctam_mu, sd = ctam_sigma), nrow = n_contaminated)
  }
  X_test <- matrix(rnorm(n * p), nrow = n)
  y_test <- X_test %*% beta + rnorm(n, sd = sigma)
  return (list(data.frame(y, X), data.frame(y_test, X_test)))
}

# Function to perform OLS, robust regression, and OLS with diagnostics (Cook Distance)
run_simulation <- function(data, cook_cutoff) {
  # Ordinary Least Squares (OLS)
  ols_fit <- lm(y ~ ., data = data)
  # M-estimator (Robust regression)
  robust_fit_m <- rlm(y ~ ., data = data, k2 = 1.345)
  # robust_coef_m <- coef(robust_fit_m)

  # MM-estimator (Robust regression)
  robust_fit_mm <- lmrob(y ~ ., data = data, method = "MM")
  # robust_coef_mm <- coef(robust_fit_mm)

  # S-estimator (Robust regression)
  robust_fit_s <- lmrob(y ~ ., data = data, method = "S")
  # robust_coef_s <- coef(robust_fit_s)

  # OLS with Diagnostics (Cook's distance)
  cooks_distances <- cooks.distance(ols_fit)
  cut_off <- 4 / nrow(data)
  filtered_data <- data %>% filter(cooks_distances < cut_off)
  ols_diag_fit <- lm(y ~ ., data = filtered_data)

  return(list(
    ols_fit = ols_fit,
    robust_fit_m = robust_fit_m,
    robust_fit_mm = robust_fit_mm,
    robust_fit_s = robust_fit_s,
    ols_diag_fit = ols_diag_fit
  ))
}

# Simulation loop over contamination levels and multiple runs
run_log <- data.frame()

results <- data.frame(
  n = numeric(),
  p_n = numeric(),
  p = numeric(),
  contamination_level = numeric(), 
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  cook_cutoff = numeric(),
  robust_m_est_bias = numeric(),
  robust_m_est_variance = numeric(),
  robust_m_est_mse = numeric(),
  robust_mm_est_bias = numeric(),
  robust_mm_est_variance = numeric(),
  robust_mm_est_mse = numeric(),
  robust_s_est_bias = numeric(),
  robust_s_est_variance = numeric(),
  robust_s_est_mse = numeric(),
  ols_est_bias = numeric(),
  ols_est_variance = numeric(),
  ols_est_mse = numeric(),
  ols_diag_est_bias = numeric(),
  ols_diag_est_variance = numeric(),
  ols_diag_est_mse = numeric(),
  ols_mspe = numeric(),
  robust_m_mspe = numeric(),
  robust_mm_mspe = numeric(),
  robust_s_mspe = numeric(),
  ols_diag_mspe = numeric()
)

calculate_mse <- function(beta, real_beta){
  # beta is a k*d matrix, real_beta is a d*1 vector
  return(mean(apply(beta, 1, function(x) mean((x - real_beta)^2))))
}

for (n in n_list) {
  for (p_n in p_n_list) {
    p <- round(n * p_n)
    for (contamination_level in contamination_levels) {
      for (ctam_sigma in ctam_sigma_levels) {
        for (ctam_mu in ctam_mu_levels) {
          for (cook_cutoff in cook_cutoff_levels) {
            print (paste("n:", n, "p_n:", p_n, "p:", p, "contamination_level:", contamination_level, "ctam_sigma:", ctam_sigma, "ctam_mu:", ctam_mu, "cook_cutoff:", cook_cutoff))
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
            for (run_id in 1:num_runs) {
              beta <- rnorm(p)
              # Generate data unpack from list
              gen_data <- generate_data(n, p, beta, contamination_level, ctam_sigma, ctam_mu)
              train_data <- gen_data[[1]]
              test_data <- gen_data[[2]]
              # Run simulation
              simulation_results <- run_simulation(train_data, cook_cutoff)
              # Extract coefficients
              ols_coef <- coef(simulation_results$ols_fit)[-1]
              robust_m_coef <- coef(simulation_results$robust_fit_m)[-1]
              robust_mm_coef <- coef(simulation_results$robust_fit_mm)[-1]
              robust_s_coef <- coef(simulation_results$robust_fit_s)[-1]
              ols_diag_coef <- coef(simulation_results$ols_diag_fit)[-1]
              # Calculate bias and variance
              robust_m_est <- rbind(robust_m_est, robust_m_coef)
              robust_mm_est <- rbind(robust_mm_est, robust_mm_coef)
              robust_s_est <- rbind(robust_s_est, robust_s_coef)
              ols_est <- rbind(ols_est, ols_coef)
              ols_diag_est <- rbind(ols_diag_est, ols_diag_coef)
              # Calculate MSPE
              y_test <- as.matrix(test_data$y_test)
              # X_test is from the 2nd column in the test_data
              X_test <- as.matrix(test_data[,-1])
              X_test_pred <- cbind(1, X_test)
              
              # Calculate predictions error
              robust_m_pred_squared_error <- c(robust_m_pred_squared_error, mean((y_test - X_test_pred %*% as.vector(simulation_results$robust_fit_m$coef))^2))
              robust_mm_pred_squared_error <- c(robust_mm_pred_squared_error, mean((y_test - X_test_pred %*% as.vector(simulation_results$robust_fit_mm$coef))^2))
              robust_s_pred_squared_error <- c(robust_s_pred_squared_error, mean((y_test - X_test_pred %*% as.vector(simulation_results$robust_fit_s$coef))^2))
              ols_pred_squared_error <- c(ols_pred_squared_error, mean((y_test - X_test_pred %*% as.vector(simulation_results$ols_fit$coef))^2))
              ols_diag_pred_squared_error <- c(ols_diag_pred_squared_error, mean((y_test - X_test_pred %*% as.vector(simulation_results$ols_diag_fit$coef))^2))
             
            }
            results <- rbind(results, data.frame(
              n = n,
              p_n = p_n,
              p = p,
              contamination_level = contamination_level,
              ctam_sigma = ctam_sigma,
              ctam_mu = ctam_mu,
              cook_cutoff = cook_cutoff,
              # first mean over all runs, subtract true beta, then mean over all coefficients
              # bias squared
              robust_m_est_bias = mean((colMeans(robust_m_est) - beta)^2),
              robust_m_est_variance = mean(apply(robust_m_est, 2, var)),
              robust_m_est_mse = calculate_mse(robust_m_est, beta),
              robust_mm_est_bias = mean((colMeans(robust_mm_est) - beta)^2),
              robust_mm_est_variance = mean(apply(robust_mm_est, 2, var)),
              robust_mm_est_mse = calculate_mse(robust_mm_est, beta),
              robust_s_est_bias = mean((colMeans(robust_s_est) - beta)^2),
              robust_s_est_variance = mean(apply(robust_s_est, 2, var)),
              robust_s_est_mse = calculate_mse(robust_s_est, beta),
              ols_est_bias = mean((colMeans(ols_est) - beta)^2),
              ols_est_variance = mean(apply(ols_est, 2, var)),
              ols_est_mse = calculate_mse(ols_est, beta),
              ols_diag_est_bias = mean((colMeans(ols_diag_est) - beta)^2),
              ols_diag_est_variance = mean(apply(ols_diag_est, 2, var)),
              ols_diag_est_mse = calculate_mse(ols_diag_est, beta),
              ols_mspe = mean(ols_pred_squared_error),
              robust_m_mspe = mean(robust_m_pred_squared_error),
              robust_mm_mspe = mean(robust_mm_pred_squared_error),
              robust_s_mspe = mean(robust_s_pred_squared_error),
              ols_diag_mspe = mean(ols_diag_pred_squared_error)
            ))
          }
        }
      }
    }
  }
}


# # Save results to RDS
saveRDS(results, "nd_fixed_results.RDS")
# saveRDS(run_log, "run_log.RDS")

# # Save results to CSV, discarding the index column
write.csv(results, "nd_fixed_results.csv", row.names = FALSE)
# write.csv(run_log, "run_log.csv", row.names = FALSE)


# # Show results
print(results)
