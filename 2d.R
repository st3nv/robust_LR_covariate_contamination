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
beta_0 <- 5                    # Intercept
beta_1 <- 2                    # Slope1
beta_2 <- 3                    # Slope2
sigma <- 1                     # Standard deviation of noise

n_list <- c(50)
contamination_levels <- c(0.4)
ctam_sigma_levels <- c(1)
ctam_phos_levels <- c(0.5)
ctam_mu_levels <- c(0)
cook_cutoff_levels <- c(4)
num_runs <- 10                # Number of runs for each paramter combination

# n_list <- c(50, 100, 200)
# contamination_levels <- seq(from = 0, to = 0.6, by = 0.1)
# ctam_sigma_levels <- c(1)
# ctam_phos_levels <- seq(from = -1, to = 1, by = 0.2)
# ctam_mu_levels <- c(0)
# cook_cutoff_levels <- c(4)
# num_runs <- 500                # Number of runs for each paramter combination

# Function to generate contaminated data
generate_data <- function(n, contamination_level, ctam_sigma, ctam_mu, ctam_rho){
  x1 <- rnorm(n)
  x2 <- rnorm(n)             # Generate real predictor variable
  x1_test <- rnorm(n)
  x2_test <- rnorm(n)
  y <- beta_0 + beta_1 * x1 + beta_2 * x2 + rnorm(n, sd = sigma)  # Response without contamination
  y_test <- beta_0 + beta_1 * x1_test + beta_2 * x2_test + rnorm(n, sd = sigma)
  x1_real <- x1
  x2_real <- x2
  
  # Apply contamination
  Sigma <- matrix(c(ctam_sigma^2, ctam_sigma^2*ctam_rho, ctam_sigma^2*ctam_rho, ctam_sigma^2), nrow = 2)
  n_contaminated <- round(n * contamination_level)
  if (n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    contamination <- mvrnorm(n_contaminated, mu = c(0, 0), Sigma = Sigma)
    x1[contaminated_indices] <- x1[contaminated_indices] + contamination[, 1]
    x2[contaminated_indices] <- x2[contaminated_indices] + contamination[, 2]
  }
  return(data.frame(x1 = x1, x2 = x2, y = y, x1_real = x1_real, x2_real = x2_real, y_test = y_test, x1_test = x1_test, x2_test = x2_test))
}

calculate_mse <- function(beta, real_beta){
  # beta is a k*d matrix, real_beta is a d*1 vector
  return(mean(apply(beta, 1, function(x) mean((x - real_beta)^2))))
}

# Function to perform OLS, robust regression, and OLS with diagnostics (Cook Distance)
run_simulation <- function(data, cook_cutoff) {
  # Ordinary Least Squares (OLS)
  ols_fit <- lm(y ~ x1 + x2, data = data)
  ols_coef <- coef(ols_fit)

  # M-estimator (Robust regression)
  robust_fit_m <- rlm(y ~ x1 + x2, data = data, k2 = 1.345)
  # robust_coef_m <- coef(robust_fit_m)

  # MM-estimator (Robust regression)
  robust_fit_mm <- lmrob(y ~ x1 + x2, data = data, method = "MM")
  # robust_coef_mm <- coef(robust_fit_mm)

  # S-estimator (Robust regression)
  robust_fit_s <- lmrob(y ~ x1 + x2, data = data, method = "S")
  # robust_coef_s <- coef(robust_fit_s)

  # OLS with Diagnostics (Cook's distance)
  cooks_distances <- cooks.distance(ols_fit)
  cut_off <- 4 / nrow(data)
  filtered_data <- data %>% filter(cooks_distances < cut_off)
  ols_diag_fit <- lm(y ~ x1 + x2, data = filtered_data)
  ols_diag_coef <- coef(ols_diag_fit)

  return(list(ols_coef = ols_coef, ols_diag_coef = ols_diag_coef, robust_fit_m = robust_fit_m, robust_fit_mm = robust_fit_mm, robust_fit_s = robust_fit_s))
}

# Simulation loop over contamination levels and multiple runs
run_log <- data.frame()


results <- data.frame(
  contamination_level = numeric(), 
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  cook_cutoff = numeric(),
  ctam_rho = numeric(),  # Add ctam_rho to results
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

for (n in n_list) {
  for (contamination_level in contamination_levels) {
    for (ctam_sigma in ctam_sigma_levels) {
      for (ctam_mu in ctam_mu_levels) {
        for (cook_cutoff in cook_cutoff_levels) {
          for (ctam_rho in ctam_phos_levels) {
            # show a progress message
            print(paste("n:", n, "contamination_level:", contamination_level, "ctam_sigma:", ctam_sigma, "ctam_mu:", ctam_mu, "cook_cutoff:", cook_cutoff, "ctam_rho:", ctam_rho))
            # Define metrics
            robust_m_est <- matrix(NA, nrow = num_runs, ncol = 3)
            robust_mm_est <- matrix(NA, nrow = num_runs, ncol = 3)
            robust_s_est <- matrix(NA, nrow = num_runs, ncol = 3)
            ols_est <- matrix(NA, nrow = num_runs, ncol = 3)
            ols_diag_est <- matrix(NA, nrow = num_runs, ncol = 3)
            robust_m_pred_squared_error <- numeric(num_runs)
            robust_mm_pred_squared_error <- numeric(num_runs)
            robust_s_pred_squared_error <- numeric(num_runs)
            ols_pred_squared_error <- numeric(num_runs)
            ols_diag_pred_squared_error <- numeric(num_runs)
            for (run in 1:num_runs) {
              data <- generate_data(n, contamination_level, ctam_sigma, ctam_mu, ctam_rho)
              results_list <- run_simulation(data, cook_cutoff)
              run_log_row <- data.frame(
                contamination_level = contamination_level, 
                ctam_sigma = ctam_sigma,
                ctam_mu = ctam_mu,
                n = n,
                run_id = run,
                cook_cutoff = cook_cutoff,
                ctam_rho = ctam_rho,  # Add ctam_rho to run_log_row
                # save x_real, y, x, y_test, x_test as string of list like [1,2,3,...]
                x1_real = paste(data$x1_real, collapse = ","),
                x2_real = paste(data$x2_real, collapse = ","),
                y = paste(data$y, collapse = ","),
                x1 = paste(data$x1, collapse = ","),
                x2 = paste(data$x2, collapse = ","),
                y_test = paste(data$y_test, collapse = ","),
                x1_test = paste(data$x1_test, collapse = ","),
                x2_test = paste(data$x2_test, collapse = ","),
                ols_coef = paste(results_list$ols_coef, collapse = ","),
                robust_m_coef = paste(results_list$robust_fit_m$coefficients, collapse = ","),
                robust_mm_coef = paste(results_list$robust_fit_mm$coefficients, collapse = ","),
                robust_s_coef = paste(results_list$robust_fit_s$coefficients, collapse = ","),
                ols_diag_coef = paste(results_list$ols_diag_coef, collapse = ","),
                ols_pred_squared_error = mean((data$y_test - data$x1_test * results_list$ols_coef[2] - data$x2_test * results_list$ols_coef[3] - results_list$ols_coef[1])^2),
                robust_m_pred_squared_error = mean((data$y_test - data$x1_test * results_list$robust_fit_m$coefficients[2] - data$x2_test * results_list$robust_fit_m$coefficients[3] - results_list$robust_fit_m$coefficients[1])^2),
                robust_mm_pred_squared_error = mean((data$y_test - data$x1_test * results_list$robust_fit_mm$coefficients[2] - data$x2_test * results_list$robust_fit_mm$coefficients[3] - results_list$robust_fit_mm$coefficients[1])^2),
                robust_s_pred_squared_error = mean((data$y_test - data$x1_test * results_list$robust_fit_s$coefficients[2] - data$x2_test * results_list$robust_fit_s$coefficients[3] - results_list$robust_fit_s$coefficients[1])^2),
                ols_diag_pred_squared_error = mean((data$y_test - data$x1_test * results_list$ols_diag_coef[2] - data$x2_test * results_list$ols_diag_coef[3] - results_list$ols_diag_coef[1])^2)
              )
              run_log <- rbind(run_log, run_log_row)

              # Calculate metrics
              robust_m_est[run, ] <- results_list$robust_fit_m$coefficients
              robust_mm_est[run, ] <- results_list$robust_fit_mm$coefficients
              robust_s_est[run, ] <- results_list$robust_fit_s$coefficients
              ols_est[run, ] <- results_list$ols_coef
              ols_diag_est[run, ] <- results_list$ols_diag_coef
              robust_m_pred_squared_error[run] <- run_log_row$robust_m_pred_squared_error
              robust_mm_pred_squared_error[run] <- run_log_row$robust_mm_pred_squared_error
              robust_s_pred_squared_error[run] <- run_log_row$robust_s_pred_squared_error
              ols_pred_squared_error[run] <- run_log_row$ols_pred_squared_error
              ols_diag_pred_squared_error[run] <- run_log_row$ols_diag_pred_squared_error
            }
            # Calculate metrics
            robust_m_est_bias <- colMeans(robust_m_est[, 2:3]) - c(beta_1, beta_2)
            robust_m_est_variance <- apply(robust_m_est[, 2:3], 2, var)
            robust_m_est_mse <- calculate_mse(robust_m_est[, 2:3], c(beta_1, beta_2))
            robust_mm_est_bias <- colMeans(robust_mm_est[, 2:3]) - c(beta_1, beta_2)
            robust_mm_est_variance <- apply(robust_mm_est[, 2:3], 2, var)
            robust_mm_est_mse <- calculate_mse(robust_mm_est[, 2:3], c(beta_1, beta_2))
            robust_s_est_bias <- colMeans(robust_s_est[, 2:3]) - c(beta_1, beta_2)
            robust_s_est_variance <- apply(robust_s_est[, 2:3], 2, var)
            robust_s_est_mse <- calculate_mse(robust_s_est[, 2:3], c(beta_1, beta_2))
            ols_est_bias <- colMeans(ols_est[, 2:3]) - c(beta_1, beta_2)
            ols_est_variance <- apply(ols_est[, 2:3], 2, var)
            ols_est_mse <- calculate_mse(ols_est[, 2:3], c(beta_1, beta_2))
            ols_diag_est_bias <- colMeans(ols_diag_est[, 2:3]) - c(beta_1, beta_2)
            ols_diag_est_variance <- apply(ols_diag_est[, 2:3], 2, var)
            ols_diag_est_mse <- calculate_mse(ols_diag_est[, 2:3], c(beta_1, beta_2))
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
              ctam_rho = ctam_rho,  # Add ctam_rho to results_row
              robust_m_est_bias = mean(robust_m_est_bias^2),
              robust_m_est_variance = mean(robust_m_est_variance),
              robust_m_est_mse = robust_m_est_mse,
              robust_mm_est_bias = mean(robust_mm_est_bias^2),
              robust_mm_est_variance = mean(robust_mm_est_variance),
              robust_mm_est_mse = robust_mm_est_mse,
              robust_s_est_bias = mean(robust_s_est_bias^2),
              robust_s_est_variance = mean(robust_s_est_variance),
              robust_s_est_mse = robust_s_est_mse,
              ols_est_bias = mean(ols_est_bias^2),
              ols_est_variance = mean(ols_est_variance),
              ols_est_mse = ols_est_mse,
              ols_diag_est_bias = mean(ols_diag_est_bias^2),
              ols_diag_est_variance = mean(ols_diag_est_variance),
              ols_diag_est_mse = ols_diag_est_mse,
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
}


# # Save results to RDS
# saveRDS(results, "results.RDS")
# saveRDS(run_log, "run_log.RDS")

# # Save results to CSV, discarding the index column
# write.csv(results, "results.csv", row.names = FALSE)
# write.csv(run_log, "run_log.csv", row.names = FALSE)


# # Show results
print(results)
