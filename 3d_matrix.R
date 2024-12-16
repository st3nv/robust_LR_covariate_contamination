# Load necessary libraries
library(MASS)       # For basic linear regression diagnostics
library(robustbase) # For lmrob (robust regression)
library(dplyr)      # For data manipulation

# Simulation parameters for fixed betas
beta_0 <- 5                    # Intercept
p <- 3
beta <- c(2, 3, 1)

sigma <- 1                     # Standard deviation of noise
ctam_sigma_levels <- c(1)
cook_cutoff_levels <- c(4)

n_list <- c(50)
contamination_levels <- c(0.4)
num_runs <- 20                # Number of runs for each parameter combination
num_random_matrices <- 1

# n_list <- c(50, 100, 200, 400)
# contamination_levels <- seq(from = 0, to = 0.6, by = 0.1)
# num_runs <- 250                # Number of runs for each parameter combination
# num_random_matrices <- 200

# Function to generate a random covariance matrix
generate_corr_matrix <- function(p) {
  # Create p x p matrix
  corr_mat <- matrix(0, nrow = p, ncol = p)
  
  # Fill diagonal with 1
  diag(corr_mat) <- 1
  
  # Calculate the number of elements in the lower triangle
  num_lower_tri_elements <- sum(lower.tri(corr_mat))
  
  # Fill lower triangle with random values
  corr_mat[lower.tri(corr_mat)] <- runif(n = num_lower_tri_elements, min = -1, max = 1)
  
  # Fill upper triangle with lower triangle values
  corr_mat[upper.tri(corr_mat)] <- t(corr_mat)[upper.tri(corr_mat)]
  
  return(corr_mat)
}

calculate_mse <- function(beta, real_beta){
  # beta is a k*d matrix, real_beta is a d*1 vector
  return(mean(apply(beta, 1, function(x) mean((x - real_beta)^2))))
}

# Function to generate contaminated data
generate_data <- function(n, sigma, ctam_sigma, beta0, beta, contamination_level, corr_matrix) {
  # X will be 3d random normal with covariance matrix diag(1, 1, 1)
  X_sigma <- matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = X_sigma)
  y <- X %*% beta + beta0 + rnorm(n, sd = sigma)
  n_contaminated <- round(n * contamination_level)
  cov_matrix  <- corr_matrix * ctam_sigma
  if(n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    X[contaminated_indices, ] <- X[contaminated_indices, ] + mvrnorm(n = n_contaminated, mu = rep(0, p), Sigma = cov_matrix)
  }
  X_test <- mvrnorm(n = n, mu = rep(0, p), Sigma = X_sigma)
  y_test <- X_test %*% beta + beta0 + rnorm(n, sd = sigma)
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


results <- data.frame(
  n = numeric(),
  corr_matrix = numeric(),
  contamination_level = numeric(), 
  ctam_sigma = numeric(),
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

is_positive_definite <- function(mat) {
  # Compute eigenvalues
  eigenvalues <- eigen(mat)$values
  
  # Check if all eigenvalues are positive
  return(all(eigenvalues > 0))
}

for (matrix_id in 1:num_random_matrices) {
  success <- FALSE
  while(!success) {
    corr_matrix <- generate_corr_matrix(3)
    # check if the matrix is positive definite
    success <- is_positive_definite(corr_matrix)
  }
  print(corr_matrix)
  for (n in n_list) {
    for (contamination_level in contamination_levels) {
      for (ctam_sigma in ctam_sigma_levels) {
        for (cook_cutoff in cook_cutoff_levels) {
          print (paste("corr_matrix:", corr_matrix, "n:", n, "contamination_level:", contamination_level, "ctam_sigma:", ctam_sigma, "cook_cutoff:", cook_cutoff))
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
            gen_data <- generate_data(n, sigma, ctam_sigma, beta_0, beta, contamination_level, corr_matrix)
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
            # save corr_matrix as a list like [1,2,3,4,5,6,7,8,9] by row
            corr_matrix = paste(corr_matrix, collapse = ","),
            n = n,
            contamination_level = contamination_level,
            ctam_sigma = ctam_sigma,
            cook_cutoff = cook_cutoff,
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




# # Save results to RDS
saveRDS(results, "3d_matrix_results.RDS")
# # saveRDS(run_log, "run_log.RDS")

# # # Save results to CSV, discarding the index column
write.csv(results, "3d_matrix_results.csv", row.names = FALSE)
# # write.csv(run_log, "run_log.csv", row.names = FALSE)


# Show results
print(results)

