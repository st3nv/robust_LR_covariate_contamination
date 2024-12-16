# Load necessary libraries
library(MASS)       # For basic linear regression diagnostics
library(robustbase) # For robust regression
library(dplyr)      # For data manipulation

# Function to generate contaminated data for multi-dimensional predictors
generate_multidim_data <- function(n, p, contamination_level, ctam_sigma, ctam_mu) {
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)  # Generate predictors
  x_test <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  beta <- c(5, rep(2, p))  # Set coefficients: intercept = 5, slopes = 2 for each predictor
  y <- 5 + rowSums(x * 2) + rnorm(n, sd = 1)  # Response without contamination
  y_test <- 5 + rowSums(x_test * 2) + rnorm(n, sd = 1)
  x_real <- x
  
  # Apply contamination
  n_contaminated <- round(n * contamination_level)
  if (n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    x[contaminated_indices, ] <- x[contaminated_indices, ] + 
      matrix(rnorm(n_contaminated * p, mean = ctam_mu, sd = ctam_sigma), 
             nrow = n_contaminated, ncol = p)
  }
  
  return(list(x = x, y = y, x_real = x_real, y_test = y_test, x_test = x_test))
}

# Function to measure runtime for each method
measure_runtime_multidim <- function(data) {
  timings <- list()
  df <- data.frame(y = data$y, data$x)  # Convert to data.frame for model fitting
  
  # Measure OLS runtime
  start_time <- Sys.time()
  ols_fit <- lm(y ~ ., data = df)
  timings$ols <- Sys.time() - start_time
  
  # Measure OLS with diagnostics (Cook's Distance)
  start_time <- Sys.time()
  cooks_distances <- cooks.distance(ols_fit)
  cut_off <- 4 / nrow(df)
  filtered_data <- df[cooks_distances < cut_off, ]
  lm(y ~ ., data = filtered_data)
  timings$ols_diag <- Sys.time() - start_time
  
  # Measure M-estimator runtime
  start_time <- Sys.time()
  rlm(y ~ ., data = df, k2 = 1.345)
  timings$m_estimator <- Sys.time() - start_time
  
  # Measure MM-estimator runtime
  start_time <- Sys.time()
  lmrob(y ~ ., data = df, method = "MM")
  timings$mm_estimator <- Sys.time() - start_time
  
  # Measure S-estimator runtime
  start_time <- Sys.time()
  lmrob(y ~ ., data = df, method = "S")
  timings$s_estimator <- Sys.time() - start_time
  
  return(timings)
}

# Parameters for simulation
n <- 1000
p_list <- c(1, 5, 10, 20, 50)
contamination_levels <- seq(0, 0.6, by = 0.2)
ctam_sigma_levels <- c(1)
ctam_mu_levels <- c(0)
num_runs <- 40  # Reduced for quick testing

# Results to store runtimes
runtime_results <- data.frame(
  p = numeric(),
  contamination_level = numeric(),
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  ols_avg_runtime = numeric(),
  ols_diag_avg_runtime = numeric(),
  m_estimator_avg_runtime = numeric(),
  mm_estimator_avg_runtime = numeric(),
  s_estimator_avg_runtime = numeric()
)

# Simulation loop
for (p in p_list) {
  for (contamination_level in contamination_levels) {
    for (ctam_sigma in ctam_sigma_levels) {
      for (ctam_mu in ctam_mu_levels) {
        print(paste("Running simulation for p =", p, ", contamination =", contamination_level, 
                    ", ctam_sigma =", ctam_sigma, ", ctam_mu =", ctam_mu))
        # Collect runtimes
        ols_times <- c()
        ols_diag_times <- c()
        m_estimator_times <- c()
        mm_estimator_times <- c()
        s_estimator_times <- c()
        
        for (run in 1:num_runs) {
          data <- generate_multidim_data(n, p, contamination_level, ctam_sigma, ctam_mu)
          timings <- measure_runtime_multidim(data)
          
          ols_times <- c(ols_times, as.numeric(timings$ols, units = "secs"))
          ols_diag_times <- c(ols_diag_times, as.numeric(timings$ols_diag, units = "secs"))
          m_estimator_times <- c(m_estimator_times, as.numeric(timings$m_estimator, units = "secs"))
          mm_estimator_times <- c(mm_estimator_times, as.numeric(timings$mm_estimator, units = "secs"))
          s_estimator_times <- c(s_estimator_times, as.numeric(timings$s_estimator, units = "secs"))
        }
        
        # Calculate average runtimes
        runtime_results <- rbind(runtime_results, data.frame(
          p = p,
          contamination_level = contamination_level,
          ctam_sigma = ctam_sigma,
          ctam_mu = ctam_mu,
          ols_avg_runtime = mean(ols_times),
          ols_diag_avg_runtime = mean(ols_diag_times),
          m_estimator_avg_runtime = mean(m_estimator_times),
          mm_estimator_avg_runtime = mean(mm_estimator_times),
          s_estimator_avg_runtime = mean(s_estimator_times)
        ))
      }
    }
  }
}

# Save the results to CSV
# write.csv(runtime_results, "runtime_comparison_multidim.csv", row.names = FALSE)

# Display runtime results
print(runtime_results)
