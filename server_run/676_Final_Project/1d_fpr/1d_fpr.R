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
# False positive rate in HT b = 0

# Simulation
# 1. Generate data with 50% of chance of b=0
# 2. Fit OLS, OLS with Cook's distance, M-estimator, MM-estimator, S-estimator
# 3. Calculate false positive rate, true positive rate

# Load necessary libraries
library(MASS)       # For basic linear regression diagnostics
library(robustbase) # For lmrob (robust regression)
library(dplyr)      # For data manipulation
library(sfsmisc)

# Simulation parameters for fixed betas
beta_0 <- 5                    # Intercept
sigma <- 1                     # Standard deviation of noise

# n_list <- c(200)
# contamination_levels <- seq(from = 0, to = 0.5, by = 0.1)
# beta_1_levels <- seq(from = 0, to = 2, by = 0.5)
# ctam_sigma_levels <- c(0, 0.1, 0.5, 1, 2)
# ctam_mu_levels <- c(0, 0.1, 0.5, 1, 2)
# cook_cutoff_levels <- c(4)
# num_runs <- 250                # Number of runs for each paramter combination

n_list <- c(200)
contamination_levels <- c(0.4)
beta_1_levels <- seq(from = 0, to = 1, by = 0.1)
ctam_sigma_levels <- c(0, 0.1, 0.2, 0.5, 1, 2)
ctam_mu_levels <- c(-1, -0.5, 0, 0.5, 1)
cook_cutoff_levels <- c(4)
num_runs <- 200                # Number of runs for each paramter combination

# Function to generate contaminated data
generate_data <- function(n, beta_1, contamination_level, ctam_sigma, ctam_mu){
  x <- rnorm(n)                # Generate real predictor variable
  y <- beta_0 + beta_1 * x + rnorm(n, sd = sigma)  # Response without contamination
  # Apply contamination
  n_contaminated <- round(n * contamination_level)
  if (n_contaminated > 0) {
    contaminated_indices <- sample(1:n, n_contaminated)
    x[contaminated_indices] <- x[contaminated_indices] + rnorm(n_contaminated, mean = ctam_mu, sd = ctam_sigma)
  }
  return(data.frame(x = x, y = y))
}

# Function to perform OLS, robust regression, and OLS with diagnostics (Cook Distance)
run_simulation <- function(data, cook_cutoff) {
  # Ordinary Least Squares (OLS)
  ols_fit <- lm(y ~ x, data = data)
  ols_p_value <- summary(ols_fit)$coefficients[2, 4]
  ols_significance <- ols_p_value < 0.05

  # M-estimator (Robust regression)
  robust_fit_m <-  rlm(y ~ x, data = data, k2 = 1.345)
  robust_m_p_value <- f.robftest(robust_fit_m, var = "x")$p.value
  robust_m_significance <- robust_m_p_value < 0.05
  # robust_coef_m <- coef(robust_fit_m)

  # MM-estimator (Robust regression)
  robust_fit_mm <- lmrob(y ~ x, data = data, method = "MM")
  robust_fit_mm_p_value <- summary(robust_fit_mm)$coefficients[2, 4]
  robust_fit_mm_significance <- robust_fit_mm_p_value < 0.05

  # S-estimator (Robust regression)
  robust_fit_s <- lmrob(y ~ x, data = data, method = "S")
  robust_fit_s_p_value <- summary(robust_fit_s)$coefficients[2, 4]
  robust_fit_s_significance <- robust_fit_s_p_value < 0.05
  # robust_coef_s <- coef(robust_fit_s)

  # OLS with Diagnostics (Cook's distance)
  cooks_distances <- cooks.distance(ols_fit)
  cut_off <- cook_cutoff / nrow(data)
  filtered_data <- data %>% filter(cooks_distances < cut_off)
  ols_diag_fit <- lm(y ~ x, data = filtered_data)
  ols_diag_coef <- coef(ols_diag_fit)
  ols_diag_p_value <- summary(ols_diag_fit)$coefficients[2, 4]
  ols_diag_significance <- ols_diag_p_value < 0.05

  return(list(
    ols_coef = coef(ols_fit),
    robust_fit_m = robust_fit_m,
    robust_fit_mm = robust_fit_mm,
    robust_fit_s = robust_fit_s,
    ols_diag_coef = ols_diag_coef,
    ols_p_value = ols_p_value,
    ols_significance = ols_significance,
    ols_diag_p_value = ols_diag_p_value,
    ols_diag_significance = ols_diag_significance,
    robust_m_p_value = robust_m_p_value,
    robust_m_significance = robust_m_significance,
    robust_fit_mm_p_value = robust_fit_mm_p_value,
    robust_fit_mm_significance = robust_fit_mm_significance,
    robust_fit_s_p_value = robust_fit_s_p_value,
    robust_fit_s_significance = robust_fit_s_significance
  ))
}

# Simulation loop over contamination levels and multiple runs
run_log <- data.frame(
  contamination_level = numeric(),
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  n = numeric(),
  run_id = numeric(),
  cook_cutoff = numeric(),
  y = numeric(),
  x = numeric(),
  ols_coef = numeric(),
  robust_m_coef = numeric(),
  robust_mm_coef = numeric(),
  robust_s_coef = numeric(),
  ols_diag_coef = numeric(),
  ols_p_value = numeric(),
  ols_significance = logical(),
  ols_diag_p_value = numeric(),
  ols_diag_significance = logical(),
  robust_m_p_value = numeric(),
  robust_m_significance = logical(),
  robust_fit_mm_p_value = numeric(),
  robust_fit_mm_significance = logical(),
  robust_fit_s_p_value = numeric(),
  robust_fit_s_significance = logical()
)


results <- data.frame(
  true_beta_1 = numeric(),
  contamination_level = numeric(),
  ctam_sigma = numeric(),
  ctam_mu = numeric(),
  cook_cutoff = numeric(),
  ols_pr = numeric(),
  ols_diag_pr = numeric(),
  robust_m_pr = numeric(),
  robust_mm_pr = numeric(),
  robust_s_pr = numeric()
)

# get the power of the tests so that we can compare the false positive rate and true positive rate at the same time
for (beta_1 in beta_1_levels) {
  for (contamination_level in contamination_levels) {
    for (ctam_sigma in ctam_sigma_levels) {
      for (ctam_mu in ctam_mu_levels) {
        for (n in n_list) {
          for (cook_cutoff in cook_cutoff_levels) {
            # show progress
            print(paste("beta_1:", beta_1, "contamination_level:", contamination_level, "ctam_sigma:", ctam_sigma, "ctam_mu:", ctam_mu, "n:", n, "cook_cutoff:", cook_cutoff))
            ols_num_rejects <- 0
            ols_diag_num_rejects <- 0
            robust_m_num_rejects <- 0
            robust_mm_num_rejects <- 0
            robust_s_num_rejects <- 0
            for (run_id in 1:num_runs) {
              success <- 0
              while(success == 0){ 
                data <- generate_data(n, beta_1, contamination_level, ctam_sigma, ctam_mu)
                results_list <- run_simulation(data, cook_cutoff)
                run_log <- rbind(run_log, data.frame(
                  contamination_level = contamination_level,
                  ctam_sigma = ctam_sigma,
                  ctam_mu = ctam_mu,
                  n = n,
                  run_id = run_id,
                  cook_cutoff = cook_cutoff,
                  y = data$y,
                  x = data$x,
                  ols_coef = results_list$ols_coef,
                  robust_m_coef = coef(results_list$robust_fit_m),
                  robust_mm_coef = coef(results_list$robust_fit_mm),
                  robust_s_coef = coef(results_list$robust_fit_s),
                  ols_diag_coef = results_list$ols_diag_coef,
                  ols_p_value = results_list$ols_p_value,
                  ols_significance = results_list$ols_significance,
                  ols_diag_p_value = results_list$ols_diag_p_value,
                  ols_diag_significance = results_list$ols_diag_significance,
                  robust_m_p_value = results_list$robust_m_p_value,
                  robust_m_significance = results_list$robust_m_significance,
                  robust_fit_mm_p_value = results_list$robust_fit_mm_p_value,
                  robust_fit_mm_significance = results_list$robust_fit_mm_significance,
                  robust_fit_s_p_value = results_list$robust_fit_s_p_value,
                  robust_fit_s_significance = results_list$robust_fit_s_significance
                ))
                # check if all significance are not NA
                if(!is.na(results_list$ols_significance) & !is.na(results_list$ols_diag_significance) & !is.na(results_list$robust_m_significance) & !is.na(results_list$robust_fit_mm_significance) & !is.na(results_list$robust_fit_s_significance)){
                  success <- 1
                }
              }
              ols_num_rejects <- ols_num_rejects + results_list$ols_significance
              ols_diag_num_rejects <- ols_diag_num_rejects + results_list$ols_diag_significance
              robust_m_num_rejects <- robust_m_num_rejects + results_list$robust_m_significance
              robust_mm_num_rejects <- robust_mm_num_rejects + results_list$robust_fit_mm_significance
              robust_s_num_rejects <- robust_s_num_rejects + results_list$robust_fit_s_significance
            }
              
            results <- rbind(results, data.frame(
              true_beta_1 = beta_1,
              contamination_level = contamination_level,
              ctam_sigma = ctam_sigma,
              ctam_mu = ctam_mu,
              cook_cutoff = cook_cutoff,
              ols_pr = ols_num_rejects / num_runs,
              ols_diag_pr = ols_diag_num_rejects / num_runs,
              robust_m_pr = robust_m_num_rejects / num_runs,
              robust_mm_pr = robust_mm_num_rejects / num_runs,
              robust_s_pr = robust_s_num_rejects / num_runs
            ))
          }
        }
      }
    }
  }
}

print(results)

# ## save results to dir called 1d_fpr if not exists create it
saveRDS(results, file = "1d_fpr_results.rds") 
# # save csv
write.csv(results, file = "1d_fpr_results.csv")
