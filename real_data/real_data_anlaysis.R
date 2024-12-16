# Load required libraries
library(dplyr)
library(quantreg)
library(sfsmisc)
library(robustbase)  # For lmrob
library(MASS)        # For rlm

# Read and prepare the data
data <- read.csv("/Users/songtengyu/Documents/2024Fall/Linear Models/final project/real_data/ObesityDataSet.csv")
data <- data %>% select(-c(NObeyesdad, CALC))  # Remove the target column if needed

# Set seed for reproducibility
set.seed(42)

cv <- 5

# Define a function to perform k-fold cross-validation
perform_cv <- function(data) {
  # Shuffle the dataset
  data <- data[sample(nrow(data)), ]
  
  # Create folds
  folds <- cut(seq(1, nrow(data)), breaks = cv, labels = FALSE)
  
  # Initialize vectors to store MSPE for each model
  ols_mspe <- numeric(cv)
  ols_diag_mspe <- numeric(cv)
  robust_mspe_s <- numeric(cv)
  m_mspe <- numeric(cv)
  robust_mspe_mm <- numeric(cv)
  
  # Perform CV
  for (i in 1:cv) {
    # Split data into training and testing sets
    test_indices <- which(folds == i, arr.ind = TRUE)
    test_data <- data[test_indices, ]
    train_data <- data[-test_indices, ]
    
    # Fit OLS model
    ols_fit <- lm(Weight ~ ., data = train_data)
    
    # Fit OLS with Cook's distance diagnostics
    cooks_distances <- cooks.distance(ols_fit)
    cut_off <- 4 / nrow(train_data)
    filtered_data <- train_data %>% filter(cooks_distances < cut_off)
    ols_diag_fit <- lm(Weight ~ ., data = filtered_data)
    
    # Fit robust regression models
    # Robust S-regression
    robust_fit_s <- lmrob(Weight ~ ., data = train_data, method = "S", setting = "KS2014")
    
    # M-regression
    m_fit <- rlm(Weight ~ ., data = train_data, method = "M")
    
    # MM-regression
    robust_fit_mm <- lmrob(Weight ~ ., data = train_data, method = "MM")
    
    # Predict on the testing set
    ols_pred <- predict(ols_fit, newdata = test_data)
    ols_diag_pred <- predict(ols_diag_fit, newdata = test_data)
    robust_pred_s <- predict(robust_fit_s, newdata = test_data)
    m_pred <- predict(m_fit, newdata = test_data)
    robust_pred_mm <- predict(robust_fit_mm, newdata = test_data)
    
    # Calculate MSPE for each model
    ols_mspe[i] <- mean((test_data$Weight - ols_pred)^2)
    ols_diag_mspe[i] <- mean((test_data$Weight - ols_diag_pred)^2)
    robust_mspe_s[i] <- mean((test_data$Weight - robust_pred_s)^2)
    m_mspe[i] <- mean((test_data$Weight - m_pred)^2)
    robust_mspe_mm[i] <- mean((test_data$Weight - robust_pred_mm)^2)
  }
  
  # Return average MSPE for each model
  return(list(
    ols = mean(ols_mspe),
    ols_diag = mean(ols_diag_mspe),
    robust_s = mean(robust_mspe_s),
    m = mean(m_mspe),
    robust_mm = mean(robust_mspe_mm)
  ))
}

# Perform k-fold CV and display results
cv_results <- perform_cv(data)

cat("5-Fold CV MSPE Results:\n")
cat("OLS MSPE:", cv_results$ols, "\n")
cat("OLS with Diagnosis MSPE:", cv_results$ols_diag, "\n")
cat("Robust S-Regression MSPE:", cv_results$robust_s, "\n")
cat("M-Regression MSPE:", cv_results$m, "\n")
cat("MM-Regression MSPE:", cv_results$robust_mm, "\n")
