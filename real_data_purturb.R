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

# Function to contaminate only the training set during cross-validation
contaminate_data <- function(data) {
  n <- nrow(data)
  p <- ncol(data) - 1  # Exclude target variable
  
  # Identify numeric and categorical columns
  numeric_cols <- sapply(data, is.numeric)
  categorical_cols <- !numeric_cols
  
  # Randomly select 10% of cells to contaminate
  num_contaminated <- ceiling(0.1 * n * sum(numeric_cols))
  rows <- sample(1:n, size = num_contaminated, replace = TRUE)
  cols <- sample(which(numeric_cols), size = num_contaminated, replace = TRUE)
  
  for (i in seq_along(rows)) {
    data[rows[i], cols[i]] <- data[rows[i], cols[i]] + rnorm(1, mean = 0, sd = 2)
  }
  
  # Contaminate categorical columns
  for (col in which(categorical_cols)) {
    contaminated_rows <- sample(1:n, size = ceiling(0.1 * n), replace = FALSE)
    unique_vals <- unique(data[[col]])
    for (row in contaminated_rows) {
      data[row, col] <- sample(unique_vals[unique_vals != data[row, col]], 1)
    }
  }
  
  return(data)
}

# Define a function to perform 5-fold cross-validation
perform_cv <- function(data) {
  # Shuffle the dataset
  data <- data[sample(nrow(data)), ]
  
  # Create 5 folds
  folds <- cut(seq(1, nrow(data)), breaks = 5, labels = FALSE)
  
  # Initialize vectors to store MSPE for each model
  ols_mspe <- numeric(5)
  ols_diag_mspe <- numeric(5)
  robust_mspe_s <- numeric(5)
  m_mspe <- numeric(5)
  robust_mspe_mm <- numeric(5)
  
  # Perform 5-fold CV
  for (i in 1:5) {
    # Split data into training and testing sets
    test_indices <- which(folds == i, arr.ind = TRUE)
    test_data <- data[test_indices, ]
    train_data <- data[-test_indices, ]
    
    # Contaminate the training data only
    train_data <- contaminate_data(train_data)
    
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

# Perform 5-fold CV on original data and display results
cv_results <- perform_cv(data)

cat("5-Fold CV MSPE Results (Original Data):\n")
cat("OLS MSPE:", cv_results$ols, "\n")
cat("OLS with Diagnosis MSPE:", cv_results$ols_diag, "\n")
cat("Robust S-Regression MSPE:", cv_results$robust_s, "\n")
cat("M-Regression MSPE:", cv_results$m, "\n")
cat("MM-Regression MSPE:", cv_results$robust_mm, "\n")
