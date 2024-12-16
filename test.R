library(MASS)

print(matrix(c(1^2, 0.5^2*1, 0.5^2*1, 1^2), nrow = 2))

test <- mvrnorm(10, mu = c(0, 0), Sigma = matrix(c(1^2, 0.5^2*1, 0.5^2*1, 1^2), nrow = 2))

test