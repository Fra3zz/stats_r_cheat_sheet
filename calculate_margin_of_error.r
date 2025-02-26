# Given values
n <- 30  # sample size
sigma_squared <- 3.6  # population variance
confidence_level <- 0.95  # confidence level for the interval

# Calculate population standard deviation
sigma <- sqrt(sigma_squared)

# Calculate critical value from standard normal distribution
critical_value <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate margin of error
margin_of_error <- critical_value * (sigma / sqrt(n))

# Round to three decimal places
margin_of_error <- round(margin_of_error, 3)

# Print the result
print(margin_of_error)