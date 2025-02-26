# Function to calculate point estimate
point_estimate <- function(sample_mean, population_std_dev, n) {
  return(sample_mean)
}

# Function to calculate 95% margin of error
margin_of_error <- function(sample_mean, population_std_dev, n) {
  # Calculate standard error
  std_error <- population_std_dev / sqrt(n)
  
  # Calculate critical value for 95% confidence interval (Z-score for normal distribution)
  critical_value <- qnorm(0.975)
  
  # Calculate margin of error
  me <- critical_value * std_error
  
  return(round(me, 3))
}

# Function to calculate the margin of error for a confidence interval
calculate_margin_of_error <- function(n, sigma_squared, confidence_level) {
  # Check if input values are valid
  if (n <= 0 || sigma_squared <= 0 || confidence_level < 0 || confidence_level > 1) {
    stop("Invalid input values. Please check the following: n > 0, sigma_squared > 0, 0 < confidence_level < 1")
  }

  # Calculate population standard deviation
  sigma <- sqrt(sigma_squared)

  # Calculate critical value from standard normal distribution
  critical_value <- qnorm(1 - (1 - confidence_level) / 2)

  # Calculate margin of error
  margin_of_error <- critical_value * (sigma / sqrt(n))

  # Return the result
  return(margin_of_error)
}