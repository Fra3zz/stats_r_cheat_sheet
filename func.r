
# Function to calculate point estimate
point_estimate <- function(sample_mean) {
  return(sample_mean)
}

# General function to calculate the margin of error for any confidence interval
calculate_margin_of_error <- function(n, sigma_squared, confidence_level = 0.95) { 
  if (n <= 0 || sigma_squared <= 0 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Invalid input values. Ensure n > 0, sigma_squared > 0, and 0 < confidence_level < 1.")
  }
  sigma <- sqrt(sigma_squared)
  critical_value <- qnorm(1 - (1 - confidence_level) / 2)
  margin_of_error <- critical_value * (sigma / sqrt(n))
  return(margin_of_error)
}

# Function to calculate sigma squared
calculate_sigma_squared <- function(population_std_dev) {
  return(population_std_dev ^ 2)
}

# Function to calculate confidence interval using Z-distribution
confidence_interval_z <- function(sample_mean, population_std_dev, n, confidence_level = 0.95) {
  if (n <= 0 || population_std_dev <= 0 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Invalid input values. Ensure n > 0, population_std_dev > 0, and 0 < confidence_level < 1.")
  }
  me <- calculate_margin_of_error(n, population_std_dev^2, confidence_level)
  lower_bound <- sample_mean - me
  upper_bound <- sample_mean + me
  return(list(lower_limit = lower_bound, upper_limit = upper_bound, margin_of_error = me))
}

# Function to calculate the width of a confidence interval
confidence_interval_width <- function(n, population_std_dev, confidence_level) {
  if (n <= 0 || population_std_dev <= 0 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Invalid input values. Ensure n > 0, population_std_dev > 0, and 0 < confidence_level < 1.")
  }
  std_error <- population_std_dev / sqrt(n)
  critical_value <- qnorm(1 - (1 - confidence_level) / 2)
  margin_of_error <- critical_value * std_error
  width <- 2 * margin_of_error
  return(width)
}

# Function to calculate confidence interval using t-distribution
confidence_interval_t <- function(sample_mean, sample_std_dev, n, confidence_level) {
  if (n <= 1 || sample_std_dev <= 0 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Invalid input values. Ensure n > 1, sample_std_dev > 0, and 0 < confidence_level < 1.")
  }
  df <- n - 1
  std_error <- sample_std_dev / sqrt(n)
  critical_value <- qt(1 - (1 - confidence_level) / 2, df)
  margin_of_error <- critical_value * std_error
  lower_limit <- sample_mean - margin_of_error
  upper_limit <- sample_mean + margin_of_error
  return(list(lower_limit = lower_limit, upper_limit = upper_limit, margin_of_error = margin_of_error, critical_value = critical_value))
}

# Master function to calculate various statistics
calculate_statistics <- function(data, confidence_level = NULL, population_std_dev = NULL, p = NULL, n_trials = NULL) {
  # Validate input data
  if (!is.numeric(data) || length(data) == 0) {
    stop("Data must be a non-empty numeric vector.")
  }
  
  # Calculate basic statistics
  mean_value <- mean(data)
  median_value <- median(data)
  
  # Calculate mode manually
  mode_value <- as.numeric(names(sort(table(data), decreasing = TRUE)[1]))
  
  range_value <- diff(range(data))
  sd_value <- sd(data)
  
  # Prepare results list
  results <- list(
    mean = mean_value,
    median = median_value,
    mode = mode_value,
    range = range_value,
    standard_deviation = sd_value
  )
  
  # Calculate confidence interval if confidence level and population standard deviation are provided
  if (!is.null(confidence_level) && !is.null(population_std_dev)) {
    n <- length(data)
    std_error <- population_std_dev / sqrt(n)
    critical_value <- qnorm(1 - (1 - confidence_level) / 2)
    margin_of_error <- critical_value * std_error
    lower_limit <- mean_value - margin_of_error
    upper_limit <- mean_value + margin_of_error
    
    results$confidence_interval <- list(
      lower_limit = lower_limit,
      upper_limit = upper_limit,
      margin_of_error = margin_of_error
    )
  }
  
  # Calculate normal probability if population standard deviation is provided
  if (!is.null(population_std_dev)) {
    z_scores <- (data - mean_value) / population_std_dev
    normal_probabilities <- pnorm(z_scores)
    results$normal_probabilities <- normal_probabilities
  }
  
  # Calculate binomial probabilities if p and n_trials are provided
  if (!is.null(p) && !is.null(n_trials)) {
    binomial_probabilities <- dbinom(data, size = n_trials, prob = p)
    results$binomial_probabilities <- binomial_probabilities
  }
  
  return(results)
}



calculate_t_critical_value_With_only_confidenceLevel_and_n <- function(confidence_level, n) {
  # Validate inputs
  if (n <= 1 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Invalid input values. Ensure n > 1 and 0 < confidence_level < 1.")
  }
  
  # Calculate degrees of freedom
  df <- n - 1
  
  # Calculate t-critical value
  t_critical_value <- qt(1 - (1 - confidence_level) / 2, df)
  
  # Return the t-critical value rounded to two decimal places
  return(t_critical_value)
}



#p-test test statistic
make_test_statistic <- function(meanX, h0, s, number_of_samplesS){
  return(((mean - h0)/(s/sqrt(number_of_samplesS)))) # nolint
}

make_p_value <- function(t_score, degrees_of_freedom, lower_tail) {
  # Calculate p-value using t-distribution
  result <- pt(q = t_score, df = degrees_of_freedom, lower.tail = lower_tail)
  return(result)
}
