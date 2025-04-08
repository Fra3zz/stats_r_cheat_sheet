
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



# Function to calculate the t-test statistic
calculate_t_statistic <- function(sample_mean, hypothesized_mean, standard_deviation, sample_size) {
  # Calculate the t-statistic using the formula: (sample_mean - hypothesized_mean) / (standard_deviation / sqrt(sample_size))
  return(((sample_mean - hypothesized_mean) / (standard_deviation / sqrt(sample_size))))
}

# Function to calculate the p-value for a t-test
calculate_p_value <- function(t_score, degrees_of_freedom, is_lower_tail) {
  # Calculate p-value using t-distribution
  result <- pt(q = t_score, df = degrees_of_freedom, lower.tail = is_lower_tail)
  return(result)
}


# Function to perform a one-sample t-test
one_sample_t_test <- function(sample_mean, population_mean, sample_std_dev, sample_size, alpha = 0.05) {
  # Calculate the t-statistic
  t_statistic = (sample_mean - population_mean) / (sample_std_dev / sqrt(sample_size))
  
  # Calculate the degrees of freedom
  df = sample_size - 1
  
  # Calculate the p-value
  p_value = 2 * pt(-abs(t_statistic), df, lower.tail = FALSE)
  
  # Determine the rejection region
  if (t_statistic > 0) {
    # One-tailed test with t > t-critical
    t_upper = round(qt(1 - alpha/2, df), 3)
    t_lower = "NONE"
  } else if (t_statistic < 0) {
    # One-tailed test with t < t-critical
    t_upper = "NONE"
    t_lower = round(-qt(1 - alpha/2, df), 3)
  } else {
    # Two-tailed test with t between -t-critical and t-critical
    t_upper = round(qt(1 - alpha/2, df), 3)
    t_lower = round(-qt(1 - alpha/2, df), 3)
  }
  
  # Print the results
  print(paste("t-statistic:", t_statistic))
  print(paste("p-value:", p_value))
  print(paste("Reject the null hypothesis if p-value <", alpha))
  print(paste("Rejection region: t <", t_lower, "or t >", t_upper))
  
  # Check if the p-value is less than alpha
  if (p_value < alpha) {
    print("Reject the null hypothesis. The sample mean is different from the population mean.")
  } else {
    print("Fail to reject the null hypothesis. The sample mean is not different from the population mean.")
  }
}

calculate_t_from_population_mean_differance <- function(population_mean_difference, standard_deviation, sample_size_n){
  result = population_mean_difference/(standard_deviation/sqrt(sample_size_n))
  return(result)
}

calculate_f_statistic <- function(sd1, n1, sd2, n2) {
  # Calculate the F-statistic
  f_statistic <- (sd1^2) / (sd2^2)
  
  # Return the F-statistic and degrees of freedom
  return(f_statistic)
}
calculate_p_value_from_f <- function(f_statistic, df1, df2, alternative = "two.sided") {
  # Calculate the p-value based on the alternative hypothesis
  if (alternative == "two.sided") {
    # Two-tailed test
    p_value <- 2 * min(pf(f_statistic, df1, df2), 1 - pf(f_statistic, df1, df2))
  } else if (alternative == "greater") {
    # One-tailed test (greater)
    p_value <- 1 - pf(f_statistic, df1, df2)
  } else if (alternative == "less") {
    # One-tailed test (less)
    p_value <- pf(f_statistic, df1, df2)
  } else {
    stop("Invalid alternative hypothesis. Choose 'two.sided', 'greater', or 'less'.")
  }
  
  # Return the p-value rounded to four decimal places
  return(p_value)
}

calculate_t_test_two_data_sets <- function(mean1, sd1, n1, mean2, sd2, n2, alternative = "two.sided") {
  # Calculate the pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Calculate the standard error
  se <- pooled_sd * sqrt(1/n1 + 1/n2)
  
  # Calculate the t-statistic
  t_statistic <- (mean1 - mean2) / se
  
  # Calculate degrees of freedom
  df <- n1 + n2 - 2
  
  # Calculate the p-value based on the alternative hypothesis
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t_statistic), df)
  } else if (alternative == "greater") {
    p_value <- pt(-t_statistic, df)
  } else if (alternative == "less") {
    p_value <- pt(t_statistic, df)
  } else {
    stop("Invalid alternative hypothesis. Choose 'two.sided', 'greater', or 'less'.")
  }
  
  # Return the t-statistic and p-value rounded to two and four decimal places, respectively
  return(t_statistic)
}

calculate_t_statistic_two_sample <- function(mean1, sd1, n1, mean2, sd2, n2) {
  # Calculate the pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * (sd1^2) + (n2 - 1) * (sd2^2)) / (n1 + n2 - 2))
  
  # Calculate the standard error
  se <- pooled_sd * sqrt((1/n1) + (1/n2))
  
  # Calculate the t-statistic
  t_statistic <- (mean1 - mean2) / se
  
  # Return the t-statistic rounded to two decimal places
  return(t_statistic)
}

determine_rejection_region <- function(test_statistic, df1, df2 = NULL, alpha = 0.05, test_type = "t", alternative = "two.sided") {
  # Determine the critical values based on the test type and alternative hypothesis
  if (test_type == "t") {
    if (alternative == "two.sided") {
      critical_value_low <- qt(alpha / 2, df1)
      critical_value_high <- qt(1 - alpha / 2, df1)
    } else if (alternative == "greater") {
      critical_value_low <- qt(1 - alpha, df1)
      critical_value_high <- "NONE"
    } else if (alternative == "less") {
      critical_value_low <- "NONE"
      critical_value_high <- qt(alpha, df1)
    } else {
      stop("Invalid alternative hypothesis. Choose 'two.sided', 'greater', or 'less'.")
    }
  } else if (test_type == "f") {
    if (alternative == "two.sided") {
      critical_value_low <- qf(alpha / 2, df1, df2)
      critical_value_high <- qf(1 - alpha / 2, df1, df2)
    } else if (alternative == "greater") {
      critical_value_low <- qf(1 - alpha, df1, df2)
      critical_value_high <- "NONE"
    } else if (alternative == "less") {
      critical_value_low <- "NONE"
      critical_value_high <- qf(alpha, df1, df2)
    } else {
      stop("Invalid alternative hypothesis. Choose 'two.sided', 'greater', or 'less'.")
    }
  } else {
    stop("Invalid test type. Choose 't' or 'f'.")
  }
  
  # Print the rejection region
  cat("Rejection region:\n")
  cat("test statistic <", ifelse(critical_value_low == "NONE", "NONE", round(critical_value_low, 6)), "\n")
  cat("test statistic >", ifelse(critical_value_high == "NONE", "NONE", round(critical_value_high, 6)), "\n")
}

calculate_p_value_t_test_for_two_data_sets <- function(t_statistic, n1, n2, alternative = "two.sided") {
  # Calculate degrees of freedom for equal variances
  df <- n1 + n2 - 2
  
  # Calculate the p-value based on the alternative hypothesis
  if (alternative == "two.sided") {
    p_value <- 2 * pt(-abs(t_statistic), df)
  } else if (alternative == "greater") {
    p_value <- pt(-t_statistic, df)
  } else if (alternative == "less") {
    p_value <- pt(t_statistic, df)
  } else {
    stop("Invalid alternative hypothesis. Choose 'two.sided', 'greater', or 'less'.")
  }
  
  # Return the p-value rounded to four decimal places
  return(p_value)
}

calculate_sample_variance <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sum_squares <- sum((x - mean_x)^2)
  Si <- (1 / (n - 1)) * sum_squares
  return(Si)
}
