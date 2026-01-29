## Functions for cmd-paper-heaping.R

smooth_age_heaping <- function(age_counts) {
  # Ensure the age_counts is named
  if (is.null(names(age_counts))) {
    stop("age_counts must be a named vector with ages as names.")
  }
  
  # Convert names to numeric for indexing
  ages <- as.numeric(names(age_counts))
  
  # Create a vector to store the smoothed counts
  smoothed_counts <- numeric(length(age_counts))
  
  # Apply the Carrier-Farrag formula
  for (i in 2:(length(age_counts) - 1)) {
    smoothed_counts[i] <- (age_counts[i-1] + 2 * age_counts[i] + age_counts[i+1]) / 4
  }
  
  # Handle the edge cases for the first and last ages
  smoothed_counts[1] <- (2 * age_counts[1] + age_counts[2]) / 3
  smoothed_counts[length(age_counts)] <- (age_counts[length(age_counts) - 1] + 2 * age_counts[length(age_counts)]) / 3
  
  # Return the smoothed counts with original names
  names(smoothed_counts) <- names(age_counts)
  
  return(smoothed_counts)
}



simulate_covariates <- function(age, strangeness = 1) {
  # Number of observations
  n_samp <- length(age)
  
  # Normalize age for scaling
  age_scaled <- scale(age)
  
  # Simulate income as a function of age with added noise
  income <- exp(3 + 0.5 * age_scaled + rnorm(n_samp, 0, strangeness))
  
  # Simulate marital status as a binary variable influenced by age
  marital_status <- ifelse(plogis(-3 + 0.05 * age_scaled + rnorm(n_samp, 0, strangeness)) > 0.5, 1, 0)
  
  # Simulate another covariate, e.g., education level (categorical)
  education_levels <- c("None", "Primary", "Secondary", "Tertiary")
  #education <- factor(sample(education_levels, n_samp, replace = TRUE, prob = plogis(-3 + 0.05 * age_scaled + rnorm(n_samp, 0, strangeness))))
  # Define a function to assign probabilities based on age
  get_education_prob <- function(age) {
    if (age < 20) {
      # Younger people are more likely to have no education or just primary/secondary education
      return(c(0.6, 0.3, 0.1, 0.0))  # None, Primary, Secondary, Tertiary
    } else if (age >= 20 & age < 35) {
      # For people in this age group, secondary and tertiary education are more likely
      return(c(0.1, 0.2, 0.5, 0.2))
    } else if (age >= 35 & age < 50) {
      # For mid-aged people, secondary and tertiary education are more common
      return(c(0.05, 0.1, 0.5, 0.35))
    } else {
      # Older people are more likely to have higher education (Secondary or Tertiary)
      return(c(0.05, 0.05, 0.3, 0.6))
    }
  }
  # Create an empty vector to store education levels
  education <- vector("character", length = n_samp)
  # Assign education levels to each person based on their age
  for (i in 1:n_samp) {
    myage <- age[i]
    prob <- get_education_prob(myage)
    education[i] <- sample(education_levels, size = 1, prob = prob)
  }
  
  
  # Combine into a data frame
  data <- data.frame(age = age, income = income, marital_status = marital_status, 
                     education = education)
  
  return(data)
}

fill_missing_ages <- function(age_counts, max_age = 95) {
  # Get the minimum age from the age counts
  min_age <- min(as.numeric(names(age_counts)))
  
  # Create a sequence of ages from the minimum age to 95
  all_ages <- seq(min_age, max_age)
  
  # Create a named vector with all ages, initialized to zero
  complete_age_counts <- setNames(rep(0, length(all_ages)), all_ages)
  
  # Fill in the counts from the original data
  complete_age_counts[names(age_counts)] <- age_counts
  
  return(complete_age_counts)
}

# Now introduce_heaping rounds to age5 randomly. The model can thus no better.
# When rounding almost always up, then the model might correct this?
# But: first is rounded, then the values are changed, and afterwards
# it is decided if the change was in a good direction. 


introduce_heaping <- function(df, age_var, heap_ratio, interval = 5) {
  # Convert to data.table if not already
  is_data_table <- is.data.table(df)
  # if (!is_data_table) {
  #   setDT(df)
  # }
  df <- data.frame(df)
  
  # Check if the age variable exists in the data frame
  if (!age_var %in% names(df)) {
    stop("The specified age variable does not exist in the data frame.")
  }
  
  # Ensure the age variable is numeric
  if (!is.numeric(df[[age_var]])) {
    stop("The specified age variable must be numeric.")
  }
  
  # Determine the number of values to change
  n <- nrow(df)
  n_to_heap <- round(heap_ratio * n)
  
  # Sample indices to heap
  indices_to_heap <- sample(n, n_to_heap)
  
  # Round ages to the nearest specified interval
  #df[indices_to_heap, (age_var) := round(df[[age_var]][indices_to_heap] / interval) * interval]
  df[indices_to_heap, age_var] <- round(df[indices_to_heap, age_var] / interval) * interval
  
  # Convert back to data.frame if necessary
  # if (!is_data_table) {
  #   setDF(df)
  # }
  if (!is_data_table) {
    df <- data.table(df)
  }
  return(df)
}

# Define the function
handle_missing_values <- function(df) {
  df <- df %>%
    mutate(across(where(is.factor), ~ ifelse(is.na(.), "not applicable", as.character(.))),
           across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))
  
  # Convert modified factor columns back to factor type
  df <- df %>%
    mutate(across(where(is.character), ~ as.factor(.)))
  
  return(df)
}

