# Load required packages
library(tscount)

# Load our custom functions
# Ensure these file paths are correct relative to your working directory
# If they are in a subdirectory named "custom":
source("./ACTS/ingarch.sim.R") # This should be the path to your NEWLY MODIFIED ingarch.sim.R
source("./ACTS/myingarch.R")
source("./ACTS/auto.ingarch.R")


# Set a random seed for reproducibility
set.seed(123)

# Define the new parameters with vector coefficients
new_params <- list(
  intercept = 1.5, # Baseline level
  past_obs = c(0.05, 0.08), # Coefficients for Y(t-1), Y(t-2)
  past_mean = c(0.005, 0.0025, 0.002, 0.0035, 0.003, 0.0045) # Coefficients for mu(t-1)...mu(t-6)
  # If you were using external regressors, this list would contain an 'xreg' element with coefficients
  # e.g., xreg = c(0.1, -0.05)
)

# Define the corresponding model structure (lags)
# Assuming contiguous lags starting from 1
p_order <- length(new_params$past_obs) # Will be 2
q_order <- length(new_params$past_mean) # Will be 6

new_model_structure <- list(
  past_obs = if (p_order > 0) 1:p_order else NULL, # Will be c(1,2)
  past_mean = if (q_order > 0) 1:q_order else NULL, # Will be c(1,2,3,4,5,6)
  external = FALSE # No external regressors. If TRUE, you'd also need to provide 'xreg' data matrix to ingarch.sim
)

# Generate data with the new parameters and model structure
simulated_data <- ingarch.sim(
  n = 1000,                   # Generate 1000 observations
  param = new_params,
  model = new_model_structure,
  link = "identity",            # Use identity link function for simulation
  distr = "nbinom",             # Use Negative Binomial distribution for simulation
  distrcoefs = 2                # MODIFIED: Changed 'size = 2' to 'distrcoefs = 2'
  # If you were using external regressors, you would also pass the 'xreg = your_xreg_matrix' argument here.
)

# Convert to time series object for model fitting
ts_data <- ts(simulated_data$ts)

# Print the true parameters we used to generate the data
# Access them via $parameters_used and $distrcoefs_used as stored by your new ingarch.sim
cat("\nTrue parameters used in simulation:\n")
cat("Intercept:", simulated_data$parameters_used$intercept, "\n")

if (!is.null(simulated_data$parameters_used$past_obs) && length(simulated_data$parameters_used$past_obs) > 0) {
  cat("Past observations coefficients (p=", length(simulated_data$parameters_used$past_obs), ", lags=", paste(simulated_data$model_structure$past_obs, collapse=","), "): ", paste(simulated_data$parameters_used$past_obs, collapse=", "), "\n")
} else {
  cat("Past observations coefficients (p=0)\n")
}

if (!is.null(simulated_data$parameters_used$past_mean) && length(simulated_data$parameters_used$past_mean) > 0) {
  cat("Past means coefficients (q=", length(simulated_data$parameters_used$past_mean), ", lags=", paste(simulated_data$model_structure$past_mean, collapse=","), "): ", paste(simulated_data$parameters_used$past_mean, collapse=", "), "\n")
} else {
  cat("Past means coefficients (q=0)\n")
}

# MODIFIED: Accessing the NB dispersion parameter (now distrcoefs_used)
cat("True NB dispersion parameter (distrcoefs):", simulated_data$distrcoefs_used, "\n")
if (!is.null(simulated_data$distrcoefs_used)) {
  cat("Implied true alpha (1/distrcoefs) for NB dist (if var=mu+alpha*mu^2):", 1/simulated_data$distrcoefs_used, "\n")
}

# Now let's try to recover these parameters using auto.ingarch
# The max.p and max.q should be appropriate for the new true orders (p=2, q=6)
cat("\nFitting model with auto.ingarch...\n")
fitted_model <- try({
  auto.ingarch(
    y = ts_data,
    max.p = 7, # Maximum AR order (true is 2)
    max.q = 7, # Maximum MA order (true is 6)
    max.order = 14,
    ic = "aic",
    distribution = "nbinom",
    stepwise = TRUE,
    link = "log", # Note: Simulation used "identity", fitting uses "log". This is a valid choice.
    trace = TRUE,
    nmodels = 100, # You might want to increase this if max.p/max.q are larger
    parallel = FALSE,
    show_warnings = FALSE
  )
})

# --- ADDED SECTION TO EXTRACT MODEL INFORMATION INTO A DATAFRAME ---

# Check if model fitting was successful
if (inherits(fitted_model, "try-error")) {
  cat("\nError in model fitting:\n")
  print(fitted_model)
  model_info_df <- data.frame(Error = as.character(fitted_model))
} else if (is.null(fitted_model)) {
  cat("\nModel fitting returned NULL (no model selected or other issue).\n")
  model_info_df <- data.frame(Message = "No model selected or auto.ingarch returned NULL")
} else {
  cat("\nExtracting information from the fitted model...\n")
  
  # Assuming fitted_model is a 'tsglm' object or similar that summary() works on
  model_summary <- summary(fitted_model)
  
  # Initialize a list to store data frame rows
  df_rows <- list()
  
  # Extract coefficients
  if (!is.null(model_summary$coefficients) && nrow(model_summary$coefficients) > 0) {
    coef_table <- as.data.frame(model_summary$coefficients)
    expected_colnames <- c("Estimate", "Std. Error", "statistic", "p-value")
    if(ncol(coef_table) == length(expected_colnames)){
      colnames(coef_table) <- expected_colnames 
    } else {
      warning("Coefficient table has unexpected number of columns. Using existing names.")
    }
    
    param_names <- rownames(model_summary$coefficients)
    
    for (i in 1:nrow(coef_table)) {
      df_rows[[length(df_rows) + 1]] <- data.frame(
        Item = param_names[i],
        Estimate = coef_table[i, "Estimate"],
        Std.Error = coef_table[i, "Std. Error"],
        Statistic = coef_table[i, "statistic"],
        p.value = coef_table[i, "p-value"],
        stringsAsFactors = FALSE
      )
    }
  } else {
    cat("No coefficients found in model summary.\n")
  }
  
  # Extract dispersion parameter for Negative Binomial distribution
  if (tolower(fitted_model$distr) %in% c("nbinom", "nbinom2") && !is.null(model_summary$dispersion)) {
    disp_estimate <- NA
    disp_se <- NA
    
    if ("estimate" %in% names(model_summary$dispersion)) disp_estimate <- model_summary$dispersion["estimate"]
    if ("se" %in% names(model_summary$dispersion)) disp_se <- model_summary$dispersion["se"]
    
    df_rows[[length(df_rows) + 1]] <- data.frame(
      Item = "Overdispersion_alpha_for_NB (1/distrcoefs)", # Changed comment from size to distrcoefs
      Estimate = disp_estimate,
      Std.Error = disp_se,
      Statistic = NA, 
      p.value = NA,
      stringsAsFactors = FALSE
    )
  } else if (tolower(fitted_model$distr) == "nbinom1" && !is.null(model_summary$dispersion)) {
    disp_estimate <- NA
    disp_se <- NA
    if ("estimate" %in% names(model_summary$dispersion)) disp_estimate <- model_summary$dispersion["estimate"]
    if ("se" %in% names(model_summary$dispersion)) disp_se <- model_summary$dispersion["se"]
    df_rows[[length(df_rows) + 1]] <- data.frame(
      Item = "Dispersion_alpha_for_NB1",
      Estimate = disp_estimate,
      Std.Error = disp_se,
      Statistic = NA,
      p.value = NA,
      stringsAsFactors = FALSE
    )
  }
  
  # Extract other model statistics
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "LogLikelihood", Estimate = as.numeric(logLik(fitted_model)), Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "AIC", Estimate = AIC(fitted_model), Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "BIC", Estimate = BIC(fitted_model), Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  
  selected_p <- 0 
  if (!is.null(fitted_model$model$past_obs) && length(fitted_model$model$past_obs) > 0) {
    if(all(is.numeric(fitted_model$model$past_obs))) { 
      selected_p <- max(fitted_model$model$past_obs)
    } else if (is.logical(fitted_model$model$past_obs) && !fitted_model$model$past_obs) { 
      selected_p <- 0
    } else if(is.list(fitted_model$model$past_obs) && length(fitted_model$model$past_obs) > 0 && is.numeric(fitted_model$model$past_obs[[1]])) {
      selected_p <- max(unlist(fitted_model$model$past_obs))
    }
  }
  
  selected_q <- 0 
  if (!is.null(fitted_model$model$past_mean) && length(fitted_model$model$past_mean) > 0) {
    if(all(is.numeric(fitted_model$model$past_mean))) {
      selected_q <- max(fitted_model$model$past_mean)
    } else if (is.logical(fitted_model$model$past_mean) && !fitted_model$model$past_mean) {
      selected_q <- 0
    } else if(is.list(fitted_model$model$past_mean) && length(fitted_model$model$past_mean) > 0 && is.numeric(fitted_model$model$past_mean[[1]])) {
      selected_q <- max(unlist(fitted_model$model$past_mean))
    }
  }
  
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Selected_p_order_past_obs", Estimate = selected_p, Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Selected_q_order_past_mean", Estimate = selected_q, Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  
  if (isS3method("nobs", class(fitted_model))) { 
    df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Number_of_Observations", Estimate = nobs(fitted_model), Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  } else if (!is.null(fitted_model$nobs)) { 
    df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Number_of_Observations", Estimate = fitted_model$nobs, Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  } else {
    df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Number_of_Observations", Estimate = length(ts_data), Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE) 
  }
  
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Distribution", Estimate = fitted_model$distr, Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  df_rows[[length(df_rows) + 1]] <- data.frame(Item = "Link_Function", Estimate = fitted_model$link, Std.Error = NA, Statistic = NA, p.value = NA, stringsAsFactors = FALSE)
  
  if (length(df_rows) > 0) {
    model_info_df <- do.call(rbind, df_rows)
    rownames(model_info_df) <- NULL 
  } else {
    model_info_df <- data.frame(Message="No information extracted from the model.")
  }
  
  cat("\n--- Fitted Model Information Data Frame ---\n")
  print(model_info_df)
  
  tester <- tsglm(ts = ts_data, model = list(past_obs = c(1,2), past_mean = c(1,2,3,4,5,6), external = FALSE), distr = "nbinom", link = "log", xreg = NULL)
  print(summary(tester)$AIC)
}