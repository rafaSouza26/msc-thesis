# Ensure these paths are correct for your environment
source("./ACTS/newmodel.R")       # Load helper function for checking if a model is new
source("./ACTS/ingarch.string.R") # Load helper function for formatting model strings
source("./ACTS/search.ingarch.R") # Load function for non-stepwise model search

# Function to automatically find the best INGARCH(p,q) model for a time series.
auto.ingarch <- function(y,                          # Input time series (count data)
                         max.p = 7,                  # Maximum AR (past observations) order
                         max.q = 7,                  # Maximum MA (past conditional means) order
                         max.order = 5,              # Maximum sum of p+q for non-stepwise search (search.ingarch)
                         start.p = 2,                # Initial p for stepwise search
                         start.q = 2,                # Initial q for stepwise search
                         distribution = c("poisson", "nbinom"), # Distribution for the counts
                         link = c("log", "identity"),    # Link function for the conditional mean
                         xreg = NULL,                # Optional matrix of external regressors
                         ic = c("aicc", "aic", "bic", "qic"), # Information criterion for model selection
                         stepwise = TRUE,            # If TRUE, use stepwise search; otherwise, search all models
                         nmodels = 94,               # Maximum number of models to evaluate in stepwise search
                         trace = FALSE,              # If TRUE, print progress during model search
                         show_warnings = FALSE,      # If TRUE, display all warnings from model fits at the end
                         method = NULL,              # Estimation method (passed to myingarch if supported)
                         parallel = FALSE,           # If TRUE, use parallel processing (only for non-stepwise)
                         num.cores = 2,              # Number of cores for parallel processing
                         x = y,                      # Allow 'x' as an alias for 'y' for internal use
                         ...) {                      # Other arguments passed to myingarch
  
  # --- 1. Initial Validations & Setup ---
  
  # Validate parallel processing settings
  if (stepwise && parallel) {
    warning("Parallel computing is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE # Disable parallel if stepwise is TRUE
  }
  
  # Warning for trace output in parallel non-stepwise mode
  if (trace && parallel) { # parallel can only be TRUE here if stepwise is FALSE
    message("Tracing model searching in parallel (stepwise=FALSE) might produce interleaved output.")
  }
  
  # Store series name and ensure input is a time series object
  series <- deparse(substitute(y)) # Get the name of the input series
  # x_original_arg <- x # Original argument x (or y) before ts conversion
  x <- as.ts(x)    # Convert input to a time series object; 'x' is the working series
  orig.x <- x     # Store the initial ts object (before NA trim) for the output 'fit$x'
  
  # Validate max.order argument
  if (!is.numeric(max.order) || max.order < 0 || length(max.order) != 1) {
    warning("max.order must be a single non-negative integer. Setting to default (5).")
    max.order <- 5
  }
  
  # Validate start.p and start.q against max.p and max.q
  if (!is.null(start.p) && start.p > max.p) {
    warning("start.p cannot be greater than max.p, setting start.p = max.p")
    start.p <- max.p
  }
  if (!is.null(start.q) && start.q > max.q) {
    warning("start.q cannot be greater than max.q, setting start.q = max.q")
    start.q <- max.q
  }
  
  # Check for univariate series
  if (NCOL(x) > 1) {
    stop("auto.ingarch can only handle univariate time series")
  }
  
  # Check for non-negative integer (count) data
  if (!all(x[!is.na(x)] >= 0) || !all(x[!is.na(x)] == floor(x[!is.na(x)]))) {
    stop("INGARCH models require non-negative integer (count) data")
  }
  
  # --- 2. Handle Leading NAs & Determine Effective Series Length ---
  
  firstnonmiss_idx <- which(!is.na(x))[1] # Find the index of the first non-missing value
  
  # Stop if all data are NAs or the series is empty
  if (is.na(firstnonmiss_idx)) {
    stop("Time series 'y' contains only NAs or is empty.")
  }
  # Calculate effective series length for p,q capping (count of non-NAs from first non-NA point)
  current_serieslength <- sum(!is.na(x[firstnonmiss_idx:length(x)]))
  
  # Trim leading NAs from 'x' and adjust 'xreg' accordingly
  if (firstnonmiss_idx > 1) {
    x <- window(x, start = time(x)[firstnonmiss_idx]) # Trim 'x'
    if (!is.null(xreg)) {
      xreg_orig_nrows <- NROW(xreg)
      xreg <- as.matrix(xreg) # Ensure xreg is a matrix before subsetting
      if (xreg_orig_nrows >= firstnonmiss_idx) {
        xreg <- xreg[firstnonmiss_idx:xreg_orig_nrows, , drop = FALSE] # Trim 'xreg'
      } else {
        # This case implies xreg is shorter than the original series up to the first non-NA.
        # The subsequent length check (NROW(xreg) != length(x)) should catch mismatches.
        stop("xreg is too short to be trimmed to align with leading NAs in y.")
      }
    }
  }
  
  # Match and validate key arguments (ic, distribution, link)
  ic <- match.arg(ic) # 'ic' variable now holds the validated IC name (e.g., "aic")
  distribution_matched <- match.arg(distribution)
  link_matched <- match.arg(link)
  
  # --- 3. External Regressors (xreg) Validation (post NA-trimming) ---
  if (!is.null(xreg)) {
    if (!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg) # Ensure xreg is a matrix
    
    # Check if number of rows in xreg matches the (potentially trimmed) length of x
    if (NROW(xreg) != length(x))
      stop(paste("Number of rows in xreg (", NROW(xreg), ") must match length of time series (", length(x), ") after NA handling.", sep = ""))
    
    # Check for NAs in xreg (after trimming)
    if (any(is.na(xreg)))
      stop("Missing values in external regressors are not allowed after initial NA trim of y.")
    
    # Assign default column names if not present
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("xreg", 1:NCOL(xreg))
    }
    
    # Validate xreg values if using identity link
    if (link_matched == "identity" && any(xreg < 0)) {
      stop("When using 'identity' link, all values in xreg must be non-negative")
    }
  }
  
  # --- 4. Handle Constant Data Series ---
  # Check if all non-NA values in the (trimmed) series 'x' are the same
  first_val_for_const_check <- x[which(!is.na(x))[1]]
  if (all(x == first_val_for_const_check, na.rm = TRUE)) {
    if(all(is.na(x))) # Should be caught earlier if x became all NAs
      stop("All data are missing")
    
    const_value <- first_val_for_const_check
    # Create a simplified fit object for constant series
    fit <- list(
      ts = orig.x, # The original, untrimmed time series object
      model = list(past_obs = 0, past_mean = 0, external = FALSE), # external is FALSE as per original older script
      distribution = distribution_matched,
      link = link_matched,
      coefficients = if (is.null(xreg)) setNames(const_value, "(Intercept)") else setNames(const_value, "(Intercept)"), # Simplified coefficients
      final_estimates = const_value, # From older script
      fitted.values = rep(const_value, length(orig.x)), # Fitted values over original length
      call = match.call(), # Store the function call
      x = orig.x,          # Store the original data (as a ts object)
      series = series      # Store the original series name
    )
    if (fit$distribution == "nbinom") { # If nbinom, dispersion parameter (size) is effectively Inf for Poisson-like constant
      fit$size <- Inf
    }
    class(fit) <- "tsglm" # Assign class for compatibility
    return(fit)
  }
  
  # --- 5. Adjust Max Orders & IC for Short Series ---
  # Cap max.p and max.q based on the effective series length
  max.p <- min(max.p, floor(current_serieslength / 3))
  max.q <- min(max.q, floor(current_serieslength / 3))
  
  # For very short series, force AIC as other criteria might be unstable
  if (current_serieslength <= 3L) {
    ic <- "aic" # Note: 'ic' variable is updated here if condition met
  }
  
  # --- 6. Model Search Strategy ---
  
  # --- 6a. Non-Stepwise (Brute-Force) Search ---
  if (!stepwise) {
    if (trace) cat("Searching for best INGARCH model using non-stepwise (brute-force) method...\n")
    bestfit <- search.ingarch(
      x = x,                           # Pass the (NA-trimmed) series
      max.p = max.p,
      max.q = max.q,
      max.order = max.order,           # Max sum of p+q for search.ingarch
      distribution = distribution_matched,
      link = link_matched,
      xreg = xreg,
      ic = ic,                         # Pass the chosen information criterion
      trace = trace,                   # Pass trace argument to search.ingarch
      show_warnings = show_warnings,
      parallel = parallel,
      num.cores = num.cores,
      ...                              # Pass other arguments
    )
    
    # Finalize the bestfit object from search.ingarch
    bestfit$x <- orig.x        # Attach original, untrimmed ts data
    bestfit$series <- series
    bestfit$ic <- NULL         # Remove IC field as per older script; IC value itself is in model output
    bestfit$call <- match.call()
    
    # Process results from search.ingarch to ensure correct matrix format
    # (This logic adopted from the "newer" script for robustness)
    if (!is.null(bestfit$results) && is.list(bestfit$results) && !is.matrix(bestfit$results)) {
      temp_results_df <- tryCatch(do.call(rbind, lapply(bestfit$results, function(item) {
        if(is.list(item) && all(c("p", "q", "ic") %in% names(item))) {
          data.frame(p = item$p, q = item$q, ic = item$ic)
        } else { NULL } # If item is not structured as expected, return NULL
      })), error = function(e) NULL) # Catch errors during rbind/lapply
      if(!is.null(temp_results_df) && nrow(temp_results_df) > 0) {
        bestfit$results <- as.matrix(temp_results_df)
      } else { # If no valid results, create an empty matrix
        bestfit$results <- matrix(numeric(0), ncol = 3, dimnames = list(NULL, c("p", "q", "ic")))
      }
    } else if (is.null(bestfit$results) || !is.matrix(bestfit$results)){ # If results are NULL or not already a matrix
      bestfit$results <- matrix(numeric(0), ncol = 3, dimnames = list(NULL, c("p", "q", "ic")))
    }
    
    # Trace output for the best model found by search.ingarch
    if (trace && !is.null(bestfit$model)) { # Check if a model was actually found
      cat("\n\n Best model from search.ingarch:", ingarch.string(bestfit, padding = TRUE), "\n\n")
    }
    
    # Adjust very small coefficients to prevent issues (safer version)
    if (!is.null(bestfit$coefficients)) {
      mask <- abs(bestfit$coefficients) < 1e-4 & bestfit$coefficients != 0 # Avoid affecting actual zeros
      bestfit$coefficients[mask] <- 1e-4 * sign(bestfit$coefficients[mask])
    }
    return(bestfit) # Return the best model from non-stepwise search
  }
  
  # --- 6b. Stepwise Search ---
  # Adjust starting p, q for very short series if using stepwise
  if (length(x) < 10L) { # length(x) is the (NA-trimmed) series
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  # Initialize current best p and q
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  
  # Initialize structures for storing warnings and results
  model_warnings <- list() # List to store warnings for each (p,q) combination
  results <- matrix(NA, nrow = nmodels, ncol = 3) # Matrix to store p, q, and IC for each model
  colnames(results) <- c("p", "q", "ic")
  model_count <- 0 # Counter for total models evaluated
  
  # --- Fit Initial Models for Stepwise Search ---
  if (trace) cat("Starting stepwise model search...\n")
  
  # 1. Initial model (start.p, start.q)
  warnings_captured <- character(0) # Temp store for warnings of current model
  bestfit <- withCallingHandlers({    # Fit model and capture warnings
    model_count <- model_count + 1 # Increment model counter
    myingarch(x, order = c(p, q), ic = ic, trace = trace, xreg = xreg, # Call internal fitting function
              distr = distribution_matched, link = link_matched, ...)
  }, warning = function(w) { # Warning handler
    warnings_captured <<- c(warnings_captured, conditionMessage(w)) # Store warning message
    invokeRestart("muffleWarning") # Suppress a R's default warning display
  })
  results[model_count, ] <- c(p, q, bestfit$ic) # Store results
  model_warnings[[paste0("p", p, "_q", q)]] <- warnings_captured # Store warnings
  k <- model_count # 'k' tracks the number of rows filled in the 'results' matrix
  
  # 2. INGARCH(0,0) model (intercept only, or with xreg)
  if (k < nmodels && (p != 0 || q != 0)) { # Check if (0,0) is different from initial
    warnings_captured <- character(0)
    fit <- withCallingHandlers({
      model_count <- model_count + 1
      myingarch(x, order = c(0, 0), ic = ic, trace = trace, xreg = xreg,
                distr = distribution_matched, link = link_matched, ...)
    }, warning = function(w) {
      warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    k <- k + 1
    results[k, ] <- c(0, 0, fit$ic)
    model_warnings[[paste0("p0_q0")]] <- warnings_captured
    if (fit$ic < bestfit$ic) { bestfit <- fit; p <- 0; q <- 0; } # Update if better
  }
  
  # 3. INGARCH(1,0) model
  if (max.p > 0 && (p != 1 || q != 0)) { # If p can be 1, and it's not the current best
    if (k < nmodels) {
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(1, 0), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      k <- k + 1
      results[k, ] <- c(1, 0, fit$ic)
      model_warnings[[paste0("p1_q0")]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p <- 1; q <- 0; }
    }
  }
  
  # 4. INGARCH(0,1) model
  if (max.q > 0 && (p != 0 || q != 1)) { # If q can be 1, and it's not the current best
    if (k < nmodels) {
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(0, 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      k <- k + 1
      results[k, ] <- c(0, 1, fit$ic)
      model_warnings[[paste0("p0_q1")]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p <- 0; q <- 1; }
    }
  }
  
  # 5. INGARCH(1,1) model (explicitly requested addition)
  if (max.p >= 1 && max.q >= 1 && (p != 1 || q != 1)) { # If p,q can be 1, and it's not current best
    if (k < nmodels) {
      warnings_captured <- character(0)
      fit_11 <- withCallingHandlers({ # Use a distinct variable name for clarity
        model_count <- model_count + 1
        myingarch(x, order = c(1, 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      k <- k + 1
      results[k, ] <- c(1, 1, fit_11$ic)
      model_warnings[[paste0("p1_q1")]] <- warnings_captured
      if (fit_11$ic < bestfit$ic) { bestfit <- fit_11; p <- 1; q <- 1; }
    }
  }
  
  startk <- 0 # Control variable for the main stepwise loop iteration
  if (trace && k > 0) { # Context message before starting the main search loop
    cat("\nFitting through step-wise search now (exploring neighbors of current best)...\n")
  }
  
  # --- Main Stepwise Search Loop (Adhering to older script's explicit structure) ---
  # This loop iteratively explores neighbors of the current best (p,q) model.
  # 'p' and 'q' store the orders of the current best model found so far.
  # 'k' tracks the number of models evaluated and stored in 'results'.
  # 'startk' is used to check if any new, better model was found in an iteration.
  # 'nmodels' is the overall limit on models to fit.
  while (startk < k && k < nmodels) {
    startk <- k # Record 'k' at the start of this iteration
    
    # Candidate 1: Try (p-1, q)
    if (p > 0 && newmodel(p - 1, q, results[1:k, , drop=FALSE])) { # Check if valid and new
      k <- k + 1 # Increment count of models in results table
      if(k > nmodels) break # Stop if model limit reached
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1 # Increment total models evaluated
        myingarch(x, order = c(p - 1, q), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p - 1, q, fit$ic) # Store result
      model_warnings[[paste0("p", p - 1, "_q", q)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p <- p - 1; next } # If better, update best and restart loop from new best
    }
    
    # Candidate 2: Try (p, q+1)
    if (q < max.q && newmodel(p, q + 1, results[1:k, , drop=FALSE])) {
      k <- k + 1
      if(k > nmodels) break
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(p, q + 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p, q + 1, fit$ic)
      model_warnings[[paste0("p", p, "_q", q + 1)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; q <- q + 1; next }
    }
    
    # Candidate 3: Try (p-1, q-1)
    if (q > 0 && p > 0 && newmodel(p - 1, q - 1, results[1:k, , drop=FALSE])) {
      k <- k + 1
      if(k > nmodels) break
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(p - 1, q - 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p - 1, q - 1, fit$ic)
      model_warnings[[paste0("p", p - 1, "_q", q - 1)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; q <- q - 1; p <- p - 1; next }
    }
    
    # Candidate 4: Try (p-1, q+1)
    if (q < max.q && p > 0 && newmodel(p - 1, q + 1, results[1:k, , drop=FALSE])) {
      k <- k + 1
      if(k > nmodels) break
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(p - 1, q + 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p - 1, q + 1, fit$ic)
      model_warnings[[paste0("p", p - 1, "_q", q + 1)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; q <- q + 1; p <- p - 1; next }
    }
    
    # Candidate 5: Try (p+1, q-1)
    if (q > 0 && p < max.p && newmodel(p + 1, q - 1, results[1:k, , drop=FALSE])) {
      k <- k + 1
      if(k > nmodels) break
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(p + 1, q - 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p + 1, q - 1, fit$ic)
      model_warnings[[paste0("p", p + 1, "_q", q - 1)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; q <- q - 1; p <- p + 1; next }
    }
    
    # Candidate 6: Try (p+1, q+1)
    if (q < max.q && p < max.p && newmodel(p + 1, q + 1, results[1:k, , drop=FALSE])) {
      k <- k + 1
      if(k > nmodels) break
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        myingarch(x, order = c(p + 1, q + 1), ic = ic, trace = trace, xreg = xreg,
                  distr = distribution_matched, link = link_matched, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      results[k, ] <- c(p + 1, q + 1, fit$ic)
      model_warnings[[paste0("p", p + 1, "_q", q + 1)]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; q <- q + 1; p <- p + 1; next }
    }
  } # End of while loop for stepwise search
  
  # Warning if model limit was reached during an iteration where models were still being added
  if (k >= nmodels && startk < k) { # startk < k indicates loop made progress before hitting nmodels
    warning(sprintf("Stepwise search was stopped after %d models due to reaching the model number limit: `nmodels = %i`", model_count, nmodels))
  }
  
  # --- 7. Finalize and Return Best Model (from Stepwise Search) ---
  
  # Check if any suitable model was found (i.e., IC is not Inf)
  if (is.infinite(bestfit$ic)) { # Note: bestfit$ic holds the IC value of the current best model
    stop("No suitable INGARCH model found (all resulted in Inf IC or errors).")
  }
  
  # Store the actual IC value of the best model before nullifying bestfit$ic field
  # This is used for trace and show_warnings if those rely on this value directly.
  actual_best_ic_value <- bestfit$ic
  
  # Populate final fields in the bestfit object
  bestfit$x <- orig.x        # Original, untrimmed ts data
  bestfit$series <- series    # Name of the original series
  bestfit$ic <- NULL          # Remove 'ic' field as per original script's practice (value was used for selection)
  bestfit$call <- match.call()# Store the function call
  bestfit$results <- results[1:k, , drop = FALSE] # Table of (p,q,IC) for all models considered
  bestfit$model_warnings <- model_warnings      # List of warnings for all models
  
  # Determine p and q orders of the final best model for fetching its specific warnings
  best_p_val <- NA_integer_ ; best_q_val <- NA_integer_
  if (!is.null(bestfit$model)) { # Check if model component exists
    past_obs_component <- bestfit$model$past_obs
    past_mean_component <- bestfit$model$past_mean
    # Determine orders based on length of coefficient vectors for past_obs/past_mean
    best_p_val <- if (is.null(past_obs_component) || !is.numeric(past_obs_component)) 0 else length(past_obs_component)
    best_q_val <- if (is.null(past_mean_component) || !is.numeric(past_mean_component)) 0 else length(past_mean_component)
    
    if (is.numeric(best_p_val) && is.numeric(best_q_val)) { # Ensure p,q are numeric
      bestfit$best_model_warnings <- model_warnings[[paste0("p", best_p_val, "_q", best_q_val)]] # Warnings for the best model
    } else {
      bestfit$best_model_warnings <- character(0) # Fallback
    }
  } else {
    bestfit$best_model_warnings <- character(0) # Fallback if model structure is unexpected
  }
  bestfit$n_total_models <- model_count # Total number of models evaluated
  
  # Trace output for the final best model from stepwise search
  if (trace) {
    cat("\n") # Newline for clarity
    if (!is.null(bestfit$model) && exists("ingarch.string", mode = "function")) {
      cat("\nBest model:", ingarch.string(bestfit, padding = TRUE), "\n")
      # Use 'ic' (function arg for IC name) and 'actual_best_ic_value' for the trace display
    } else if (!is.null(bestfit$model) && !is.na(best_p_val) && !is.na(best_q_val) &&
               !is.null(actual_best_ic_value) && !is.infinite(actual_best_ic_value) && !is.null(ic)) {
      cat("\nBest model: INGARCH(", best_p_val, ",", best_q_val, ") with ", ic, " = ", sprintf("%.4f", actual_best_ic_value), "\n", sep = "")
    }
    cat("Total models evaluated:", model_count, "\n")
  }
  
  # --- 8. Display All Warnings (if requested) ---
  # This section uses the more detailed warning display logic from the "newer" script.
  if (show_warnings && k > 0) { # 'k' is the number of models in the results table
    cat("\nWarnings for all models tested:\n")
    cat("------------------------------\n")
    valid_results_df <- as.data.frame(results[1:k, , drop = FALSE]) # Convert results to data frame
    valid_results_df <- valid_results_df[!is.na(valid_results_df$ic), , drop = FALSE] # Remove models with NA IC
    
    if (nrow(valid_results_df) > 0) {
      # Order models by their IC value
      ordered_results <- valid_results_df[order(valid_results_df$ic), , drop = FALSE]
      
      # 'actual_best_ic_value', 'best_p_val', 'best_q_val' are already determined for the best model
      # 'ic' (function argument) holds the name of the criterion used (e.g., "aic")
      
      for (i_row in 1:nrow(ordered_results)) {
        p_val_loop <- ordered_results[i_row, "p"]
        q_val_loop <- ordered_results[i_row, "q"]
        ic_val_loop <- ordered_results[i_row, "ic"]
        
        # Check if this model from the results table is the overall best model
        is_best <- (!is.na(best_p_val) && !is.na(best_q_val) && # Ensure best p,q are valid
                      p_val_loop == best_p_val && q_val_loop == best_q_val &&
                      !is.na(actual_best_ic_value) && !is.na(ic_val_loop) && # Ensure ICs are valid
                      abs(ic_val_loop - actual_best_ic_value) < 1e-9) # Compare IC values
        model_label <- if (is_best) "--- ★ BEST MODEL" else ""
        
        # Print model info and its IC; use 'ic' (function arg) for the IC name
        cat(sprintf("INGARCH(%d,%d) with %s = %.4f %s\n",
                    p_val_loop, q_val_loop, ic, ic_val_loop, model_label))
        
        # Print warnings for this specific model
        warnings_key <- paste0("p", p_val_loop, "_q", q_val_loop)
        if (warnings_key %in% names(model_warnings) && length(model_warnings[[warnings_key]]) > 0) {
          for (w_msg in model_warnings[[warnings_key]]) {
            cat("  - ", w_msg, "\n", sep = "")
          }
        } else {
          cat("  No warnings\n")
        }
        cat("\n") # Blank line for readability
      }
    } else {
      cat("No models with valid ICs to display warnings for.\n")
    }
  }
  
  # Final adjustment for very small coefficients (safer version)
  if (!is.null(bestfit$coefficients)) {
    mask <- abs(bestfit$coefficients) < 1e-4 & bestfit$coefficients != 0 # Avoid affecting actual zeros
    bestfit$coefficients[mask] <- 1e-4 * sign(bestfit$coefficients[mask])
  }
  
  return(bestfit) # Return the final best model object
}