source("./ACTS/myingarch.R") # Ensure this path is correct

search.ingarch <- function(x,
                           max.p = 3,                     # Maximum AR order
                           max.q = 3,                     # Maximum MA order
                           max.order = 5,                 #<< ADDED: Maximum sum of p and q
                           distribution = c("poisson", "nbinom"), # Distribution family
                           link = c("log", "identity"),   # Link function
                           xreg = NULL,                   # External covariates
                           ic = c("aicc", "aic", "bic", "qic"), # Information criterion
                           trace = FALSE,                 # Whether to print model search progress
                           show_warnings = FALSE,         # Whether to show warnings at the end
                           parallel = FALSE,              # Whether to use parallel processing
                           num.cores = 2,                 # Number of cores for parallel processing
                           ...) {
  
  # Basic argument matching
  distribution <- match.arg(distribution)
  link <- match.arg(link)
  ic <- match.arg(ic)
  
  # Initialize model warnings collection
  model_warnings <- list()
  
  # Choose model orders
  if (parallel == FALSE) {
    best.ic <- Inf
    model_count <- 0 # model_count for serial
    results <- list()
    
    for (i in 0:max.p) {
      for (j in 0:max.q) {
        if (i + j <= max.order) { # << ADDED: max.order condition for serial
          model_count <- model_count + 1 # Increment only if model is considered
          
          # Create a unique key for this model configuration
          model_key <- sprintf("p%d_q%d", i, j)
          
          # Capture warnings during model fitting
          warnings_captured <- character(0)
          fit <- withCallingHandlers({
            myingarch(x, order = c(i, j), ic = ic,
                      trace = trace, xreg = xreg,
                      distr = distribution, link = link, ...)
          }, warning = function(w) {
            warnings_captured <<- c(warnings_captured, w$message)
            invokeRestart("muffleWarning")
          })
          
          # Store warnings for this model
          model_warnings[[model_key]] <- warnings_captured
          
          # Store results for display
          results[[model_count]] <- list( # Use model_count as index, or model_key
            p = i, q = j,
            ic = fit$ic,
            model_key = model_key
          )
          
          if (fit$ic < best.ic) {
            best.ic <- fit$ic
            bestfit <- fit
            best_model_key <- model_key
          }
        }
      }
    }
  } else if (parallel == TRUE) {
    # Setup parallel processing functions similar to search.arima
    WhichModels <- function(max.p, max.q) {
      expand.grid(p = 0:max.p, q = 0:max.q)
    }
    
    UndoWhichModels <- function(row) {
      c(p = row$p, q = row$q)
    }
    
    to.check <- WhichModels(max.p, max.q)
    
    # Function to fit one model in parallel
    par.all.ingarch <- function(row) { # max.order will be in lexical scope
      params <- UndoWhichModels(row)
      i <- params["p"]
      j <- params["q"]
      
      if (i + j <= max.order) { # << ADDED: max.order condition for parallel
        # Create model key
        model_key <- sprintf("p%d_q%d", i, j)
        
        # Capture warnings
        warnings_captured <- character(0)
        fit <- tryCatch({
          withCallingHandlers({
            myingarch(x, order = c(i, j), ic = ic,
                      trace = trace, xreg = xreg, # Note: trace in parallel might be problematic for output interleaving
                      distr = distribution, link = link, ...)
          }, warning = function(w) {
            warnings_captured <<- c(warnings_captured, w$message)
            invokeRestart("muffleWarning")
          })
        }, error = function(e) {
          # Return a model with Inf IC on error to ensure it's not chosen
          list(ic = Inf, error_message = conditionMessage(e)) # Include error message if needed
        })
        
        return(list(
          fit = fit,
          params = c(i, j),
          warnings = warnings_captured,
          model_key = model_key
        ))
      } else {
        return(NULL) # << ADDED: Return NULL if model exceeds max.order
      }
    }
    
    if (is.null(num.cores)) {
      # num.cores <- parallel::detectCores() # Ensure parallel package is loaded or use as is if already loaded
      # For CRAN compatibility, it's often better to let user specify num.cores explicitly
      # or default to a known safe value like 2 if parallel::detectCores is not available/desired.
      # Your original code has parallel::detectCores(), so keeping it if `parallel` package is a dependency.
      if(requireNamespace("parallel", quietly = TRUE)) {
        num.cores <- parallel::detectCores()
      } else {
        warning("parallel package not available, setting num.cores to 1.")
        num.cores <- 1 # Fallback if parallel isn't available
      }
    }
    
    if (num.cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not found, running in serial. Please install 'parallel'.")
      # Re-run in serial or stop. For simplicity here, assuming it would run with num.cores=1 or error.
      # A full implementation might switch to the serial path.
      # Given the structure, mclapply will likely handle num.cores=1 fine if parallel isn't fully setup.
    }
    
    # Apply the function to all models in parallel
    # Using lapply for simplicity if not on Unix-like for mclapply, or use parLapply
    # Your original code uses mclapply, assuming Unix-like or appropriate setup
    if (.Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
      all.models <- parallel::mclapply(
        X = split(to.check, 1:nrow(to.check)), # split creates a list of rows
        FUN = par.all.ingarch,
        mc.cores = num.cores
      )
    } else {
      if (num.cores > 1) warning("mclapply is not reliably available on this OS or 'parallel' pkg missing; using lapply (serial within this block).")
      all.models <- lapply(
        X = split(to.check, 1:nrow(to.check)),
        FUN = par.all.ingarch
      )
    }
    
    
    all.models <- all.models[!sapply(all.models, is.null)] # << ADDED: Remove NULLs from skipped models
    
    # Process results
    best.ic <- Inf
    model_count <- 0 # model_count for parallel (number of actual fits processed)
    results <- list() # results list for parallel
    
    # Check if any models were actually fitted
    if (length(all.models) == 0) {
      if (!exists("bestfit")) { # Ensure bestfit is defined if no models run
        bestfit <- list(ic = Inf) # Default to an infinite IC if no models are viable
      }
    } else {
      for (i in 1:length(all.models)) {
        model_result <- all.models[[i]] # Renamed 'model' to 'model_result' to avoid conflict
        
        # Check if model_result itself or its fit component is NULL (should be filtered by now, but defensive)
        if (is.null(model_result) || is.null(model_result$fit)) next
        
        model_count <- model_count + 1 # Increment for models processed after filtering
        
        # Store warnings
        model_warnings[[model_result$model_key]] <- model_result$warnings
        
        # Store results for display
        results[[model_count]] <- list( # Use model_count or model_result$model_key
          p = model_result$params[1],
          q = model_result$params[2],
          ic = model_result$fit$ic, # Make sure fit$ic is present
          model_key = model_result$model_key
        )
        
        if (!is.null(model_result$fit$ic) && model_result$fit$ic < best.ic) {
          bestfit <- model_result$fit
          best.ic <- bestfit$ic
          best_model_key <- model_result$model_key
        }
      }
    }
  }
  
  if (!exists("bestfit") || bestfit$ic == Inf) {
    stop("No suitable INGARCH model found (all models resulted in IC=Inf or failed to fit).")
  }
  
  # Add useful information to the best model
  bestfit$model_warnings <- model_warnings
  if (exists("best_model_key") && best_model_key %in% names(model_warnings)) {
    bestfit$best_model_warnings <- model_warnings[[best_model_key]]
  } else {
    bestfit$best_model_warnings <- character(0) # No warnings if key doesn't exist
  }
  bestfit$results <- results # Store the list of all model ICs
  bestfit$ic_value_for_selection <- best.ic # Store the IC value that led to selection
  bestfit$ic <- NULL # Remove the specific IC field as per auto.arima style
  bestfit$n_models_evaluated <- model_count
  
  # Display warnings for all models if requested
  if (show_warnings) {
    cat("\nWarnings for INGARCH models tested:\n")
    cat("--------------------------------\n")
    
    # Convert results list to data frame for easier sorting and display
    if (length(results) > 0) {
      results_df <- do.call(rbind, lapply(results, function(r) {
        # Ensure all components are present and correctly typed for data.frame
        data.frame(
          p = as.integer(r$p), 
          q = as.integer(r$q), 
          ic = as.numeric(r$ic), 
          model_key = as.character(r$model_key),
          stringsAsFactors = FALSE
        )
      }))
      results_df <- results_df[!is.na(results_df$ic), ] # Remove models where IC might be NA
      
      if (nrow(results_df) > 0) {
        results_df <- results_df[order(results_df$ic), ]
        
        # Get best model parameters from the actual bestfit object
        # This assumes myingarch returns a model object with $model$past_obs etc.
        best_p_val <- NA
        best_q_val <- NA
        if(!is.null(bestfit$model)){
          best_p_val <- if(is.null(bestfit$model$past_obs) || length(bestfit$model$past_obs) == 0) 0 else max(bestfit$model$past_obs)
          best_q_val <- if(is.null(bestfit$model$past_mean) || length(bestfit$model$past_mean) == 0) 0 else max(bestfit$model$past_mean)
        } else if (exists("best_model_key")) { # Fallback if $model structure is not there
          # Try to parse from best_model_key if needed, though ideally bestfit$model is standard
          parsed_orders <- sscanf(best_model_key, "p%dq%d")
          if(length(parsed_orders)==2){
            best_p_val <- parsed_orders[1]
            best_q_val <- parsed_orders[2]
          }
        }
        
        
        for (i_row in 1:nrow(results_df)) {
          p_val <- results_df$p[i_row]
          q_val <- results_df$q[i_row]
          ic_val <- results_df$ic[i_row]
          model_key_val <- results_df$model_key[i_row] # Already character
          
          is_best <- (!is.na(best_p_val) && !is.na(best_q_val) && p_val == best_p_val && q_val == best_q_val)
          # Refined best check: also consider if this IC is the minimum found
          if(!is_best && !is.na(ic_val) && abs(ic_val - best.ic) < 1e-9) is_best <- TRUE
          
          
          model_label <- if(is_best) "★ BEST MODEL (or equivalent)" else ""
          
          cat(sprintf("INGARCH(%d,%d) with IC=%.4f %s\n", 
                      p_val, q_val, ic_val, model_label))
          
          if (model_key_val %in% names(model_warnings) && length(model_warnings[[model_key_val]]) > 0) {
            for (w_msg in model_warnings[[model_key_val]]) { # Renamed w to w_msg
              cat("  - ", w_msg, "\n", sep="")
            }
          } else {
            cat("  No warnings\n")
          }
          cat("\n")
        }
      } else {
        cat("No models were successfully fitted to display warnings for.\n")
      }
    } else {
      cat("No results to display warnings for.\n")
    }
  }
  
  if (trace) {
    # Assuming ingarch.string function exists and is compatible
    # if(exists("ingarch.string") && is.function(ingarch.string)){
    #   cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
    # } else {
    cat("\n\n Best model orders: p=", if(!is.null(bestfit$model$past_obs)) max(bestfit$model$past_obs) else 0,
        ", q=", if(!is.null(bestfit$model$past_mean)) max(bestfit$model$past_mean) else 0,
        " with ", ic, "=", best.ic, "\n\n", sep="")
    # }
  }
  
  return(bestfit)
}