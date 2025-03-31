source("./custom/myingarch.R")

search.ingarch <- function(x, 
                           max.p = 3,           # Maximum AR order
                           max.q = 3,           # Maximum MA order
                           distribution = c("poisson", "nbinom"),  # Distribution family
                           link = c("log", "identity"),        # Link function
                           xreg = NULL,         # External covariates
                           ic = c("aicc", "aic", "bic", "qic"),  # Information criterion
                           trace = FALSE,       # Whether to print model search progress
                           show_warnings = FALSE, # Whether to show warnings at the end
                           parallel = FALSE,    # Whether to use parallel processing
                           num.cores = 2,       # Number of cores for parallel processing
                           ...) {
  
  # Basic argument matching
  distribution <- match.arg(distribution)
  link <- match.arg(link)
  ic <- match.arg(ic)
  
  # Initialize model warnings collection
  model_warnings <- list()
  
  # Choose model orders
  # Serial - technically could be combined with the code below
  if (parallel == FALSE) {
    best.ic <- Inf
    model_count <- 0
    results <- list()
    
    for (i in 0:max.p) {
      for (j in 0:max.q) {
        model_count <- model_count + 1
        
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
        results[[model_count]] <- list(
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
    par.all.ingarch <- function(row) {
      params <- UndoWhichModels(row)
      i <- params["p"]
      j <- params["q"]
      
      # Create model key
      model_key <- sprintf("p%d_q%d", i, j)
      
      # Capture warnings
      warnings_captured <- character(0)
      fit <- tryCatch({
        withCallingHandlers({
          myingarch(x, order = c(i, j), ic = ic,
                    trace = trace, xreg = xreg,
                    distr = distribution, link = link, ...)
        }, warning = function(w) {
          warnings_captured <<- c(warnings_captured, w$message)
          invokeRestart("muffleWarning")
        })
      }, error = function(e) {
        # Return a model with Inf IC on error
        list(ic = Inf)
      })
      
      return(list(
        fit = fit,
        params = c(i, j),
        warnings = warnings_captured,
        model_key = model_key
      ))
    }
    
    if (is.null(num.cores)) {
      num.cores <- parallel::detectCores()
    }
    
    # Apply the function to all models in parallel
    all.models <- parallel::mclapply(
      X = split(to.check, 1:nrow(to.check)),
      FUN = par.all.ingarch,
      mc.cores = num.cores
    )
    
    # Process results
    best.ic <- Inf
    model_count <- 0
    results <- list()
    
    for (i in 1:length(all.models)) {
      model <- all.models[[i]]
      model_count <- model_count + 1
      
      # Store warnings
      model_warnings[[model$model_key]] <- model$warnings
      
      # Store results for display
      results[[model_count]] <- list(
        p = model$params[1],
        q = model$params[2],
        ic = model$fit$ic,
        model_key = model$model_key
      )
      
      if (!is.null(model$fit$ic) && model$fit$ic < best.ic) {
        bestfit <- model$fit
        best.ic <- bestfit$ic
        best_model_key <- model$model_key
      }
    }
  }
  
  if (!exists("bestfit") || bestfit$ic == Inf) {
    stop("No suitable INGARCH model found")
  }
  
  # Add useful information to the best model
  bestfit$model_warnings <- model_warnings
  bestfit$best_model_warnings <- model_warnings[[best_model_key]]
  bestfit$results <- results
  bestfit$ic <- NULL
  
  # Display warnings for all models if requested
  if (show_warnings) {
    cat("\nWarnings for INGARCH models tested:\n")
    cat("--------------------------------\n")
    
    # Convert results list to data frame
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    
    # Order by IC
    if (nrow(results_df) > 0) {
      results_df <- results_df[order(results_df$ic), ]
      
      # Get best model parameters
      best_p <- if(is.null(bestfit$model$past_obs)) 0 else length(bestfit$model$past_obs)
      best_q <- if(is.null(bestfit$model$past_mean)) 0 else length(bestfit$model$past_mean)
      
      # Display model warnings
      for (i in 1:nrow(results_df)) {
        p_val <- results_df$p[i]
        q_val <- results_df$q[i]
        ic_val <- results_df$ic[i]
        model_key <- as.character(results_df$model_key[i])
        
        # Mark best model with a star
        is_best <- (p_val == best_p && q_val == best_q)
        model_label <- if(is_best) "★ BEST MODEL" else ""
        
        cat(sprintf("INGARCH(%d,%d) with IC=%.4f %s\n", 
                    p_val, q_val, ic_val, model_label))
        
        # Display warnings for this model
        if (model_key %in% names(model_warnings) && length(model_warnings[[model_key]]) > 0) {
          for (w in model_warnings[[model_key]]) {
            cat("  - ", w, "\n", sep="")
          }
        } else {
          cat("  No warnings\n")
        }
        cat("\n")
      }
    }
  }
  
  if (trace) {
    cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
  }
  
  return(bestfit)
}