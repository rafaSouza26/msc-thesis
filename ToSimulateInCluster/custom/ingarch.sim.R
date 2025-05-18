ingarch.sim <- function(n,
                        param = list(intercept = 1, past_obs = NULL, past_mean = NULL, external = NULL),
                        model = list(past_obs = NULL, past_mean = NULL, xreg = FALSE),
                        xreg = NULL, # Covariate matrix
                        link = c("identity", "log"),
                        distr = c("poisson", "nbinom"),
                        size = NULL,  # Use 'size' directly for nbinom
                        fit = NULL,   # Option to simulate from a fitted model
                        n_start = 50) {
  
  # --- Input Validation ---
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(n_start) || length(n_start) != 1 || n_start < 0 || n_start != floor(n_start)) {
    stop("n_start must be a non-negative integer")
  }
  
  # --- Parameters from fit object ---
  if (!is.null(fit)) {
    if (!inherits(fit, "tsglm")) {
      stop("fit must be an object of class 'tsglm'")
    }
    tryCatch({
      model <- fit$model
      link <- fit$link
      distr <- fit$distr
      # Extract 'size' if nbinom
      size <- if(distr == "nbinom") fit$distrcoefs else NULL
      
      # Determine coefficient indices carefully
      coefs <- fit$coefficients
      names_coefs <- names(coefs)
      intercept_ind <- which(names_coefs == "intercept")
      past_obs_inds <- which(grepl("^beta", names_coefs))
      past_mean_inds <- which(grepl("^alpha", names_coefs))
      external_inds <- which(! 1:length(coefs) %in% c(intercept_ind, past_obs_inds, past_mean_inds))
      
      param <- list(
        intercept = coefs[intercept_ind],
        past_obs = if(length(past_obs_inds) > 0) coefs[past_obs_inds] else NULL,
        past_mean = if(length(past_mean_inds) > 0) coefs[past_mean_inds] else NULL,
        external = if(length(external_inds) > 0) coefs[external_inds] else NULL
      )
      # If simulating from fit, xreg should ideally be provided matching fit$xreg
      if (is.null(xreg) && length(external_inds) > 0) {
        warning("Simulating from 'fit' object with external regressors, but 'xreg' was not provided.")
        # Attempt to use xreg from fit object if available, otherwise error
        if (!is.null(fit$xreg)) {
          if (nrow(fit$xreg) >= n) {
            xreg <- fit$xreg[1:n, , drop = FALSE] # Use first n rows from fit$xreg
            #cat("Using first", n, "rows of xreg from the provided 'fit' object.\n")
          } else {
            stop("Provided 'fit' object has fewer xreg rows than requested simulation length 'n'.")
          }
        } else {
          stop("Cannot simulate external regressors effects from 'fit' object without providing 'xreg'.")
        }
      }
      
    }, error = function(e) {
      stop("Failed to extract parameters from fit object: ", conditionMessage(e))
    })
  }
  
  # --- Validate link and distribution ---
  link <- match.arg(link)
  distr <- match.arg(distr)
  
  # --- Validate parameters (basic checks) ---
  if (!is.numeric(param$intercept) || length(param$intercept) != 1) {
    stop("param$intercept must be a single numeric value")
  }
  if (!is.null(param$past_obs) && !is.numeric(param$past_obs)) {
    stop("param$past_obs coefficients must be numeric")
  }
  if (!is.null(param$past_mean) && !is.numeric(param$past_mean)) {
    stop("param$past_mean coefficients must be numeric")
  }
  if (!is.null(param$external) && !is.numeric(param$external)) {
    stop("param$external coefficients must be numeric")
  }
  
  # --- Validate distribution parameters ---
  if (distr == "nbinom") {
    if (is.null(size)) {
      stop("size parameter must be provided for negative binomial distribution")
    }
    if (!is.numeric(size) || length(size) != 1 || size <= 0) {
      stop("size parameter must be a single positive numeric value for nbinom distribution")
    }
  } else {
    size <- NULL # Ensure size is NULL if not nbinom
  }
  
  # --- Validate external regressors xreg ---
  n_external_coefs <- length(param$external)
  has_external_model <- model$external || (n_external_coefs > 0) # Check if external effects expected
  
  if (has_external_model) {
    if (is.null(xreg)) {
      stop("External regressors defined in model/parameters, but 'xreg' matrix was not provided.")
    }
    if (!is.numeric(xreg)) {
      stop("xreg must be a numeric matrix or vector")
    }
    xreg <- as.matrix(xreg) # Ensure it's a matrix
    if (nrow(xreg) != n) {
      stop(sprintf("Number of rows in xreg (%d) must match n (%d)", nrow(xreg), n))
    }
    if (ncol(xreg) != n_external_coefs) {
      stop(sprintf("Number of columns in xreg (%d) must match number of external coefficients (%d)", ncol(xreg), n_external_coefs))
    }
    if(anyNA(xreg)){
      stop("Missing values (NA) are not allowed in xreg.")
    }
  } else {
    # If no external effects expected, ensure xreg is NULL
    if (!is.null(xreg)) {
      warning("xreg provided but no external coefficients found in parameters/model. xreg will be ignored.")
      xreg <- NULL
    }
    if (!is.null(param$external)) {
      warning("External coefficients provided but model$external=FALSE. Coefficients will be ignored.")
      param$external <- NULL
    }
  }
  
  # --- Prepare parameters and model lists for tscount::tsglm.sim ---
  # Ensure param list structure matches tscount::tsglm.sim expectations
  ts_param <- list(
    intercept = param$intercept,
    past_obs = param$past_obs,
    past_mean = param$past_mean,
    xreg = param$external  # Pass external coefficients directly
  )
  
  # model list mainly defines the lags
  # Using seq_along assumes contiguous lags matching the number of coefficients
  # A more robust approach might pass model$past_obs/mean directly if they contain actual lags
  ts_model <- list(
    past_obs = if(!is.null(model$past_obs) && length(model$past_obs)>0) model$past_obs else NULL,
    past_mean = if(!is.null(model$past_mean) && length(model$past_mean)>0) model$past_mean else NULL
    # 'external = model$external' flag in model list is less critical here,
    # as presence of param$external and xreg dictates behavior in tsglm.sim
  )
  
  # --- Call tscount::tsglm.sim ---
  #cat("Calling tscount::tsglm.sim...\n") # Debug message
  result <- tryCatch({
    tscount::tsglm.sim( # Explicitly call the tscount version
      n = n,
      param = ts_param,    # Pass the correctly structured parameters
      model = ts_model,    # Pass the model structure (lags)
      xreg = xreg,         # Pass the original xreg matrix DIRECTLY
      link = link,
      distr = distr,
      distrcoefs = size,   # Pass 'size' to 'distrcoefs' for nbinom
      n_start = n_start
    )
  }, error = function(e) {
    # Provide more context on failure
    err_msg <- paste("tscount::tsglm.sim failed:", conditionMessage(e))
    #cat("Error details:\n")
    #cat(" n =", n, "\n")
    #cat(" link =", link, "\n")
    #cat(" distr =", distr, "\n")
    #cat(" size =", size, "\n")
    #cat(" n_start =", n_start, "\n")
    #print("Parameters (ts_param):")
    #print(ts_param)
    #print("Model lags (ts_model):")
    #print(ts_model)
    #print("xreg dimensions:")
    #print(dim(xreg))
    stop(err_msg, call. = FALSE) # Stop execution
  })
  # cat("tscount::tsglm.sim call successful.\n") # Debug message
  
  # --- Store additional information ---
  # Store the original input parameters/model for reference
  result$parameters_used <- param
  result$model_structure <- model
  result$link <- link
  result$distr <- distr
  result$size_used <- size
  result$n <- n
  result$n_start <- n_start
  
  # --- Set class and return --
  class(result) <- c("ingarch.sim", "tsglm.sim") # More specific classes
  return(result)
  
}