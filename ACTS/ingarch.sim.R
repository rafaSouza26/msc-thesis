ingarch.sim <- function(n,
                        param = list(intercept = 1, past_obs = NULL, past_mean = NULL, xreg = NULL),
                        model = list(past_obs = NULL, past_mean = NULL, external = FALSE),
                        xreg = NULL, 
                        link = c("identity", "log"),
                        distr = c("poisson", "nbinom"),
                        distrcoefs = NULL, # Changed: size to distrcoefs
                        fit = NULL,
                        n_start = 50) {
  
  # Store original arguments for checking if they were supplied, before 'fit' overrides them
  arg_param_missing <- missing(param)
  original_param_arg <- param
  arg_model_missing <- missing(model)
  original_model_arg <- model
  arg_link_missing <- missing(link)
  original_link_arg <- link
  arg_distr_missing <- missing(distr)
  original_distr_arg <- distr
  arg_distrcoefs_missing <- missing(distrcoefs)
  original_distrcoefs_arg <- distrcoefs
  arg_xreg_missing <- missing(xreg)
  original_xreg_arg <- xreg
  
  # --- Input Validation (basic) ---
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(n_start) || length(n_start) != 1 || n_start < 0 || n_start != floor(n_start)) {
    stop("n_start must be a non-negative integer")
  }
  
  # --- Parameters from fit object ---
  if (!is.null(fit)) {
    if (!inherits(fit, "tsglm")) { # Assuming 'tsglm' is the class of the fit object
      stop("fit must be an object of class 'tsglm'")
    }
    tryCatch({
      # Model
      if (!arg_model_missing && !identical(original_model_arg, fit$model)) {
        warning("Argument 'model' is overridden by 'fit$model'")
      }
      model <- fit$model
      
      # Link
      if (!arg_link_missing && !identical(original_link_arg, fit$link)) {
        warning("Argument 'link' is overridden by 'fit$link'")
      }
      link <- fit$link
      
      # Distr
      if (!arg_distr_missing && !identical(original_distr_arg, fit$distr)) {
        warning("Argument 'distr' is overridden by 'fit$distr'")
      }
      distr <- fit$distr
      
      # Distrcoefs
      if (!arg_distrcoefs_missing && !is.null(original_distrcoefs_arg) && !identical(original_distrcoefs_arg, fit$distrcoefs)) {
        warning("Argument 'distrcoefs' is overridden by 'fit$distrcoefs'")
      }
      distrcoefs <- fit$distrcoefs # Directly use from fit
      
      # Coefficients to form 'param' list
      if (!arg_param_missing && !is.null(original_param_arg)) {
        warning("Argument 'param' is overridden by parameters from 'fit$coefficients'")
      }
      coefs <- stats::coefficients(fit) # Use stats::coefficients for safety
      names_coefs <- names(coefs)
      intercept_ind <- which(names_coefs == "(Intercept)" | names_coefs == "intercept") # Allow for "(Intercept)"
      past_obs_inds <- which(grepl("^beta", names_coefs))
      past_mean_inds <- which(grepl("^alpha", names_coefs))
      
      # Determine external regressor coefficient indices
      # These are coefficients not matching intercept, beta, or alpha
      known_inds <- c(intercept_ind, past_obs_inds, past_mean_inds)
      external_inds <- which(! (1:length(coefs) %in% known_inds) )
      
      param <- list(
        intercept = if(length(intercept_ind) > 0) coefs[intercept_ind] else 0, # Default to 0 if no intercept in fit
        past_obs = if(length(past_obs_inds) > 0) coefs[past_obs_inds] else NULL,
        past_mean = if(length(past_mean_inds) > 0) coefs[past_mean_inds] else NULL,
        xreg = if(length(external_inds) > 0) coefs[external_inds] else NULL # Changed: external to xreg
      )
      if (length(intercept_ind) == 0 && ("(Intercept)" %in% names_coefs || "intercept" %in% names_coefs) ){
        # Should have found intercept if it's named commonly
        warning("Could not identify intercept coefficient from 'fit' object by name.")
      }
      
      
      # xreg data matrix handling when 'fit' is provided
      # param$xreg are coefficients from fit. model$external is structure from fit.
      # model$external might be a logical vector or single logical.
      # The presence of coefficients in param$xreg from fit indicates external regressors were used.
      if (length(param$xreg) > 0 || (is.logical(model$external) && any(model$external)) ) { # If fit implies external regressors
        if (arg_xreg_missing || is.null(original_xreg_arg)) { # and user did not provide xreg data
          if (!is.null(fit$xreg)) {
            if (nrow(fit$xreg) >= n) {
              xreg <- fit$xreg[1:n, , drop = FALSE]
            } else {
              stop("Provided 'fit' object's xreg has fewer rows (", nrow(fit$xreg), ") than requested simulation length 'n' (", n, ").")
            }
          } else {
            stop("'fit' object indicates external regressors but 'fit$xreg' is NULL.")
          }
        } # else: user provided xreg, it will be used. Validation will occur later.
      } else { # No external regressors from fit's parameters/model
        if (!arg_xreg_missing && !is.null(original_xreg_arg)){
          warning("xreg data was provided by user, but 'fit' object's parameters/model do not indicate use of external regressors. User-supplied xreg will be checked for consistency.")
        } # If user didn't supply xreg either, xreg remains NULL, which is fine.
      }
    }, error = function(e) {
      stop("Failed to extract parameters/data from fit object: ", conditionMessage(e))
    })
  }
  
  # --- Validate link and distribution (match.arg ensures they are one of the allowed options) ---
  link <- match.arg(link)
  distr <- match.arg(distr)
  
  # --- Validate parameters (basic checks) ---
  if (!is.list(param)) stop("param must be a list")
  if (!is.numeric(param$intercept) || length(param$intercept) != 1) {
    stop("param$intercept must be a single numeric value")
  }
  if (!is.null(param$past_obs) && !is.numeric(param$past_obs)) {
    stop("param$past_obs coefficients must be numeric")
  }
  if (!is.null(param$past_mean) && !is.numeric(param$past_mean)) {
    stop("param$past_mean coefficients must be numeric")
  }
  if (!is.null(param$xreg) && !is.numeric(param$xreg)) { # Changed: external to xreg
    stop("param$xreg coefficients must be numeric")
  }
  
  # --- Validate distribution parameters (distrcoefs) ---
  if (distr == "nbinom") {
    if (is.null(distrcoefs)) {
      stop("distrcoefs parameter must be provided for negative binomial distribution")
    }
    if (!is.numeric(distrcoefs) || length(distrcoefs) != 1 || distrcoefs <= 0) {
      stop("distrcoefs parameter must be a single positive numeric value for nbinom distribution")
    }
  } else {
    if (!is.null(distrcoefs)) {
      warning(paste0("distrcoefs was provided for '", distr, "' distribution, but it is only used for 'nbinom'. Ignoring distrcoefs."))
    }
    distrcoefs <- NULL # Ensure distrcoefs is NULL if not nbinom
  }
  
  # --- Validate model structure (basic checks) ---
  if(!is.list(model)) stop("model must be a list")
  if(!is.null(model$past_obs) && (!is.numeric(model$past_obs) || any(model$past_obs <= 0) || any(model$past_obs != floor(model$past_obs)))){
    stop("model$past_obs must be NULL or a vector of positive integers (lags)")
  }
  if(!is.null(model$past_mean) && (!is.numeric(model$past_mean) || any(model$past_mean <= 0) || any(model$past_mean != floor(model$past_mean)))){
    stop("model$past_mean must be NULL or a vector of positive integers (lags)")
  }
  # model$external should be a logical. tscount::tsglm.sim handles vector logicals for model$external.
  # For ingarch.sim's own validation, we check based on coefficients and xreg matrix.
  if(!is.null(model$external) && !is.logical(model$external)){
    stop("model$external must be NULL or a logical value/vector.")
  }
  
  
  # --- Validate external regressors xreg ---
  n_external_coefs <- length(param$xreg) # Uses param$xreg
  # Determine if external regressors are effectively part of the model
  # model$external could be a single logical or a vector from a 'fit' object.
  # 'any(model$external)' handles vector case. If model$external is NULL (e.g. not in fit$model), it's FALSE here.
  model_defines_external <- !is.null(model$external) && any(model$external) 
  
  if (model_defines_external || n_external_coefs > 0) {
    if (is.null(xreg)) {
      stop("External regressors defined in model/parameters (model$external is TRUE or param$xreg has coefficients), but 'xreg' matrix was not provided.")
    }
    if (!is.numeric(xreg) && !is.data.frame(xreg)) { # Allow data.frame for xreg
      stop("xreg must be a numeric matrix, vector, or data.frame")
    }
    xreg <- as.matrix(xreg) # Ensure it's a matrix
    if (nrow(xreg) != n) {
      stop(sprintf("Number of rows in xreg (%d) must match n (%d)", nrow(xreg), n))
    }
    if (ncol(xreg) != n_external_coefs) {
      stop(sprintf("Number of columns in xreg (%d) must match number of external coefficients in param$xreg (%d)", ncol(xreg), n_external_coefs))
    }
    if(anyNA(xreg)){
      stop("Missing values (NA) are not allowed in xreg.")
    }
  } else { # No external effects expected from parameters or model structure
    if (!is.null(xreg)) {
      warning("xreg provided but no external coefficients found in param$xreg and model$external is FALSE. xreg will be ignored.")
      xreg <- NULL # Set to NULL as it won't be used by tscount::tsglm.sim if param$xreg is NULL and model$external is FALSE
    }
    if (!is.null(param$xreg)) { # param$xreg had coefficients but model$external was FALSE
      warning("param$xreg coefficients provided but model$external=FALSE. These coefficients might be ignored by the underlying simulation if model$external strictly dictates usage.")
      # tscount::tsglm.sim primarily uses length(param$xreg) to determine 'r'.
      # If model$external is FALSE but param$xreg is not NULL, tscount still expects matching ncol(xreg).
      # For clarity, if model$external is FALSE, we should ensure param$xreg is also NULL or that xreg is NULL.
      # The safest is to ensure xreg is NULL if no external component.
      # Or rely on tscount's internal logic: if r > 0 (from param$xreg), it expects xreg.
      # The current logic: if n_external_coefs > 0, xreg must be provided. If model_defines_external is also false, this is a contradiction.
      # Let's refine: if model$external is FALSE, then param$xreg should ideally be NULL.
      if (!model_defines_external && n_external_coefs > 0){
        warning("param$xreg coefficients are present, but model$external is FALSE. Forcing param$xreg to NULL and ignoring xreg to ensure no external component as per model$external=FALSE.")
        param$xreg <- NULL
        xreg <- NULL
      }
    }
  }
  
  # --- Prepare parameters and model lists for tscount::tsglm.sim ---
  ts_param <- list(
    intercept = param$intercept,
    past_obs = param$past_obs,
    past_mean = param$past_mean,
    xreg = param$xreg  # Changed: Pass param$xreg (coefficients)
  )
  
  # model list for tscount::tsglm.sim
  ts_model <- list(
    past_obs = if(!is.null(model$past_obs) && length(model$past_obs)>0) model$past_obs else NULL,
    past_mean = if(!is.null(model$past_mean) && length(model$past_mean)>0) model$past_mean else NULL,
    external = if(!is.null(model$external)) model$external else FALSE # Pass model$external flag/vector
  )
  
  # --- Call tscount::tsglm.sim ---
  result_sim <- tryCatch({
    tscount::tsglm.sim(
      n = n,
      param = ts_param,
      model = ts_model,
      xreg = xreg,         # Pass the xreg data matrix
      link = link,
      distr = distr,
      distrcoefs = distrcoefs, # Changed: Pass 'distrcoefs'
      n_start = n_start
    )
  }, error = function(e) {
    err_msg <- paste("tscount::tsglm.sim failed:", conditionMessage(e))
    # Optional: Add more debug info here if needed
    # print(list(n=n, param=ts_param, model=ts_model, xreg_dim=dim(xreg), link=link, distr=distr, distrcoefs=distrcoefs, n_start=n_start))
    stop(err_msg, call. = FALSE)
  })
  
  # --- Store additional information ---
  # Create a list to return, similar to tsglm.sim but can add more specific ingarch.sim info
  final_result <- list()
  final_result$ts <- result_sim$ts
  final_result$linear.predictors <- result_sim$linear.predictors
  final_result$xreg.effects <- result_sim$xreg.effects # If tscount::tsglm.sim returns it
  
  final_result$parameters_used <- param # param now has $xreg
  final_result$model_structure <- model # model now has $external
  final_result$link_used <- link
  final_result$distr_used <- distr
  final_result$distrcoefs_used <- distrcoefs # Changed: size_used to distrcoefs_used
  final_result$n <- n
  final_result$n_eff <- length(result_sim$ts) # Effective number of simulated observations
  final_result$n_start <- n_start
  final_result$xreg_input <- xreg # Store the xreg matrix used for simulation
  if (!is.null(fit)) final_result$fit_object_used <- TRUE
  
  
  class(final_result) <- c("ingarch.sim", "tsglm.sim.result") # Use a distinct class, maybe "tsglm.sim.result" if structure is very similar
  return(final_result)
}