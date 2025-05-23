source("./ACTS/newmodel.R")      # Ensure this path is correct
source("./ACTS/ingarch.string.R") # Ensure this path is correct
source("./ACTS/search.ingarch.R") # Ensure this path is correct

auto.ingarch <- function(y,
                         max.p = 7,
                         max.q = 7,
                         max.order = 5, # Still an argument, primarily for stepwise=FALSE path
                         start.p = 2,
                         start.q = 2,
                         distribution = c("poisson", "nbinom"),
                         link = c("log", "identity"),
                         xreg = NULL,
                         ic = c("aicc", "aic", "bic", "qic"),
                         stepwise = TRUE,
                         nmodels = 94,
                         trace = FALSE,
                         show_warnings = FALSE,
                         method = NULL, 
                         parallel = FALSE,
                         num.cores = 2,
                         x = y, 
                         ...) {
  
  # Initial validation
  if (stepwise && parallel) {
    warning("Parallel computing is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  
  if (trace && parallel && stepwise == FALSE) {
    message("Tracing model searching in parallel (stepwise=FALSE) might produce interleaved output.")
  }
  
  series <- deparse(substitute(y))
  x <- as.ts(x) 
  
  # Parameter validation for max.order (still good to have for the arg itself)
  if (!is.numeric(max.order) || max.order < 0) {
    warning("max.order must be a non-negative integer. Setting to default (5).")
    max.order <- 5
  }
  
  if (!is.null(start.p) && start.p > max.p) {
    warning("start.p cannot be greater than max.p, setting start.p = max.p")
    start.p <- max.p
  }
  if (!is.null(start.q) && start.q > max.q) { 
    warning("start.q cannot be greater than max.q, setting start.q = max.q")
    start.q <- max.q
  }
  
  if (NCOL(x) > 1) {
    stop("auto.ingarch can only handle univariate time series")
  }
  
  if (!all(x[!is.na(x)] >= 0) || !all(x[!is.na(x)] == floor(x[!is.na(x)]))) {
    stop("INGARCH models require non-negative integer (count) data")
  }
  
  orig.x <- x
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing), 1)
  serieslength <- sum(!missing[firstnonmiss:length(x)]) 
  
  if(length(firstnonmiss) > 0 && firstnonmiss > 1){
    x <- window(x, start=time(x)[firstnonmiss])
    if(!is.null(xreg)){
      xreg <- xreg[firstnonmiss:NROW(xreg), , drop=FALSE] 
    }
  }
  
  if(!is.null(xreg)) {
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if(NROW(xreg) != length(x))
      stop(paste("Number of rows in xreg (", NROW(xreg), ") must match length of time series (", length(x), ") after NA handling.", sep=""))
    if(any(is.na(xreg)))
      stop("Missing values in external regressors are not allowed after initial NA trim of y.")
    if(is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("xreg", 1:NCOL(xreg))
    }
    link_matched_arg <- match.arg(link) 
    if(link_matched_arg == "identity" && any(xreg < 0)) { 
      stop("When using 'identity' link, all values in xreg must be non-negative")
    }
  }
  
  if (all(x == x[which(!is.na(x))[1]], na.rm = TRUE)) { # Check against first non-NA value
    if(all(is.na(x))) 
      stop("All data are missing")
    fit <- myingarch(x, order = c(0, 0), ic = match.arg(ic),
                     distr = match.arg(distribution), link = match.arg(link), xreg = xreg, ...)
    fit$x <- orig.x
    fit$series <- series
    fit$call <- match.call()
    return(fit)
  }
  
  ic <- match.arg(ic)
  distribution_matched <- match.arg(distribution) 
  link_matched <- match.arg(link)                 
  
  # Adjust max.p, max.q based on serieslength
  max.p <- min(max.p, floor(serieslength / 3))
  max.q <- min(max.q, floor(serieslength / 3))
  # Note: The line `max.order <- min(max.order, max.p + max.q)` was removed in previous versions
  # and stays removed, as max.order should be an independent constraint.
  
  if (serieslength <= 3L) { 
    ic <- "aic"
  }
  
  if (!stepwise) {
    bestfit <- search.ingarch(
      x = x, 
      max.p = max.p, 
      max.q = max.q,
      max.order = max.order, # << max.order is passed to search.ingarch
      distribution = distribution_matched,
      link = link_matched,
      xreg = xreg,
      ic = ic, 
      trace = trace, 
      show_warnings = show_warnings,
      parallel = parallel, 
      num.cores = num.cores, 
      ...
    )
    
    bestfit$x <- orig.x
    bestfit$series <- series
    bestfit$call <- match.call() 
    
    if (!is.null(bestfit$results) && is.list(bestfit$results) && !is.matrix(bestfit$results)) {
      temp_results_df <- tryCatch(do.call(rbind, lapply(bestfit$results, function(item) {
        if(is.list(item) && all(c("p", "q", "ic") %in% names(item))) {
          data.frame(p = item$p, q = item$q, ic = item$ic)
        } else { NULL }
      })), error = function(e) NULL)
      if(!is.null(temp_results_df) && nrow(temp_results_df) > 0) {
        bestfit$results <- as.matrix(temp_results_df)
      } else {
        bestfit$results <- matrix(numeric(0), ncol = 3, dimnames = list(NULL, c("p", "q", "ic")))
      }
    } else if (is.null(bestfit$results)){
      bestfit$results <- matrix(numeric(0), ncol = 3, dimnames = list(NULL, c("p", "q", "ic")))
    }
    
    if (trace && !is.null(bestfit$model)) { 
      if(exists("ingarch.string", mode="function")) {
        cat("\n\n Best model from search.ingarch:", ingarch.string(bestfit, padding = TRUE), "\n\n")
      } else {
        cat("\n\n Best model from search.ingarch: p=", bestfit$model$past_obs, " q=", bestfit$model$past_mean, "\n\n")
      }
    }
    if (!is.null(bestfit$coefficients)) {
      bestfit$coefficients[abs(bestfit$coefficients) < 1e-4] <- 1e-4 * sign(bestfit$coefficients[abs(bestfit$coefficients) < 1e-4])
    }
    return(bestfit)
  } 
  
  # --- Stepwise Path from here ---
  
  if (length(x) < 10L) { 
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  p_current <- min(start.p, max.p)
  q_current <- min(start.q, max.q)
  
  model_warnings <- list()
  results <- matrix(NA, nrow = nmodels, ncol = 3) 
  colnames(results) <- c("p", "q", "ic")
  model_count <- 0
  
  warnings_captured <- character(0)
  bestfit <- withCallingHandlers({
    model_count <- model_count + 1
    if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(%d,%d)", model_count, p_current, q_current))
    myingarch(x, order = c(p_current, q_current), ic = ic, trace = FALSE, xreg = xreg, 
              distr = distribution_matched, link = link_matched, ...)
  }, warning = function(w) { warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
  results[model_count, ] <- c(p_current, q_current, bestfit$ic)
  model_warnings[[paste0("p", p_current, "_q", q_current)]] <- warnings_captured
  k <- model_count
  
  if ((p_current != 0 || q_current != 0)) { # Check (0,0) if different
    if (k < nmodels) {
      k <- k + 1
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(0,0)", model_count))
        myingarch(x, order = c(0,0), ic = ic, trace = FALSE, xreg=xreg, distr=distribution_matched, link=link_matched, ...)
      }, warning = function(w) {warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
      results[k, ] <- c(0,0, fit$ic)
      model_warnings[[paste0("p0_q0")]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p_current <- 0; q_current <- 0; }
    }
  }
  
  # << REMOVED: (p_val + q_val <= max.order) from the conditions below >>
  if (max.p > 0 && (p_current != 1 || q_current != 0)) {
    if (k < nmodels) {
      k <- k + 1
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(1,0)", model_count))
        myingarch(x, order = c(1,0), ic = ic, trace = FALSE, xreg=xreg, distr=distribution_matched, link=link_matched, ...)
      }, warning = function(w) {warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
      results[k, ] <- c(1,0, fit$ic)
      model_warnings[[paste0("p1_q0")]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p_current <- 1; q_current <- 0; }
    }
  }
  
  if (max.q > 0 && (p_current != 0 || q_current != 1)) {
    if (k < nmodels) {
      k <- k + 1
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(0,1)", model_count))
        myingarch(x, order = c(0,1), ic = ic, trace = FALSE, xreg=xreg, distr=distribution_matched, link=link_matched, ...)
      }, warning = function(w) {warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
      results[k, ] <- c(0,1, fit$ic)
      model_warnings[[paste0("p0_q1")]] <- warnings_captured
      if (fit$ic < bestfit$ic) { bestfit <- fit; p_current <- 0; q_current <- 1; }
    }
  }
  
  if (max.p >= 1 && max.q >= 1 && (p_current != 1 || q_current != 1)) {
    if (k < nmodels) {
      k <- k + 1
      warnings_captured <- character(0)
      fit_11 <- withCallingHandlers({
        model_count <- model_count + 1
        if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(1,1)", model_count))
        myingarch(x, order = c(1,1), ic = ic, trace = FALSE, xreg=xreg, distr=distribution_matched, link=link_matched, ...)
      }, warning = function(w) {warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
      results[k, ] <- c(1,1, fit_11$ic)
      model_warnings[[paste0("p1_q1")]] <- warnings_captured
      if (fit_11$ic < bestfit$ic) { bestfit <- fit_11; p_current <- 1; q_current <- 1; }
    }
  }
  
  startk <- 0 
  while (startk < k && k < nmodels) {
    startk <- k 
    
    candidate_orders <- list(
      c(p_current - 1, q_current), c(p_current + 1, q_current),
      c(p_current, q_current - 1), c(p_current, q_current + 1),
      c(p_current - 1, q_current - 1), c(p_current - 1, q_current + 1),
      c(p_current + 1, q_current - 1), c(p_current + 1, q_current + 1)
    )
    
    for (cand_order in candidate_orders) {
      p_cand <- cand_order[1]
      q_cand <- cand_order[2]
      
      if (p_cand < 0 || q_cand < 0 || p_cand > max.p || q_cand > max.q) next 
      # << REMOVED: if (p_cand + q_cand > max.order) next >>
      if (!newmodel(p_cand, q_cand, results[1:k, , drop = FALSE])) next 
      
      if (k >= nmodels) break 
      
      k <- k + 1
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        model_count <- model_count + 1
        if(trace) cat(sprintf("\n%2d: Evaluating INGARCH(%d,%d)", model_count, p_cand, q_cand))
        myingarch(x, order = c(p_cand, q_cand), ic = ic, trace = FALSE, xreg=xreg, distr=distribution_matched, link=link_matched, ...)
      }, warning = function(w) {warnings_captured <<- c(warnings_captured, conditionMessage(w)); invokeRestart("muffleWarning")})
      
      results[k, ] <- c(p_cand, q_cand, fit$ic)
      model_warnings[[paste0("p", p_cand, "_q", q_cand)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p_current <- p_cand 
        q_current <- q_cand
      }
    } 
    if (k >= nmodels && startk < k) { 
      if(trace) cat(sprintf("\nModel limit %d reached during stepwise iteration.", nmodels))
      break 
    }
  } 
  
  if(k >= nmodels && trace){
    cat(sprintf("\nStepwise search was stopped after %d models due to reaching the model number limit: `nmodels = %i`", model_count, nmodels))
  }
  
  if (is.infinite(bestfit$ic)) { 
    stop("No suitable INGARCH model found (all resulted in Inf IC or errors).")
  }
  
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic_value <- bestfit$ic 
  bestfit$ic_name <- ic          
  bestfit$ic <- NULL             
  bestfit$call <- match.call()   
  
  bestfit$results_summary <- results[1:k, , drop = FALSE] 
  bestfit$all_warnings <- model_warnings 
  
  best_p_val <- NA; best_q_val <- NA
  if(!is.null(bestfit$model)) { 
    best_p_val <- if(is.null(bestfit$model$past_obs) || length(bestfit$model$past_obs)==0) 0 else max(bestfit$model$past_obs)
    best_q_val <- if(is.null(bestfit$model$past_mean) || length(bestfit$model$past_mean)==0) 0 else max(bestfit$model$past_mean)
  }
  if(!is.na(best_p_val) && !is.na(best_q_val)){
    bestfit$best_model_warnings <- model_warnings[[paste0("p", best_p_val, "_q", best_q_val)]]
  } else {
    bestfit$best_model_warnings <- character(0)
  }
  
  bestfit$n_models_evaluated <- model_count
  
  if (trace) {
    if(exists("ingarch.string", mode="function")) {
      cat("\n\nBest model:", ingarch.string(bestfit, padding = TRUE), "\n")
    } else if (!is.null(bestfit$model)) {
      cat("\n\nBest model: INGARCH(", best_p_val, ",", best_q_val, ") with ",bestfit$ic_name," = ", bestfit$ic_value, "\n", sep="")
    }
    cat("Total models evaluated:", model_count, "\n")
  }
  
  if (show_warnings && k > 0) {
    cat("\nWarnings for all models tested:\n")
    cat("------------------------------\n")
    
    valid_results_df <- as.data.frame(results[1:k, , drop = FALSE])
    valid_results_df <- valid_results_df[!is.na(valid_results_df$ic), , drop=FALSE] 
    
    if(nrow(valid_results_df) > 0) {
      ordered_results <- valid_results_df[order(valid_results_df$ic), , drop = FALSE]
      
      for (i_row in 1:nrow(ordered_results)) {
        p_val <- ordered_results[i_row, "p"]
        q_val <- ordered_results[i_row, "q"]
        ic_val <- ordered_results[i_row, "ic"]
        
        is_best <- (!is.na(best_p_val) && !is.na(best_q_val) && p_val == best_p_val && q_val == best_q_val && abs(ic_val - bestfit$ic_value) < 1e-9)
        model_label <- if(is_best) "--- ★ BEST MODEL" else ""
        
        cat(sprintf("INGARCH(%d,%d) with %s=%.4f %s\n", 
                    p_val, q_val, ic, ic_val, model_label)) 
        
        warnings_key <- paste0("p", p_val, "_q", q_val)
        if (warnings_key %in% names(model_warnings) && length(model_warnings[[warnings_key]]) > 0) {
          for (w_msg in model_warnings[[warnings_key]]) {
            cat("  - ", w_msg, "\n", sep="")
          }
        } else {
          cat("  No warnings\n")
        }
        cat("\n")
      }
    } else {
      cat("No models with valid ICs to display warnings for.\n")
    }
  }
  
  # Removed the coefficient adjustment line as its appropriateness for INGARCH is uncertain
  # if (!is.null(bestfit$coefficients)) {
  #   bestfit$coefficients[abs(bestfit$coefficients) < 1e-4] <- 1e-4 * sign(bestfit$coefficients[abs(bestfit$coefficients) < 1e-4])
  # }
  
  return(bestfit)
}