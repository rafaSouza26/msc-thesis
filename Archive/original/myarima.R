myarima <- function(x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant=TRUE, ic="aic", trace=FALSE, approximation=FALSE, offset=0, xreg=NULL, method = NULL, ...) {
  # Length of non-missing interior
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing),1)
  lastnonmiss <- tail(which(!missing),1)
  n <- sum(!missing[firstnonmiss:lastnonmiss])
  m <- frequency(x)
  use.season <- (sum(seasonal) > 0) & m > 0
  diffs <- order[2] + seasonal[2]
  if(is.null(method)){
    if (approximation) {
      method <- "CSS"
    } else {
      method <- "CSS-ML"
    }
  }
  if (diffs == 1 && constant) {
    xreg <- `colnames<-`(cbind(drift = 1:length(x), xreg),
      make.unique(c("drift", if(is.null(colnames(xreg)) && !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    if (use.season) {
      suppressWarnings(fit <- try(stats::arima(x = x, order = order, seasonal = list(order = seasonal, period = m), xreg = xreg, method = method, ...), silent = TRUE))
    } else {
      suppressWarnings(fit <- try(stats::arima(x = x, order = order, xreg = xreg, method = method, ...), silent = TRUE))
    }
  }
  else {
    if (use.season) {
      suppressWarnings(fit <- try(stats::arima(x = x, order = order, seasonal = list(order = seasonal, period = m), include.mean = constant, method = method, xreg = xreg, ...), silent = TRUE))
    } else {
      suppressWarnings(fit <- try(stats::arima(x = x, order = order, include.mean = constant, method = method, xreg = xreg, ...), silent = TRUE))
    }
  }
  if (is.null(xreg)) {
    nxreg <- 0
  } else {
    nxreg <- ncol(as.matrix(xreg))
  }
  if (!is.element("try-error", class(fit))) {
    nstar <- n - order[2] - seasonal[2] * m
    if (diffs == 1 && constant) {
      # fitnames <- names(fit$coef)
      # fitnames[length(fitnames)-nxreg] <- "drift"
      # names(fit$coef) <- fitnames
      fit$xreg <- xreg
    }
    npar <- length(fit$coef[fit$mask]) + 1
    if (method == "CSS") {
      fit$aic <- offset + nstar * log(fit$sigma2) + 2 * npar
    }
    if (!is.na(fit$aic)) {
      fit$bic <- fit$aic + npar * (log(nstar) - 2)
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
      fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic, aicc = fit$aicc)
    }
    else {
      fit$aic <- fit$bic <- fit$aicc <- fit$ic <- Inf
    }
    # Adjust residual variance to be unbiased
    fit$sigma2 <- sum(fit$residuals ^ 2, na.rm = TRUE) / (nstar - npar + 1)

    # Check for unit roots
    minroot <- 2
    if (order[1] + seasonal[1] > 0) {
      testvec <- fit$model$phi
      k <- abs(testvec) > 1e-8
      if (sum(k) > 0) {
        last.nonzero <- max(which(k))
      } else {
        last.nonzero <- 0
      }
      if (last.nonzero > 0) {
        testvec <- testvec[1:last.nonzero]
        proots <- try(polyroot(c(1,-testvec)))
        if (!is.element("try-error", class(proots))) {
          minroot <- min(minroot, abs(proots))
        }
        else fit$ic <- Inf
      }
    }
    if (order[3] + seasonal[3] > 0 & fit$ic < Inf) {
      testvec <- fit$model$theta
      k <- abs(testvec) > 1e-8
      if (sum(k) > 0) {
        last.nonzero <- max(which(k))
      } else {
        last.nonzero <- 0
      }
      if (last.nonzero > 0) {
        testvec <- testvec[1:last.nonzero]
        proots <- try(polyroot(c(1,testvec)))
        if (!is.element("try-error", class(proots))) {
          minroot <- min(minroot, abs(proots))
        }
        else fit$ic <- Inf
      }
    }
    # Avoid bad models
    if (minroot < 1 + 1e-2 | checkarima(fit)) {
      fit$ic <- Inf
    }

    fit$xreg <- xreg

    if (trace) {
      cat("\n", arima.string(fit, padding = TRUE), ":", fit$ic)
    }

    return(structure(fit, class = c("forecast_ARIMA", "ARIMA", "Arima")))
  }
  else {
    # Catch errors due to unused arguments
    if (length(grep("unused argument", fit)) > 0L) {
      stop(fit[1])
    }

    if (trace) {
      cat("\n ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
      if (use.season) {
        cat("(", seasonal[1], ",", seasonal[2], ",", seasonal[3], ")[", m, "]", sep = "")
      }
      if (constant && (order[2] + seasonal[2] == 0)) {
        cat(" with non-zero mean")
      } else if (constant && (order[2] + seasonal[2] == 1)) {
        cat(" with drift        ")
      } else if (!constant && (order[2] + seasonal[2] == 0)) {
        cat(" with zero mean    ")
      } else {
        cat("                   ")
      }
      cat(" :", Inf)
    }
    return(list(ic = Inf))
  }
}
