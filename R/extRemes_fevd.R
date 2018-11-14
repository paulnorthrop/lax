# NOTE: needs to deal with the -1 case in the formula.
#       check what this does!
#> fitPORTstdmax <- fevd(PORTw$TMX1, PORTw, scale.fun=~STDTMAX, use.phi=TRUE)
#> coef(fitPORTstdmax)
#location       phi0       phi1      shape
#14.9413877  0.1177255  0.2062192 -0.3677961
#> fitPORTstdmax <- fevd(PORTw$TMX1, PORTw, scale.fun=~STDTMAX-1, use.phi=TRUE)
#> coef(fitPORTstdmax)
#location  log.scale      shape
#14.9275962  0.2294392 -0.3732653

# Deal with names and getting correct covariate values

# ============================= extRemes::fevd ============================== #

# Methods for class "fevd"

#' @export
logLikVec.fevd <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # Use a model-dependent logLikVec method
  if (object$type == "GEV" || object$type == "Gumbel") {
    val <- fevd_gev_logLikVec(object, pars, ...)
  } else if (object$type == "GP" || object$type == "Exponential") {
    val <- fevd_gp_logLikVec(object, pars, ...)
  } else if (object$type == "PP") {
    val <- fevd_pp_logLikVec(object, pars, ...)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

fevd_gev_logLikVec <- function(object, pars = NULL, ...) {
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  # Use the datagrabber method for "fevd" to grab the data
  the_data <- distillery::datagrabber(object)
  response_data <- the_data$y
  # Determine whether or not there are covariates
  stat_model <- extRemes::is.fixedfevd(object)
  if (stat_model) {
    mu <- pars["location"]
    if (object$par.models$log.scale) {
      sig <- exp(pars["log.scale"])
    } else {
      sig <- pars["scale"]
    }
    if (object$type == "GEV") {
      xi <- pars["shape"]
    } else {
      xi <- 0
    }
  } else {
    design_matrices <- extRemes::setup.design(object)
    X.loc <- design_matrices$X.loc
    X.sc <- design_matrices$X.sc
    X.sh <- design_matrices$X.sh
    if (object$const.loc) {
      mu <- pars["location"]
    } else {
      n_mu <- object$results$num.pars$location
      if (n_mu == 1) {
        mu_pars <- "location"
      } else {
        mu_pars <- paste0("mu", 0:(ncol(X.loc) - 1))
      }
      mu <- as.vector(X.loc %*% pars[mu_pars])
    }
    if (object$const.scale) {
      if (object$par.models$log.scale) {
        sig <- exp(pars["log.scale"])
      } else {
        sig <- pars["scale"]
      }
    } else {
      n_sig <- object$results$num.pars$scale
      if (n_sig == 1) {
        sigphi_pars <- ifelse(object$par.models$log.scale, "log.scale",
                               "scale")
        sig <- as.vector(X.sc %*% pars[sigphi_pars])
        if (object$par.models$log.scale) {
          sig <- exp(sig)
        }
      } else {
        if (object$par.models$log.scale) {
          phi_pars <- paste0("phi", 0:(ncol(X.sc) - 1))
          phi <- as.vector(X.sc %*% pars[phi_pars])
          sig <- exp(phi)
        } else {
          sig_pars <- paste0("sigma", 0:(ncol(X.sc) - 1))
          sig <- as.vector(X.sc %*% pars[sig_pars])
        }
      }
    }
    if (object$type == "GEV") {
      if (object$const.shape) {
        xi <- pars["shape"]
      } else {
        n_xi <- object$results$num.pars$shape
        if (n_xi == 1) {
          xi_pars <- "shape"
        } else {
          xi_pars <- paste0("xi", 0:(ncol(X.sh) - 1))
        }
        xi <- as.vector(X.sh %*% pars[xi_pars])
      }
    } else {
      xi <- 0
    }
  }
  # Calculate the (weighted) loglikelihood contributions
  if (any(sig <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgev(response_data, loc = mu, scale = sig,
                           shape = xi, log = TRUE) * object$weights
  }
  return(val)
}

fevd_gp_logLikVec <- function(object, pars = NULL, ...) {
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  # Use the datagrabber method for "fevd" to grab the data
  the_data <- distillery::datagrabber(object)
  # Calculate the threshold exceedances (values above the threshold)
  response_data <- the_data$y[the_data$y > object$threshold]
  if (object$const.thresh) {
    the_threshold <- object$threshold
  } else {
    the_threshold <- object$threshold[the_data$y > object$threshold]
  }
  # Determine whether or not there are covariates
  stat_model <- extRemes::is.fixedfevd(object)
  if (stat_model) {
    if (object$par.models$log.scale) {
      sig <- exp(pars["log.scale"])
    } else {
      sig <- pars["scale"]
    }
    if (object$type == "GP") {
      xi <- pars["shape"]
    } else {
      xi <- 0
    }
  } else {
    design_matrices <- extRemes::setup.design(object)
    X.sc <- design_matrices$X.sc
    X.sh <- design_matrices$X.sh
    if (object$const.scale) {
      if (object$par.models$log.scale) {
        sig <- exp(pars["log.scale"])
      } else {
        sig <- pars["scale"]
      }
    } else {
      n_sig <- object$results$num.pars$scale
      if (n_sig == 1) {
        sigphi_pars <- ifelse(object$par.models$log.scale, "log.scale",
                              "scale")
        sig <- as.vector(X.sc %*% pars[sigphi_pars])
        if (object$par.models$log.scale) {
          sig <- exp(sig)
        }
      } else {
        if (object$par.models$log.scale) {
          phi_pars <- paste0("phi", 0:(ncol(X.sc) - 1))
          phi <- as.vector(X.sc %*% pars[phi_pars])
          sig <- exp(phi)
        } else {
          sig_pars <- paste0("sigma", 0:(ncol(X.sc) - 1))
          sig <- as.vector(X.sc %*% pars[sig_pars])
        }
      }
    }
    if (object$type == "GP") {
      if (object$const.shape) {
        xi <- pars["shape"]
      } else {
        n_xi <- object$results$num.pars$shape
        if (n_xi == 1) {
          xi_pars <- "shape"
        } else {
          xi_pars <- paste0("xi", 0:(ncol(X.sh) - 1))
        }
        xi <- as.vector(X.sh %*% pars[xi_pars])
      }
    } else {
      xi <- 0
    }
  }
  # Calculate the (weighted) loglikelihood contributions
  if (any(sig <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = the_threshold, scale = sig,
                          shape = xi, log = TRUE) * object$weights
  }
  return(val)
}

# See extRemes_methods.R for nobs, coef, vcov, logLik methods for class "fevd"
