# ============================= extRemes::fevd ============================== #

# Methods for class "extRemes_gev"

#' @export
logLikVec.extRemes_gev <- function(object, pars = NULL, ...) {
  # Create an "fevd" object, so that we can use functions to extract
  # information about the data and the model
  fevd_object <- object
  class(fevd_object) <- "fevd"
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object.  We use bespoke methods for which the coef
  # names of nested models are consistent
  if (is.null(pars)) {
    pars <- coef(object)
  }
  # Use the datagrabber method for "fevd" to grab the data
  the_data <- distillery::datagrabber(fevd_object)
  response_data <- the_data$y
  # Determine whether or not there are covariates
  stat_model <- extRemes::is.fixedfevd(fevd_object)
  if (stat_model) {
    mu <- pars["mu0"]
    if (object$par.models$log.scale) {
      sig <- exp(pars["phi0"])
    } else {
      sig <- pars["sigma0"]
    }
    if (object$type == "GEV") {
      xi <- pars["xi0"]
    } else {
      xi <- 0
    }
  } else {
    design_matrices <- extRemes::setup.design(fevd_object)
    X.loc <- design_matrices$X.loc
    X.sc <- design_matrices$X.sc
    X.sh <- design_matrices$X.sh
    if (object$const.loc) {
      mu <- pars["mu0"]
    } else {
      n_mu <- object$results$num.pars$location
      if (n_mu == 1) {
        mu_pars <- "mu0"
      } else {
        mu_pars <- paste0("mu", 0:(ncol(X.loc) - 1))
      }
      mu <- as.vector(X.loc %*% pars[mu_pars])
    }
    if (object$const.scale) {
      if (object$par.models$log.scale) {
        sig <- exp(pars["phi0"])
      } else {
        sig <- pars["sigma0"]
      }
    } else {
      n_sig <- object$results$num.pars$scale
      if (n_sig == 1) {
        sigphi_pars <- ifelse(object$par.models$log.scale, "phi0",
                               "sigma0")
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
        xi <- pars["xi0"]
      } else {
        n_xi <- object$results$num.pars$shape
        if (n_xi == 1) {
          xi_pars <- "xi0"
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
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.extRemes_gev <- function(object, ...) {
  return(object$n)
}

#' @export
coef.extRemes_gev <- function(object, ...) {
  fevd_names <- names(object$results$par)
  which_location <- which(fevd_names == "location")
  which_scale <- which(fevd_names == "scale")
  which_logscale <- which(fevd_names == "log.scale")
  which_shape <- which(fevd_names == "shape")
  fevd_names[which_location] <- "mu0"
  fevd_names[which_scale] <- "sigma0"
  fevd_names[which_logscale] <- "phi0"
  fevd_names[which_shape] <- "xi0"
  val <- object$results$par
  names(val) <- fevd_names
  return(val)
}

#' @export
vcov.extRemes_gev <- function(object, ...) {
  temp <- object
  class(temp) <- "fevd"
  vc <- extRemes::parcov.fevd(temp)
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.extRemes_gev <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# Methods for class "extRemes_gp"

#' @export
logLikVec.extRemes_gp <- function(object, pars = NULL, ...) {
  # Create an "fevd" object, so that we can use functions to extract
  # information about the data and the model
  fevd_object <- object
  class(fevd_object) <- "fevd"
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object.  We use bespoke methods for which the coef
  # names of nested models are consistent
  if (is.null(pars)) {
    pars <- coef(object)
  }
  # Use the datagrabber method for "fevd" to grab the data
  the_data <- distillery::datagrabber(fevd_object)
  # Calculate the threshold exceedances (values above the threshold)
  exc_ind <- the_data$y > object$threshold
  response_data <- the_data$y[exc_ind]
  if (object$const.thresh) {
    the_threshold <- object$threshold
  } else {
    the_threshold <- object$threshold[exc_ind]
  }
  # Determine whether or not there are covariates
  stat_model <- extRemes::is.fixedfevd(fevd_object)
  if (stat_model) {
    if (object$par.models$log.scale) {
      sig <- exp(pars["phi0"])
    } else {
      sig <- pars["sigma0"]
    }
    if (object$type == "GP") {
      xi <- pars["xi0"]
    } else {
      xi <- 0
    }
  } else {
    design_matrices <- extRemes::setup.design(fevd_object)
    X.sc <- design_matrices$X.sc[exc_ind, , drop = FALSE]
    X.sh <- design_matrices$X.sh[exc_ind, , drop = FALSE]
    if (object$const.scale) {
      if (object$par.models$log.scale) {
        sig <- exp(pars["phi0"])
      } else {
        sig <- pars["sigma0"]
      }
    } else {
      n_sig <- object$results$num.pars$scale
      if (n_sig == 1) {
        sigphi_pars <- ifelse(object$par.models$log.scale, "phi0",
                              "sigma0")
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
        xi <- pars["xi0"]
      } else {
        n_xi <- object$results$num.pars$shape
        if (n_xi == 1) {
          xi_pars <- "xi0"
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
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.extRemes_gp <- function(object, ...) {
  return(object$n)
}

#' @export
coef.extRemes_gp <- function(object, ...) {
  fevd_names <- names(object$results$par)
  which_scale <- which(fevd_names == "scale")
  which_logscale <- which(fevd_names == "log.scale")
  which_shape <- which(fevd_names == "shape")
  fevd_names[which_scale] <- "sigma0"
  fevd_names[which_logscale] <- "phi0"
  fevd_names[which_shape] <- "xi0"
  val <- object$results$par
  names(val) <- fevd_names
  return(val)
}

#' @export
vcov.extRemes_gp <- function(object, ...) {
  temp <- object
  class(temp) <- "fevd"
  vc <- extRemes::parcov.fevd(temp)
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.extRemes_gp <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# Methods for class "extRemes_pp"

#' @export
logLikVec.extRemes_pp <- function(object, pars = NULL, ...) {
  # Create an "fevd" object, so that we can use functions to extract
  # information about the data and the model
  fevd_object <- object
  class(fevd_object) <- "fevd"
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object.  We use bespoke methods for which the coef
  # names of nested models are consistent
  if (is.null(pars)) {
    pars <- coef(object)
  }
  # Use the datagrabber method for "fevd" to grab the data
  the_data <- distillery::datagrabber(fevd_object)
  response_data <- the_data$y
  the_threshold <- object$threshold
  # Determine whether or not there are covariates
  stat_model <- extRemes::is.fixedfevd(fevd_object)
  if (stat_model) {
    mu <- pars["mu0"]
    if (object$par.models$log.scale) {
      sig <- exp(pars["phi0"])
    } else {
      sig <- pars["sigma0"]
    }
    xi <- pars["xi0"]
  } else {
    design_matrices <- extRemes::setup.design(fevd_object)
    X.loc <- design_matrices$X.loc
    X.sc <- design_matrices$X.sc
    X.sh <- design_matrices$X.sh
    if (object$const.loc) {
      mu <- pars["mu0"]
    } else {
      n_mu <- object$results$num.pars$location
      if (n_mu == 1) {
        mu_pars <- "mu0"
      } else {
        mu_pars <- paste0("mu", 0:(ncol(X.loc) - 1))
      }
      mu <- as.vector(X.loc %*% pars[mu_pars])
    }
    if (object$const.scale) {
      if (object$par.models$log.scale) {
        sig <- exp(pars["phi0"])
      } else {
        sig <- pars["sigma0"]
      }
    } else {
      n_sig <- object$results$num.pars$scale
      if (n_sig == 1) {
        sigphi_pars <- ifelse(object$par.models$log.scale, "phi0",
                              "sigma0")
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
    if (object$const.shape) {
      xi <- pars["xi0"]
    } else {
      n_xi <- object$results$num.pars$shape
      if (n_xi == 1) {
        xi_pars <- "xi0"
      } else {
        xi_pars <- paste0("xi", 0:(ncol(X.sh) - 1))
      }
      xi <- as.vector(X.sh %*% pars[xi_pars])
    }
  }
  # Calculate the (weighted) loglikelihood contributions
  if (any(sig <= 0)) {
    val <- -Inf
  } else {
    pp_loglik_vec <- function(x, u, mu, sigma, xi) {
      logFu <- revdbayes::pgev(q = u, loc = mu, scale = sig, shape = xi,
                               log.p = TRUE)
      logFx <- revdbayes::pgev(q = x, loc = mu, scale = sig, shape = xi,
                               log.p = TRUE)
      logfx <- revdbayes::dgev(x = x, loc = mu, scale = sig,
                               shape = xi, log = TRUE)
      rate_term <-  logFu / object$npy
      exc_term <- ifelse(x > u, logfx - logFx, 0)
      return(rate_term + exc_term)
    }
    val <- pp_loglik_vec(x = response_data, u = object$threshold, mu = mu,
                         sigma = sig, xi = xi)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(pars)
  # Delete the names attribute (I'm not sure from where it came) because its
  # existence causes a problem in logLik.logLikVec()
  attr(val, "names") <- NULL
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.extRemes_pp <- function(object, ...) {
  return(object$n)
}

#' @export
coef.extRemes_pp <- function(object, ...) {
  fevd_names <- names(object$results$par)
  which_location <- which(fevd_names == "location")
  which_scale <- which(fevd_names == "scale")
  which_logscale <- which(fevd_names == "log.scale")
  which_shape <- which(fevd_names == "shape")
  fevd_names[which_location] <- "mu0"
  fevd_names[which_scale] <- "sigma0"
  fevd_names[which_logscale] <- "phi0"
  fevd_names[which_shape] <- "xi0"
  val <- object$results$par
  names(val) <- fevd_names
  return(val)
}

#' @export
vcov.extRemes_pp <- function(object, ...) {
  temp <- object
  class(temp) <- "fevd"
  vc <- extRemes::parcov.fevd(temp)
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}

#' @export
logLik.extRemes_pp <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

# See extRemes_methods.R for nobs, coef, vcov, logLik methods for class "fevd"
