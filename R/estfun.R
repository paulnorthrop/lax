# Create a default estfun method and, for safety, create individual methods
# for all the classes currently involved in lax.  At the moment they all
# use numDeriv::jacobian() to do the calculation

#' @export
estfun.default <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# eva

#' @export
estfun.laxeva_gpd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.laxeva_rlarg <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# evd

#' @export
estfun.evd_fgev <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.evd_fpot <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# evir

#' @export
estfun.evir_gev <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.evir_gpd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.evir_pot <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# extRemes

#' @export
estfun.fevd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# fExtremes

#' @export
estfun.fExtremes_gev <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.fExtremes_gpd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# ismev

#' @export
estfun.ismev_gev <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.ismev_pp <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.ismev_gpd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.ismev_rlarg <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# mev

#' @export
estfun.laxmev_gev <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.laxmev_pp <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.laxmev_gpd <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.laxmev_rlarg <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

#' @export
estfun.laxmev_egp <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# texmex

#' @export
estfun.texmex_evmOpt <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}
