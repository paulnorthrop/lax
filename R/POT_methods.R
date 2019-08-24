# nobs, vcov  methods for class "pot", returned by POT::fitgpd()
# POT provides coef and logLik methods

#' @export
nobs.pot <- function(object, ...) {
  return(object$nat)
}

#' @export
vcov.pot <- function(object, complete = FALSE, ...) {
  # Variance-covariance matrix for the free parameters
  vc <- object$var.cov
  free_pars <- names(coef(object, complete = FALSE))
  dimnames(vc) <- list(free_pars, free_pars)
  if (complete) {
    all_pars <- names(coef(object, complete = TRUE))
    np <- length(all_pars)
    dummy <- matrix(NA, np, np)
    which_free <- which(all_pars %in% free_pars)
    dummy[which_free, which_free] <- vc
    vc <- dummy
    dimnames(vc) <- list(all_pars, all_pars)
  }
  return(vc)
}
