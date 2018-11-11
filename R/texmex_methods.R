# nobs and vocv methods for class "evmOpt"
# (texmex already has coef and logLik methods)

#' @export
nobs.evmOpt <- function(object, ...) {
  return(length(object$data$y))
}

#' @export
vcov.evmOpt <- function(object, ...) {
  vc <- object$cov
  par_names <- names(coef(object))
  dimnames(vc) <- list(par_names, par_names)
  return(vc)
}
