# coef and nobs methods for class "evd"
# (evd already has vcov and logLik methods)
# In evir_methods.R there are coef, nobs, vcov and logLik methods for
# class "gev" (provided for package evir)
# evd::fgev objects have class c("gev", "uvevd", "evd")
# Therefore, the "gev" methods in evir_methods.R redirect to the evd methods
# if the object inherits from class "evd"

#' @export
nobs.evd <- function(object, ...) {
  if (inherits(object, "gev")) {
    x <- object
    class(x) <- "evd_fgev"
    val <- nobs(x, ...)
  } else if (inherits(object, "pot")) {
    x <- object
    class(x) <- "evd_fpot"
    val <- nobs(x, ...)
  } else {
    # object$n also works for forder(), fextreme(), fbvevd() and fbvpot()
    val <- object$n
  }
  return(val)
}

#' @export
coef.evd <- function(object, complete = FALSE, ...) {
  # No 'complete' argument, for consistency with the "evd" vcov method in evd
  # ... even though evd::fgev and evd::fpot do allow parameters to be fixed
  # and return all parameter values (free and fixed) in object$param
  return(object$estimate)
}
