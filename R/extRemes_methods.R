# Methods for class "fevd", returned by extRemes::fevd()

#' @export
nobs.fevd <- function(object, ...) {
  return(object$n)
}

#' @export
coef.fevd <- function(object, ...) {
  return(object$results$par)
}

#' @export
vcov.fevd <- function(object, ...) {
  return(extRemes::parcov.fevd(object))
}

#' @export
logLik.fevd <- function(object, ...) {
  val <- -object$results$value
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
