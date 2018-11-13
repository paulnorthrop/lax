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
  return(solve(object$results$hessian))
}

#' @export
logLik.fevd <- function(object, ...) {
  val <- -object$results$value
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
