#' @export
estfun.default <- function(x, ...) {
  my_fn <- paste("logLikFn.", class(x), sep = "")
  U <- numDeriv::jacobian(eval(as.name(my_fn)), x = coef(x), ...)
  return(U)
}
