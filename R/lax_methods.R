# nobs method for objects of class "lax"
# We add this because "chandwich" doesn't have a nobs method

#' @export
nobs.lax <- function(object, ...) {
  return(attr(object, "nobs"))
}
