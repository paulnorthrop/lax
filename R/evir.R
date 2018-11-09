# =================================== evir ================================== #

#' Loglikelihood adjustment of evir fits
#'
#' Description
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the evir package
#' got_evir <- requireNamespace("evir", quietly = TRUE)
#' if (got_evir) {
#'   library(evir)
#'   # An example from the evir::gev documentation
#'   data(bmw)
#'   out <- evir::gev(bmw, "month")
#' }
#' @name evir
NULL
## NULL

#' @rdname evir
#' @export
alogLik.gev <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # For an "evd" object then reverse the class list to use the "evd" method
  if (inherits(x, "evd")) {
    class(x) <- rev(class(x))
    return(alogLik(x, cluster = cluster, use_vcov = use_vcov, ...))
  }
  # List of evir or fExtremes objects supported
  # We need to
  supported_by_oolax <- list(evir_gev = c("gev"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_oolax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_oolax)[which(is_supported)]
  class(x) <- name_of_class
  # Call oola::adjust_object to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("oolax", "chandwich", "evir", "gev", "stat")
  return(res)
}

