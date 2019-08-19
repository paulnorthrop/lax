# =================================== evir ================================== #

#' Loglikelihood adjustment of evir fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment of fitted
#' extreme value model objects produced by the evir package.
#' If \code{\link[evir]{pot}} was used to produce \code{x} then the model will
#' need to be re-fitted using \code{\link{pot_refit}}.
#'
#' @inherit alogLik params details return references seealso
#' @examples
#' # We need the evir package
#' got_evir <- requireNamespace("evir", quietly = TRUE)
#' if (got_evir) {
#'   library(evir)
#'   # An example from the evir::gev documentation
#'   data(bmw)
#'   out <- evir::gev(bmw, "month")
#'   adj_out <- alogLik(out)
#'   summary(adj_out)
#'
#'   # An example from the evir::gpd documentation
#'   data(danish)
#'   out <- evir::gpd(danish, 10)
#'   adj_out <- alogLik(out)
#'   summary(adj_out)
#'
#'   # An example from the evir::pot documentation
#'   # We use lax::pot_refit() to return the input data
#'   out <- pot_refit(danish, 10)
#'   adj_out <- alogLik(out)
#'   summary(adj_out)
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
  # List of evir objects supported
  supported_by_lax <- list(evir_gev = c("gev"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_lax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_lax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_lax)[which(is_supported)]
  class(x) <- name_of_class
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("lax", "chandwich", "evir", "gev", "stat")
  return(res)
}

#' @rdname evir
#' @export
alogLik.gpd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x$method != "ml") {
    stop("Loglikelihood adjustment is only relevant when method = ''ml''")
  }
  # List of evir objects supported
  supported_by_lax <- list(evir_gpd = c("gpd"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_lax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_lax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_lax)[which(is_supported)]
  class(x) <- name_of_class
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("lax", "chandwich", "evir", "gpd", "stat")
  return(res)
}

#' @rdname evir
#' @export
alogLik.potd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evir objects supported
  supported_by_lax <- list(evir_pot = c("potd"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_lax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_lax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_lax)[which(is_supported)]
  class(x) <- name_of_class
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("lax", "chandwich", "evir", "potd", "stat")
  return(res)
}
