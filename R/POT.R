# =================================== POT =================================== #

#' Loglikelihood adjustment of POT fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment of fitted
#' extreme value model objects produced by the POT package.
#'
#' @inherit alogLik params details return references seealso
#' @examples
#' # We need the POT package
#' got_POT <- requireNamespace("POT", quietly = TRUE)
#'
#' if (got_POT) {
#'   library(POT)
#'   # An example from the POT::fitgpd documentation.
#'   set.seed(4082019)
#'   x <- POT::rgpd(200, 1, 2, 0.25)
#'   fit <- POT::fitgpd(x, 1, "mle")
#'   adj_fit <- alogLik(fit)
#' }
#' @name POT
NULL
## NULL

#' @rdname POT
#' @export
alogLik.uvpot <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x$est != "MLE") {
    stop("Loglikelihood adjustment is only relevant when est = ''mle''")
  }
  # List of evd objects supported
  supported_by_lax <- list(POT_fitgpd = c("uvpot", "pot"))
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
  class(res) <- c("lax", "chandwich", "POT", "pot", "gpd")
  return(res)
}
