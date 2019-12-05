# =================================== POT =================================== #

#' Loglikelihood adjustment for POT fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from
#' \code{\link[POT]{fitGPD}} function in the POT package.
#' The model must have been fitted using maximum likelihood estimation.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'
#'   \code{class(x)} is \code{c("lax", "chandwich", "POT", "pot", "gpd")}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the POT package
#' got_POT <- requireNamespace("POT", quietly = TRUE)
#'
#' if (got_POT) {
#'   library(POT)
#'   # An example from the POT::fitgpd documentation.
#'   set.seed(4082019)
#'   x <- POT::rgpd(200, 1, 2, 0.25)
#'   fit <- fitgpd(x, 1, "mle")
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
