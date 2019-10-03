# ================================== mev ================================== #

#' Loglikelihood adjustment for mev fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[mev]{fit.gev}}, \code{\link[mev]{fit.gpd}}, and
#' \code{\link[mev]{fit.pp}} and \code{\link[mev]{fit.rlarg}} in the
#' mev package.
#'
#' @inherit alogLik params details references
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "mev")}.
#'   The remaining 2 components depend on the model that was fitted.
#'   The 4th component is:
#'   \code{"gev"} if \code{\link[mev]{fit.gev}} was used;
#'   \code{"gpd"} if \code{\link[mev]{fit.gpd}} was used;
#'   \code{"pp"} \code{\link[mev]{fit.pp}} was used;
#'   The 5th component is \code{"stat"}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the mev package
#' got_mev <- requireNamespace("mev", quietly = TRUE)
#'
#' if (got_mev) {
#'   library(mev)
#'   # An example from the ismev::gev.fit documentation
#'   gev_mev <- fit.gev(revdbayes::portpirie, show = FALSE)
#'   adj_gev_mev <- alogLik(gev_mev)
#'   summary(adj_gev_mev)
#' }
#' @name mev
NULL
## NULL

#' @rdname mev
#' @export
alogLik.mev_gev <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of mev objects supported
  supported_by_lax <- list(mev_gev = "mev_gev")
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
  class(res) <- c("lax", "chandwich", "mev", "gev", "stat")
  return(res)
}
