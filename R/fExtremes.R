# ================================ fExtremes ================================ #

#' Loglikelihood adjustment for fExtremes fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[fExtremes:GevModelling]{gevFit}},
#' \code{\link[fExtremes:GevModelling]{gumbelFit}} and
#' \code{\link[fExtremes:GpdModelling]{gpdFit}}
#' in the \code{\link[fExtremes:00Extremes-package]{fExtremes}} package.
#' The model must have been fitted using maximum likelihood estimation.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "fExtremes")}.
#'   The remaining 2 components depend on the model that was fitted.
#'   If \code{\link[fExtremes:GevModelling]{gevFit}} or
#'   \code{\link[fExtremes:GevModelling]{gumbelFit}} was used then these
#'   components are \code{c("gev", "stat")}.
#'   If \code{\link[fExtremes:GpdModelling]{gpdFit}} was used then these
#'   components are \code{c("gpd", "stat")}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the fExtremes package
#' got_fExtremes <- requireNamespace("fExtremes", quietly = TRUE)
#' if (got_fExtremes) {
#'   library(fExtremes)
#'
#'   # GEV
#'   # An example from the fExtremes::gevFit documentation
# '  # Simulate GEV Data
#'   set.seed(4082019)
#'   x <- fExtremes::gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
#'   # Fit GEV distribution by maximum likelihood estimation
#'   fit <- fExtremes::gevFit(x)
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#'
#'   # GP
#'   # An example from the fExtremes::gpdFit documentation
#'   # Simulate GP data
#'   x <- fExtremes::gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
#'   # Fit GP distribution by maximum likelihood estimation
#'   fit <- fExtremes::gpdFit(x, u = min(x))
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#' }
#' @name fExtremes
NULL
## NULL

#' @rdname fExtremes
#' @export
alogLik.fGEVFIT <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x@method[2] != "mle") {
    stop("Loglikelihood adjustment is only relevant when type = ''mle''")
  }
  # List of fExtremes objects supported
  supported_by_lax <- list(fExtremes_gev = c("fGEVFIT"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_lax)) {
    is_supported[i] <- identical(class(x)[1], unlist(supported_by_lax[i],
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
  class(res) <- c("lax", "chandwich", "fExtremes", "gev", "stat")
  return(res)
}

#' @rdname fExtremes
#' @export
alogLik.fGPDFIT <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x@method[2] != "mle") {
    stop("Loglikelihood adjustment is only relevant when type = ''mle''")
  }
  # List of fExtremes objects supported
  supported_by_lax <- list(fExtremes_gpd = c("fGPDFIT"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_lax)) {
    is_supported[i] <- identical(class(x)[1], unlist(supported_by_lax[i],
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
  class(res) <- c("lax", "chandwich", "fExtremes", "gpd", "stat")
  return(res)
}
