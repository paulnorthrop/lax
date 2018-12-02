# ================================ fExtremes ================================ #

#' Loglikelihood adjustment of fExtremes fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment of fitted
#' extreme value model objects produced by the
#' \code{\link[fExtremes:00Extremes-package]{fExtremes}} package.
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the fExtremes package
#' got_fExtremes <- requireNamespace("fExtremes", quietly = TRUE)
#' if (got_fExtremes) {
#'   library(fExtremes)
#'
#'   # GEV
#'   # An example from the fExtremes::gevFit documentation
# '  # Simulate GEV Data
#'   x <- fExtremes::gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
#'   # Fit GEV distribution by maximum likelihood estimation
#'   fit <- gevFit(x)
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#'
#'   # GP
#'   # An example from the fExtremes::gpdFit documentation
#'   # Simulate GP data
#'   x <- gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
#'   # Fit GP distribution by maximum likelihood estimation
#'   fit <- gpdFit(x, u = min(x))
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
  supported_by_oolax <- list(fExtremes_gev = c("fGEVFIT"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x)[1], unlist(supported_by_oolax[i],
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
  class(res) <- c("oolax", "chandwich", "fExtremes", "gev", "stat")
  return(res)
}

#' @rdname fExtremes
#' @export
alogLik.fGPDFIT <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x@method[2] != "mle") {
    stop("Loglikelihood adjustment is only relevant when type = ''mle''")
  }
  # List of fExtremes objects supported
  supported_by_oolax <- list(fExtremes_gpd = c("fGPDFIT"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x)[1], unlist(supported_by_oolax[i],
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
  class(res) <- c("oolax", "chandwich", "fExtremes", "gpd", "stat")
  return(res)
}
