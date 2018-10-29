# ================================== ismev ================================== #

#' Loglikelihood adjustment of ismev fits
#'
#' Description
#'
#' @param x A fitted model object.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   This must have the same length as the vector returned by \code{loglik}.
#'   If \code{cluster} is not supplied then it is set inside
#'   \code{\link[chandwich]{adjust_loglik}} under the assumption that
#'   each observation forms its own cluster.
#'
#'   If the sandwich package
#'   \href{http://dx.doi.org/10.18637/jss.v016.i09}{(Zeleis, 2006)}
#'   is used to estimate the quantities required to adjust the loglikelihood
#'   then \code{cluster} determines whether the variance matrix \code{V}
#'   of the score vector is estimated using \code{\link[sandwich]{meat}}
#'   (\code{cluster} is \code{NULL}) or
#'   \code{\link[sandwich:vcovCL]{meatCL}}
#'   (\code{cluster} is not \code{NULL}).
#' @param use_vcov A logical scalar.  Should we use the \code{vcov} S3 method
#'   for \code{x} (if this exists) to estimate the Hessian of the independence
#'   loglikelihood to be passed as the argument \code{H} to
#'   \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \code{H} is estimated inside
#'   \code{\link[chandwich]{adjust_loglik}} using
#'   \code{\link[stats:optim]{optimHess}}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}}, if \code{cluster = NULL},
#'   or \code{\link[sandwich:vcovCL]{meatCL}}, otherwise.
#' @examples
#' # We need the evd package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
#'   # An example from the ismev::gev.fit documentation
#'   data(portpirie)
#'   gev_fit <- gev.fit(portpirie[,2])
#'   adj_gev_fit <- alogLik(gev_fit)
#'   summary(adj_gev_fit)
#' }
#' @export
alogLik.gev.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evd objects supported
  supported_by_oolax <- list(ismev_gev = c("gev.fit"))
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
  class(res) <- c("oolax", "chandwich")
  return(res)
}
