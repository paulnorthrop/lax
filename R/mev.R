# ================================== mev ================================== #

#' Loglikelihood adjustment for mev fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[mev]{fit.gev}}, \code{\link[mev]{fit.gpd}}, and
#' \code{\link[mev]{fit.pp}} and \code{\link[mev]{fit.rlarg}} in the
#' mev package.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#'
#' If \code{x} was returned from \code{\link[mev]{fit.pp}} then the data
#' \code{xdat} supplied to \code{\link[mev]{fit.pp}} must contain \emph{all}
#' the data, both threshold exceedances and non-exceedances.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "mev")}.
#'   The 4th component depends on which model was fitted.
#'   \code{"gev"} if \code{\link[mev]{fit.gev}} was used;
#'   \code{"gpd"} if \code{\link[mev]{fit.gpd}} was used;
#'   \code{"pp"} \code{\link[mev]{fit.pp}} was used;
#'   \code{"egp"} \code{\link[mev]{fit.egp}} was used;
#'   \code{"rlarg"} \code{\link[mev]{fit.rlarg}} was used;
#'   The 5th component is \code{"stat"} (for stationary).
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the mev package
#' got_mev <- requireNamespace("mev", quietly = TRUE)
#'
#' if (got_mev) {
#'   library(mev)
#'   # An example from the mev::gev.fit documentation
#'   gev_mev <- fit.gev(revdbayes::portpirie)
#'   adj_gev_mev <- alogLik(gev_mev)
#'   summary(adj_gev_mev)
#'
#'   # Use simulated data
#'   set.seed(1112019)
#'   x <- revdbayes::rgp(365 * 10, loc = 0, scale = 1, shape = 0.1)
#'   pfit <- fit.pp(x, threshold = 1, npp = 365)
#'   # (To do: delete the next two lines after new mev hits CRAN)
#'   pfit$xdat <- x
#'   pfit$npp <- 365
#'   adj_pfit <- alogLik(pfit)
#'   summary(adj_pfit)
#'
#'   # An example from the mev::fit.gpd documentation
#'   gpd_mev <- fit.gpd(eskrain, threshold = 35, method = 'Grimshaw')
#'   adj_gpd_mev <- alogLik(gpd_mev)
#'   summary(adj_gpd_mev)
#'
#'   # An example from the mev::fit.egp documentation
#'   # (model = "egp1" and model = "egp3" also work)
#'   xdat <- evd::rgpd(n = 100, loc = 0, scale = 1, shape = 0.5)
#'   fitted <- fit.egp(xdat = xdat, thresh = 1, model = "egp2", show = FALSE)
#'   adj_fitted <- alogLik(fitted)
#'   summary(adj_fitted)
#'
#'   # An example from the mev::fit.rlarg documentation
#'   set.seed(31102019)
#'   xdat <- rrlarg(n = 10, loc = 0, scale = 1, shape = 0.1, r = 4)
#'   fitr <- fit.rlarg(xdat)
#'   adj_fitr <- alogLik(fitr)
#'   summary(adj_fitr)
#' }
#' @name mev
NULL
## NULL

#' @rdname mev
#' @export
alogLik.mev_gev <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of mev objects supported
  supported_by_lax <- list(laxmev_gev = "mev_gev")
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

#' @rdname mev
#' @export
alogLik.mev_pp <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (all(x$xdat > x$threshold)) {
    stop("xdat must contain all the data: exceedances and non-exceedances.")
  }
  # List of mev objects supported
  supported_by_lax <- list(laxmev_pp = "mev_pp")
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
  class(res) <- c("lax", "chandwich", "mev", "pp", "stat")
  return(res)
}

#' @rdname mev
#' @export
alogLik.mev_gpd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of mev objects supported
  supported_by_lax <- list(laxmev_gpd = "mev_gpd")
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
  class(res) <- c("lax", "chandwich", "mev", "gpd", "stat")
  return(res)
}

#' @rdname mev
#' @export
alogLik.mev_egp <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of mev objects supported
  supported_by_lax <- list(laxmev_egp = "mev_egp")
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
  class(res) <- c("lax", "chandwich", "mev", "egp", "stat")
  return(res)
}

#' @rdname mev
#' @export
alogLik.mev_rlarg <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of mev objects supported
  supported_by_lax <- list(laxmev_rlarg = c("mev_rlarg", "mev_gev"))
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
  class(res) <- c("lax", "chandwich", "mev", "rlarg", "stat")
  return(res)
}
