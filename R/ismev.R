# ================================== ismev ================================== #

#' Loglikelihood adjustment for ismev fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[ismev]{gev.fit}}, \code{\link[ismev]{gpd.fit}},
#' \code{\link[ismev]{pp.fit}} and \code{\link[ismev]{rlarg.fit}} in the
#' \code{\link[ismev]{ismev}} package.  If regression modelling is used then
#' the model will need to be re-fitted, see \code{\link{ismev_refits}}.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#'
#' If regression modelling is used then the ismev functions
#' \code{\link[ismev]{gev.fit}}, \code{\link[ismev]{gpd.fit}},
#' \code{\link[ismev]{pp.fit}} and \code{\link[ismev]{rlarg.fit}}
#' return residuals but \code{\link{alogLik}} needs the raw data.
#' The model will need to be re-fitted, using one of the functions in
#' \code{\link{ismev_refits}}, and the user will be prompted to do this
#' by an error message produced by \code{\link{alogLik}}.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "ismev")}.
#'   The remaining 2 components depend on the model that was fitted.
#'   The 4th component is:
#'   \code{"gev"} if \code{\link[ismev]{gev.fit}}
#'   (or \code{\link{gev_refit}}) was used;
#'   \code{"gpd"} if \code{\link[ismev]{gpd.fit}}
#'   (or \code{\link{gpd_refit}}) was used;
#'   \code{"pp"} \code{\link[ismev]{pp.fit}}
#'   (or \code{\link{pp_refit}}) was used;
#'   \code{"rlarg"} \code{\link[ismev]{rlarg.fit}}
#'   (or \code{\link{rlarg_refit}}) was used.
#'   The 5th component is
#'   \code{"stat"} if \code{x$trans = FALSE} and
#'   \code{"nonstat"} if \code{x$trans = TRUE}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the ismev package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
#'
#'   # GEV model -----
#'
#'   # An example from the ismev::gev.fit documentation
#'   gev_fit <- gev.fit(revdbayes::portpirie, show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit)
#'   summary(adj_gev_fit)
#'
#'   # An example from chapter 6 of Coles (2001)
#'   data(fremantle)
#'   xdat <- fremantle[, "SeaLevel"]
#'   # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
#'   ydat <- cbind(fremantle[, "Year"] - 1896, fremantle[, "SOI"])
#'   gev_fit <- gev_refit(xdat, ydat, mul = 1:2, show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit)
#'   summary(adj_gev_fit)
#'
#'   # An example from Chandler and Bate (2007)
#'   gev_fit <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4,
#'                        show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit, cluster = ow$year)
#'   summary(adj_gev_fit)
#'   # Get closer to the values reported in Table 2 of Chandler and Bate (2007)
#'   gev_fit <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4,
#'                        show = FALSE, method = "BFGS")
#'   # Call sandwich::meatCL() with cadjust = FALSE
#'   adj_gev_fit <- alogLik(gev_fit, cluster = ow$year, cadjust = FALSE)
#'   summary(adj_gev_fit)
#'
#'   # GP model -----
#'
#'   # An example from the ismev::gpd.fit documentation
#'   \donttest{
#'   data(rain)
#'   rain_fit <- gpd.fit(rain, 10, show = FALSE)
#'   adj_rain_fit <- alogLik(rain_fit)
#'   summary(adj_rain_fit)
#'   # Continuing to the regression example on page 119 of Coles (2001)
#'   ydat <- as.matrix((1:length(rain)) / length(rain))
#'   reg_rain_fit <- gpd_refit(rain, 30, ydat = ydat, sigl = 1, siglink = exp,
#'                             show = FALSE)
#'   adj_reg_rain_fit <- alogLik(reg_rain_fit)
#'   summary(adj_reg_rain_fit)
#'   }
#'   # Binomial-GP model -----
#'
#'   # Use Newlyn seas surges data from the exdex package
#'   surges <- exdex::newlyn
#'   u <- quantile(surges, probs = 0.9)
#'   newlyn_fit <- gpd.fit(surges, u, show = FALSE)
#'   # Create 5 clusters each corresponding approximately to 1 year of data
#'   cluster <- rep(1:5, each = 579)[-1]
#'   adj_newlyn_fit <- alogLik(newlyn_fit, cluster = cluster, binom = TRUE,
#'                             cadjust = FALSE)
#'   summary(adj_newlyn_fit)
#'   summary(attr(adj_newlyn_fit, "pu_aloglik"))
#'
#'   # Add inference about the extremal index theta, using K = 1
#'   adj_newlyn_theta <- alogLik(newlyn_fit, cluster = cluster, binom = TRUE,
#'                               k = 1, cadjust = FALSE)
#'   summary(attr(adj_newlyn_theta, "theta"))
#'
#'   # PP model -----
#'
#'   # An example from the ismev::pp.fit documentation
#'   data(rain)
#'   # Start from the mle to save time
#'   init <- c(40.55755732, 8.99195409, 0.05088103)
#'   muinit <- init[1]
#'   siginit <- init[2]
#'   shinit <- init[3]
#'   rain_fit <- pp_refit(rain, 10, muinit = muinit, siginit = siginit,
#'                        shinit = shinit, show = FALSE)
#'   adj_rain_fit <- alogLik(rain_fit)
#'   summary(adj_rain_fit)
#'
#'   # An example from chapter 7 of Coles (2001).
#'   # Code from demo ismev::wooster.temps
#'   data(wooster)
#'   x <- seq(along = wooster)
#'   usin <- function(x, a, b, d) {
#'     return(a + b * sin(((x - d) * 2 * pi) / 365.25))
#'   }
#'   wu <- usin(x, -30, 25, -75)
#'   ydat <- cbind(sin(2 * pi * x / 365.25), cos(2 * pi *x / 365.25))
#'   # Start from the mle to save time
#'   init <- c(-15.3454188, 9.6001844, 28.5493828, 0.5067104, 0.1023488,
#'             0.5129783, -0.3504231)
#'   muinit <- init[1:3]
#'   siginit <- init[4:6]
#'   shinit <- init[7]
#'   wooster.pp <- pp_refit(-wooster, threshold = wu, ydat = ydat, mul = 1:2,
#'                          sigl = 1:2, siglink = exp, method = "BFGS",
#'                          muinit = muinit, siginit = siginit, shinit = shinit,
#'                          show = FALSE)
#'   adj_pp_fit <- alogLik(wooster.pp)
#'   summary(adj_pp_fit)
#'
#'   # r-largest order statistics model -----
#'
#'   # An example based on the ismev::rlarg.fit() documentation
#'   vdata <- revdbayes::venice
#'   rfit <- rlarg.fit(vdata, muinit = 120.54, siginit = 12.78,
#'                     shinit = -0.1129, show = FALSE)
#'   adj_rfit <- alogLik(rfit)
#'   summary(adj_rfit)
#'
#'   \donttest{
#'   # Adapt this example to add a covariate
#'   set.seed(30102019)
#'   ydat <- matrix(runif(nrow(vdata)), nrow(vdata), 1)
#'   rfit2 <- rlarg_refit(vdata, ydat = ydat, mul = 1,
#'                        muinit = c(120.54, 0), siginit = 12.78,
#'                        shinit = -0.1129, show = FALSE)
#'   adj_rfit2 <- alogLik(rfit2)
#'   summary(adj_rfit2)
#'   }
#' }
#' @name ismev
NULL
## NULL

#' @rdname ismev
#' @export
alogLik.gev.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_lax <- list(ismev_gev = "gev.fit")
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
  if (x$trans) {
    class(res) <- c("lax", "chandwich", "ismev", "gev", "nonstat")
  } else {
    class(res) <- c("lax", "chandwich", "ismev", "gev", "stat")
  }
  return(res)
}

#' @rdname ismev
#' @export
alogLik.pp.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_lax <- list(ismev_pp = "pp.fit")
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
  if (x$trans) {
    class(res) <- c("lax", "chandwich", "ismev", "pp", "nonstat")
  } else {
    class(res) <- c("lax", "chandwich", "ismev", "pp", "stat")
  }
  return(res)
}

#' @rdname ismev
#' @export
alogLik.gpd.fit <- function(x, cluster = NULL, use_vcov = TRUE, binom = FALSE,
                            k, inc_cens = TRUE, ...) {
  # List of ismev objects supported
  supported_by_lax <- list(ismev_gpd = "gpd.fit")
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
  # If cluster has been supplied then we need to check it's length is correct.
  # If binom = TRUE then it must have the same length as the raw data x$xdata.
  # If binom = FALSE then, alternatively, it may have the same length as the
  # vector of threshold excesses.
  # We check this and, if necessary, extract the cluster values that correspond
  # to the threshold excesses.
  # Save the full cluster vector for use later for the Bernoulli distribution
  full_cluster <- cluster
  if (!is.null(cluster)) {
    clength <- length(cluster)
    if (binom) {
      if (clength != x$n) {
        stop("binom=TRUE: ''cluster'' must be the same length as the raw data")
      }
      cluster <- cluster[x$xdata > x$threshold]
    } else {
      if (clength == x$n) {
        cluster <- cluster[x$xdata > x$threshold]
      } else if (clength == x$nexc) {
      } else {
        stop("''cluster'' must have the same length as either the raw or
             exceedance data")
      }
    }
  }
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (x$trans) {
    class(res) <- c("lax", "chandwich", "ismev", "gpd", "nonstat")
  } else {
    class(res) <- c("lax", "chandwich", "ismev", "gpd", "stat")
  }
  # For a stationary model, add the adjusted binomial log-likelihood for the
  # exceedance probability p_u if binom = TRUE
  if (!x$trans && binom) {
    exc <- x$xdata > x$threshold
    fitb <- fit_bernoulli(exc)
    afitb <- alogLik(fitb, cluster = full_cluster, ...)
    attr(res, "pu_aloglik") <- afitb
    class(res) <- c("lax", "chandwich", "ismev", "bin-gpd", "stat")
  }
  # For a stationary model, and if a valid k has been supplied, then calculate
  # the K-gaps model log-likelihood
  if (!x$trans && !missing(k) && length(k) == 1) {
    if (!is.numeric(k) || k < 0 || length(k) != 1) {
      stop("k must be a non-negative scalar")
    }
    # Call exdex::kgaps() to estimate the extremal index theta and store values
    # from which the log-likelihood can be calculated
    theta <- exdex::kgaps(data = x$xdata, u = x$threshold, k = k,
                          inc_cens = inc_cens)
    attr(res, "theta") <- theta
  }
  return(res)
}

#' @rdname ismev
#' @export
alogLik.rlarg.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_lax <- list(ismev_rlarg = "rlarg.fit")
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
  if (x$trans) {
    class(res) <- c("lax", "chandwich", "ismev", "rlarg", "nonstat")
  } else {
    class(res) <- c("lax", "chandwich", "ismev", "rlarg", "stat")
  }
  return(res)
}
