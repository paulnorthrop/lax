# ================================== ismev ================================== #

#' Loglikelihood adjustment of ismev fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment of fitted
#' extreme value model objects produced by the \code{\link[ismev]{ismev}}
#' package.
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the evd package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
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
#'   gev_fit <- oogev.fit(xdat, ydat, mul = 1:2, show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit)
#'   summary(adj_gev_fit)
#'
#'   # An example from the ismev::gpd.fit documentation
#'   data(rain)
#'   rain_fit <- gpd.fit(rain, 10, show = FALSE)
#'   adj_rain_fit <- alogLik(rain_fit)
#'   summary(adj_rain_fit)
#'   # Continuing to the regression example on page 119 of Coles (2001)
#'   ydat <- as.matrix((1:length(rain)) / length(rain))
#'   reg_rain_fit <- oogpd.fit(rain, 30, ydat = ydat, sigl = 1, siglink = exp,
#'                             show = FALSE)
#'   adj_reg_rain_fit <- alogLik(reg_rain_fit)
#'   summary(adj_reg_rain_fit)
#'
#'   # An example from the ismev::pp.fit documentation
#'   data(rain)
#'   # Start from the mle to save time
#'   init <- c(40.55755732, 8.99195409, 0.05088103)
#'   muinit <- init[1]
#'   siginit <- init[2]
#'   shinit <- init[3]
#'   rain_fit <- oopp.fit(rain, 10, muinit = muinit, siginit = siginit,
#'                        shinit = shinit, show = FALSE)
#'   adj_rain_fit <- alogLik(rain_fit)
#'   summary(adj_rain_fit)
#'
#'   # An example from chapter 7 of Coles (2001).
#'   # Code from demo ismev::wooster.temps
#'   data(wooster)
#'   x <- seq(along = wooster)
#'   usin <- function(x, a, b, d) {
#'    a + b * sin(((x - d) * 2 * pi) / 365.25)
#'   }
#'   wu <- usin(x, -30, 25, -75)
#'   ydat <- cbind(sin(2 * pi * x / 365.25), cos(2 * pi *x / 365.25))
#'   # Start from the mle to save time
#'   init <- c(-15.3454188, 9.6001844, 28.5493828, 0.5067104, 0.1023488,
#'             0.5129783, -0.3504231)
#'   muinit <- init[1:3]
#'   siginit <- init[4:6]
#'   shinit <- init[7]
#'   wooster.pp <- oopp.fit(-wooster, threshold = wu, ydat = ydat, mul = 1:2,
#'                          sigl = 1:2, siglink = exp, method = "BFGS",
#'                          muinit = muinit, siginit = siginit, shinit = shinit,
#'                          show = FALSE)
#'   adj_pp_fit <- alogLik(wooster.pp)
#'   summary(adj_pp_fit)
#'
#'   # An example from Chandler and Bate (2007)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- as.matrix(rep(c(1, -1), each = length(y) / 2))
#'   gev_fit <- oogev.fit(y, x, mul = 1, sigl = 1, shl = 1, show = FALSE)
#'   year <- rep(rownames(chandwich::owtemps), 2)
#'   adj_gev_fit <- alogLik(gev_fit, cluster = year)
#'   summary(adj_gev_fit)
#'   # Get closer to the values reported in Table 2 of Chandler and Bate (2007)
#'   gev_fit <- oogev.fit(y, x, mul = 1, sigl = 1, shl = 1, show = FALSE,
#'                        method = "BFGS")
#'   year <- rep(rownames(chandwich::owtemps), 2)
#'   # Call sandwich::meatCL() with cadjust = FALSE
#'   adj_gev_fit <- alogLik(gev_fit, cluster = year, cadjust = FALSE)
#'   summary(adj_gev_fit)
#' }
#' @name ismev
NULL
## NULL

#' @rdname ismev
#' @export
alogLik.gev.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_oolax <- list(ismev_gev = "gev.fit")
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
  if (x$trans) {
    class(res) <- c("oolax", "chandwich", "ismev", "gev", "nonstat")
  } else {
    class(res) <- c("oolax", "chandwich", "ismev", "gev", "stat")
  }
  return(res)
}

#' @rdname ismev
#' @export
alogLik.pp.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_oolax <- list(ismev_pp = "pp.fit")
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
  if (x$trans) {
    class(res) <- c("oolax", "chandwich", "ismev", "pp", "nonstat")
  } else {
    class(res) <- c("oolax", "chandwich", "ismev", "pp", "stat")
  }
  return(res)
}

#' @rdname ismev
#' @export
alogLik.gpd.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of ismev objects supported
  supported_by_oolax <- list(ismev_gpd = "gpd.fit")
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
  if (x$trans) {
    class(res) <- c("oolax", "chandwich", "ismev", "gpd", "nonstat")
  } else {
    class(res) <- c("oolax", "chandwich", "ismev", "gpd", "stat")
  }
  return(res)
}
