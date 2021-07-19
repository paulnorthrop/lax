# ================================== eva ================================== #

#' Loglikelihood adjustment for eva fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[eva]{gevrFit}} and \code{\link[eva]{gpdFit}} in the eva package.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#'
#'   In the stationary case (no covariates) the function
#'   \code{\link[eva]{gevrFit}} and \code{\link[eva]{gpdFit}} in the eva
#'   package offer standard errors based on the expected information or on the
#'   observed information, via the argument \code{information}.  In contrast,
#'   \code{alogLik()} always bases calculations on the observed information
#'   matrix. Therefore, unadjusted standard errors resulting from
#'   \code{alogLik()} may be different the corresponding standard errors
#'   from  \code{\link[eva]{gevrFit}} or \code{\link[eva]{gpdFit}}.
#'
#'   For \code{\link[eva]{gevrFit}} only GEV fits (\code{gumbel = FALSE}) are
#'   supported.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "eva")}.
#'   The 4th component depends on which model was fitted.
#'   \code{"rlarg"} if \code{\link[eva]{gevrFit}} was used;
#'   \code{"gpd"} if \code{\link[eva]{gpdFit}} was used.
#'   The 5th component is
#'   \code{"stat"} if there are no covariates in the mode and
#'   \code{"nonstat"} otherwise.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the eva package
#' got_eva <- requireNamespace("eva", quietly = TRUE)
#'
#' if (got_eva) {
#'   library(eva)
#'   # An example from the eva::gpdFit documentation
#'   set.seed(7)
#'   x <- eva::rgpd(2000, loc = 0, scale = 2, shape = 0.2)
#'   mle_fit <- eva::gpdFit(x, threshold = 4, method = "mle")
#'   adj_mle_fit <- alogLik(mle_fit)
#'   summary(adj_mle_fit)
#'
#'   # Another example from the eva::gpdFit documentation
#'   # A linear trend in the scale parameter
#'   set.seed(7)
#'   n <- 300
#'   x2 <- eva::rgpd(n, loc = 0, scale = 1 + 1:n / 200, shape = 0)
#'   covs <- as.data.frame(seq(1, n, 1))
#'   names(covs) <- c("Trend1")
#'   result1 <- eva::gpdFit(x2, threshold = 0, scalevars = covs,
#'                          scaleform = ~ Trend1)
#'   adj_result1 <- alogLik(result1)
#'   summary(adj_result1)
#'
#'   # An example from the eva::gevrFit documentation
#'   set.seed(7)
#'   x1 <- eva::rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
#'   result1 <- eva::gevrFit(x1, method = "mle")
#'   adj_result1 <- alogLik(result1)
#'   summary(adj_result1)
#'
#'   # Another example from the eva::gevrFit documentation
#'   # A linear trend in the location and scale parameter
#'   n <- 100
#'   r <- 10
#'   x2 <- eva::rgevr(n, r, loc = 100 + 1:n / 50,  scale = 1 + 1:n / 300,
#'                    shape = 0)
#'   covs <- as.data.frame(seq(1, n, 1))
#'   names(covs) <- c("Trend1")
#'   # Create some unrelated covariates
#'   covs$Trend2 <- rnorm(n)
#'   covs$Trend3 <- 30 * runif(n)
#'   result2 <- eva::gevrFit(data = x2, method = "mle", locvars = covs,
#'                           locform = ~ Trend1 + Trend2*Trend3,
#'                           scalevars = covs, scaleform = ~ Trend1)
#'   adj_result2 <- alogLik(result2)
#'   summary(adj_result2)
#' }
#' @name eva
NULL
## NULL

#' @rdname eva
#' @export
alogLik.gevrFit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x$method != "mle") {
    stop("The model must have been fitted using maximum likelihood estimation")
  }
  if (x$gumbel) {
    stop("Only GEV fits (gumbel = FALSE) are supported.")
  }
  # List of eva objects supported
  supported_by_lax <- list(laxeva_rlarg = "gevrFit")
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
  # If this is a stationary model and the expected information was used then
  # do not use the vcov method, because it will return the expected information
  if (x$stationary && x$information == "expected") {
    use_vcov <- FALSE
  }
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (x$stationary) {
    class(res) <- c("lax", "chandwich", "eva", "rlarg", "stat")
  } else {
    class(res) <- c("lax", "chandwich", "eva", "rlarg", "nonstat")
  }
  return(res)
}

#' @rdname eva
#' @export
alogLik.gpdFit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x$method != "mle") {
    stop("The model must have been fitted using maximum likelihood estimation")
  }
  # List of eva objects supported
  supported_by_lax <- list(laxeva_gpd = "gpdFit")
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
  # If this is a stationary model and the expected information was used then
  # do not use the vcov method, because it will return the expected information
  if (x$stationary && x$information == "expected") {
    use_vcov <- FALSE
  }
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (x$stationary) {
    class(res) <- c("lax", "chandwich", "eva", "gpd", "stat")
  } else {
    class(res) <- c("lax", "chandwich", "eva", "gpd", "nonstat")
  }
  return(res)
}
