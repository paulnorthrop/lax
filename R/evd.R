# =================================== evd =================================== #

#' Loglikelihood adjustment for evd fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the functions
#' \code{\link[evd]{fgev}} and \code{\link[evd]{fpot}} in the evd package.
#' If \code{x} is returned from \code{\link[evd]{fgev}} then the call must
#' have used \code{prob = NULL}.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "evd")}.
#'   The remaining 2 components depend on the model that was fitted.
#'   If \code{\link[evd]{fgev}} was used then these components are
#'   \code{c("gev", "stat")} if \code{nsloc} was \code{NULL} and
#'   \code{c("gev", "nonstat")} if \code{nsloc} was not \code{NULL}.
#'   If \code{\link[evd]{fpot}} was used then these components are
#'   \code{c("pot", "gpd")} if \code{model} was \code{"gpd"} and
#'   \code{c("pot", "pp")} if \code{model} was \code{"pp"}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the evd package
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   set.seed(3082019)
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
#'   adj_fgev <- alogLik(M1)
#'   summary(adj_fgev)
#'
#'   # An example from Chandler and Bate (2007)
#'   owfit <- fgev(ow$temp, nsloc = ow$loc)
#'   adj_owfit <- alogLik(owfit, cluster = ow$year)
#'   summary(adj_owfit)
#'
#'   # An example from the evd::fpot documentation
#'   set.seed(3082019)
#'   uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
#'   M1 <- fpot(uvdata, 1)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#'   # Fit using the pp model, rather than the gpd
#'   M1 <- fpot(uvdata, 1, model = "pp", npp = 365)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#' }
#' @name evd
NULL
## NULL

#' @rdname evd
#' @export
alogLik.evd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # Reverse the order of the classes: if they were reversed in alogLik.gev()
  if (class(x)[1] == "evd") {
    class(x) <- rev(class(x))
  }
  # List of evd objects supported
  supported_by_lax <- list(evd_fgev = c("gev", "uvevd", "evd"),
                           evd_fpot = c("pot", "uvevd", "evd"))
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
  if (name_of_class == "evd_fgev" && !is.null(x$prob)) {
    stop("x$prob must be NULL")
  }
  class(x) <- name_of_class
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (name_of_class == "evd_fgev") {
    if (is.null(x$nsloc)) {
      class(res) <- c("lax", "chandwich", "evd", "gev", "stat")
    } else {
      class(res) <- c("lax", "chandwich", "evd", "gev", "nonstat")
    }
  } else {
    all_pars <- coef(x, complete = TRUE)
    free_pars <- coef(x, complete = FALSE)
    which_free <- which(all_pars %in% free_pars)
    all_pars[which_free] <- free_pars
    n_all_pars <- length(all_pars)
    # If n_all_pars = 2 then model = "gp".
    # If n_all_pars = 3 then model = "pp".
    if (n_all_pars == 2) {
      class(res) <- c("lax", "chandwich", "evd", "pot", "gpd")
    } else {
      class(res) <- c("lax", "chandwich", "evd", "pot", "pp")
    }
  }
  return(res)
}
