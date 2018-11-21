# =================================== evd =================================== #

#' Loglikelihood adjustment of evd fits
#'
#' Description
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the evd package
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
#'   adj_fgev <- alogLik(M1)
#'   summary(adj_fgev)
#'
#'   # An example from Chandler and Bate (2007)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- rep(c(-1, 1), each = length(y) / 2)
#'   owfit <- evd::fgev(y, nsloc = x)
#'   year <- rep(rownames(chandwich::owtemps), 2)
#'   adj_owfit <- alogLik(owfit, cluster = year)
#'   summary(adj_owfit)
#'
#'   # An example from the evd::fpot documentation
#'   uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fpot(uvdata, 1)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#'   # Fit using the pp model, rather than the gpd
#'   M1 <- evd::fpot(uvdata, 1, model = "pp", npp = 365)
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
  supported_by_oolax <- list(evd_fgev = c("gev", "uvevd", "evd"),
                             evd_fpot = c("pot", "uvevd", "evd"))
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
  if (name_of_class == "evd_fgev") {
    if (is.null(x$nsloc)) {
      class(res) <- c("oolax", "chandwich", "evd", "gev", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "evd", "gev", "nonstat")
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
      class(res) <- c("oolax", "chandwich", "evd", "pot", "gpd")
    } else {
      class(res) <- c("oolax", "chandwich", "evd", "pot", "pp")
    }
  }
  return(res)
}
