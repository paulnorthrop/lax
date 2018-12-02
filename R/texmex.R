# ================================== texmex ================================= #

#' Loglikelihood adjustment of texmex fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment of fitted
#' extreme value model objects produced by the
#' \code{\link[texmex:texmex-package]{texmex}} package.
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the texmex package
#' got_texmex <- requireNamespace("texmex", quietly = TRUE)
#' if (got_texmex) {
#'   library(texmex)
#'   # Examples from the texmex::evm documentation
#'
#'   # GEV
#'   mod <- evm(SeaLevel, data = portpirie, family = gev)
#'   adj_mod <- alogLik(mod)
#'   summary(adj_mod)
#'
#'   # GP
#'   mod <- evm(rain, th = 30)
#'   adj_mod <- alogLik(mod)
#'   summary(adj_mod)
#'   mod <- evm(rain, th = 30, cov = "sandwich")
#'   mod$se
#'   vcov(adj_mod)
#'   vcov(mod)
#'
#'   # EGP3
#'   mod <- evm(rain, th = 30, family = egp3)
#'   adj_mod <- alogLik(mod)
#'   summary(adj_mod)
#'
#'   got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'   if (got_ismev) {
#'     # GEV regression
#'     # An example from page 113 of Coles (2001)
#'     data(fremantle)
#'     new_fremantle <- fremantle
#'     # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
#'     new_fremantle[, "Year"] <- new_fremantle[, "Year"] - 1896
#'     evm_fit <- evm(y = SeaLevel, data = new_fremantle, family = gev,
#'                    mu = ~ Year + SOI)
#'     adj_evm_fit <- alogLik(evm_fit)
#'     summary(adj_evm_fit)
#'
#'     # GP regression
#'     # An example from page 119 of Coles (2001)
#'     n_rain <- length(rain)
#'     rain_df <- data.frame(rain = rain, time = 1:n_rain / n_rain)
#'     evm_fit <- evm(y = rain, data = rain_df, family = gpd, th = 30,
#'                    phi = ~ time)
#'     adj_evm_fit <- alogLik(evm_fit)
#'     summary(adj_evm_fit)
#'     evm_fit <- evm(y = rain, data = rain_df, family = gpd, th = 30,
#'                    phi = ~ time, cov = "sandwich")
#'     evm_fit$se
#'     vcov(adj_evm_fit)
#'     vcov(evm_fit)
#'   }
#' }
#' @name texmex
NULL
## NULL

#' @rdname texmex
#' @export
alogLik.evmOpt <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of texmex objects supported
  supported_by_oolax <- list(texmex_evmOpt = c("evmOpt"))
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
  n_pars <- length(coef(x))
  if (x$family$name == "GEV") {
    if (n_pars == 3) {
      class(res) <- c("oolax", "chandwich", "texmex", "gev", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "texmex", "gev", "nonstat")
    }
  } else if (x$family$name == "GPD") {
    if (n_pars == 2) {
      class(res) <- c("oolax", "chandwich", "texmex", "gpd", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "texmex", "gpd", "nonstat")
    }
  } else if (x$family$name == "GPD") {
    if (n_pars == 3) {
      class(res) <- c("oolax", "chandwich", "texmex", "egp3", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "texmex", "egp3", "nonstat")
    }
  }
  return(res)
}
