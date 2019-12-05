# ================================= extRemes ================================ #

#' Loglikelihood adjustment for extRemes fits
#'
#' S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
#' extreme value model objects returned from the function
#' \code{\link[extRemes]{fevd}} in the
#' \code{\link[extRemes:extRemes-package]{extRemes}} package.
#' The model must have been fitted using maximum likelihood estimation.
#'
#' @inherit alogLik params references
#' @details See \code{\link{alogLik}} for details.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   \code{class(x)} is a vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "extRemes")}.
#'   The remaining 2 components depend on the model that was fitted.
#'   The 4th component is: \code{"gev"} if \code{x$type = "GEV"} or
#'   \code{x$type = "Gumbel"}; \code{"gp"} if \code{x$type = "GP"} or
#'   \code{x$type = "Exponential"}; \code{"pp"} if \code{x$type = "PP"}.
#'   The 5th component is
#'   \code{"stat"} if \code{\link[extRemes]{is.fixedfevd} = TRUE} and
#'   \code{"nonstat"} if \code{\link[extRemes]{is.fixedfevd} = FALSE}.
#' @seealso \code{\link{alogLik}}: loglikelihood adjustment for model fits.
#' @examples
#' # We need the extRemes and distillery packages
#' got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
#' got_distillery <- requireNamespace("distillery", quietly = TRUE)
#'
#' if (got_extRemes & got_distillery) {
#'   library(extRemes)
#'   library(distillery)
#'   # Examples from the extRemes::fevd documentation
#'   data(PORTw)
#'
#'   # GEV
#'   fit0 <- fevd(TMX1, PORTw, units = "deg C", use.phi = TRUE)
#'   adj_fit0 <- alogLik(fit0)
#'   summary(adj_fit0)
#'
#'   # GEV regression
#'   fitPORTstdmax <- fevd(TMX1, PORTw, scale.fun = ~STDTMAX, use.phi = TRUE)
#'   adj_fit1 <- alogLik(fitPORTstdmax)
#'   summary(adj_fit1)
#'   fitPORTstdmax2 <- fevd(TMX1, PORTw, location.fun = ~STDTMAX,
#'                          scale.fun = ~STDTMAX, use.phi = TRUE)
#'   adj_fit2 <- alogLik(fitPORTstdmax2)
#'   summary(adj_fit2)
#'   anova(adj_fit0, adj_fit1)
#'   anova(adj_fit1, adj_fit2)
#'   anova(adj_fit0, adj_fit2)
#'   anova(adj_fit0, adj_fit1, adj_fit2)
#'
#'   # Gumbel
#'   fit0 <- fevd(TMX1, PORTw, type = "Gumbel", units = "deg C")
#'   adj_fit0 <- alogLik(fit0)
#'   summary(adj_fit0)
#'
#'   # GP
#'   data(damage)
#'   fit1 <- fevd(Dam, damage, threshold = 6, type = "GP",
#'                time.units = "2.05/year")
#'   adj_fit1 <- alogLik(fit1)
#'   summary(adj_fit1)
#'
#'   # Exponential
#'   fit0 <- fevd(Dam, damage, threshold = 6, type="Exponential",
#'                time.units = "2.05/year")
#'   adj_fit0 <- alogLik(fit0)
#'   summary(adj_fit0)
#'
#'   # GP non-constant threshold
#'   data(Fort)
#'   fit <- fevd(Prec, Fort, threshold = 0.475,
#'               threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
#'               type = "GP")
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#'
#'   # Exponential non-constant threshold
#'   fit <- fevd(Prec, Fort, threshold = 0.475,
#'               threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
#'               type = "Exponential")
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#'
#'   # PP model
#'   fit <- fevd(Prec, Fort, threshold = 0.475, type = "PP", units = "inches")
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#'
#'   # PP non-constant threshold
#'   fit <- fevd(Prec, Fort, threshold = 0.475,
#'               threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
#'               type = "PP")
#'   adj_fit <- alogLik(fit)
#'   summary(adj_fit)
#' }
#' @name extRemes
NULL
## NULL

#' @rdname extRemes
#' @export
alogLik.fevd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  if (x$method != "MLE") {
    stop("Loglikelihood adjustment is only relevant when method = ''MLE''")
  }
  # List of extRemes objects supported
  supported_by_lax <- list(fevd = c("fevd"))
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
  if (x$type == "GEV" || x$type == "Gumbel") {
    class(x) <- "extRemes_gev"
  } else if (x$type == "GP" || x$type == "Exponential") {
    class(x) <- "extRemes_gp"
  } else if (x$type == "PP") {
    class(x) <- "extRemes_pp"
  }
  # Call adj()_object to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (x$type == "GEV" || x$type == "Gumbel") {
    if (extRemes::is.fixedfevd(x)) {
      class(res) <- c("lax", "chandwich", "extRemes", "gev", "stat")
    } else {
      class(res) <- c("lax", "chandwich", "extRemes", "gev", "nonstat")
    }
  } else if (x$type == "GP" || x$type == "Exponential") {
    if (extRemes::is.fixedfevd(x)) {
      class(res) <- c("lax", "chandwich", "extRemes", "gp", "stat")
    } else {
      class(res) <- c("lax", "chandwich", "extRemes", "gp", "nonstat")
    }
  } else if (x$type == "PP") {
    if (extRemes::is.fixedfevd(x)) {
      class(res) <- c("lax", "chandwich", "extRemes", "pp", "stat")
    } else {
      class(res) <- c("lax", "chandwich", "extRemes", "pp", "nonstat")
    }
  }
  return(res)
}
