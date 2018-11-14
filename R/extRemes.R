# ================================= extRemes ================================ #

#' Loglikelihood adjustment of extRemes fits
#'
#' Description
#'
#' @inherit adj_object params details return references seealso
#' @examples
#' # We need the extRemes and distillery packages
#' got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
#' got_distillery <- requireNamespace("distillery", quietly = TRUE)
#'
#' if (got_extRemes & got_distillery) {
#'   library(extRemes)
#'   # Examples from the extRemes::fevd documentation
#'   data(PORTw)
#'
#'   # GEV
#'   fit1 <- fevd(TMX1, PORTw, units = "deg C")
#'   adj_fit1 <- alogLik(fit1)
#'   summary(adj_fit1)
#'
#'   # GEV regression
#'   fitPORTstdmax <- fevd(PORTw$TMX1, PORTw, scale.fun=~STDTMAX, use.phi=TRUE)
#'   adj_fit <- alogLik(fitPORTstdmax)
#'   summary(adj_fit)
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
#'   # GP regression
#'   data(Fort)
#'   fit <- fevd(Prec, Fort, threshold=0.475,
#'               threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
#'               type="GP")
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
    stop("Loglikelihood adjustment is only relevant when method = ''mle''")
  }
  # List of extRemes objects supported
  supported_by_oolax <- list(fevd = c("fevd"))
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
  if (x$type == "GEV" || x$type == "Gumbel") {
    if (extRemes::is.fixedfevd(x)) {
      class(res) <- c("oolax", "chandwich", "extRemes", "gev", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "extRemes", "gev", "nonstat")
    }
  } else if (x$type == "GP" || x$type == "Exponential") {
    if (extRemes::is.fixedfevd(x)) {
      class(res) <- c("oolax", "chandwich", "extRemes", "gp", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "extRemes", "gp", "nonstat")
    }
  }
  return(res)
}
