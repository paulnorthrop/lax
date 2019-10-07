# =============================== mev::fit.gpd ============================== #

# Methods for class mev_gpd

#' @export
logLikVec.mev_gpd <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  #
  # Threshold exceedances (values that lie above the threshold)
  response_data <- object$threshold + object$exceedances
  sigma <- pars[1]
  xi <- pars[2]
  # Calculate the loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = object$threshold, scale = sigma,
                          shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}
