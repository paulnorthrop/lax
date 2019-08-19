#' Fits a Poisson point process to the data, an approach sometimes known as
#' peaks over thresholds (POT), and returns an object of class "potd".
#'
#' This is a slightly modified versions of the \code{\link[evir]{pot}}
#' function in the \code{evir} package.
#' The main modification is to add to the returned object the argument
#' \code{data} supplied by the user.  This is added to the returned
#' (list) object with the name \code{input_data}.
#' @param data numeric vector of data, which may have a times attribute
#'   containing (in an object of class \code{"POSIXct"}, or an object that can
#'   be converted to that class; see \code{\link{as.POSIXct}}) the times/dates
#'   of each observation. If no times attribute exists, the data are assumed to
#'   be equally spaced.
#' @param threshold a threshold value (either this or \code{nextremes} must be
#'   given but not both).
#' @param nextremes the number of upper extremes to be used (either this or
#'   \code{threshold} must be given but not both).
#' @param run if the data are to be declustered the run length parameter for
#'   the runs method (see \code{\link[evir]{decluster}}) should be entered
#'   here.
#' @param picture whether or not a picture should be drawn if declustering is
#'   performed.
#' @param ... arguments passed to \code{\link[stats]{optim}}.
#' @references Bernhard Pfaff and Alexander McNeil (2018). evir: Extreme
#'   Values in R. R package version 1.7-4.
#'   \url{https://CRAN.R-project.org/package=evir}.
#' @examples
#' # We need the evir package
#' got_evir <- requireNamespace("evir", quietly = TRUE)
#' if (got_evir) {
#'   library(evir)
#'   data(danish)
#'   out <- pot(danish, 10)
#'   ls(out)
#'   out <- pot_refit(danish, 10)
#'   ls(out)
#' }
#' @export
pot_refit <- function(data, threshold = NA, nextremes = NA, run = NA,
                      picture = TRUE, ...) {
  # Save the input data so that we can return them later
  input_data <- data
  n <- length(as.numeric(data))
  times <- attributes(data)$times
  if (is.null(times)) {
    times <- 1:n
    attributes(data)$times <- times
    start <- 1
    end <- n
    span <- end - start
  }
  else {
    start <- times[1]
    end <- times[n]
    span <- as.numeric(difftime(as.POSIXlt(times)[n], as.POSIXlt(times)[1],
                                units = "days"))
  }
  if (is.na(nextremes) && is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if (!is.na(nextremes) && !is.na(threshold))
    stop("Enter EITHER a threshold or the number of upper extremes")
  if (!is.na(nextremes))
    threshold <- evir::findthresh(as.numeric(data), nextremes)
  if (threshold > 10) {
    factor <- 10^(floor(log10(threshold)))
    cat(paste("If singularity problems occur divide data",
              "by a factor, perhaps", factor, "\n"))
  }
  exceedances.its <- structure(data[data > threshold], times = times[data >
                                                                       threshold])
  n.exceed <- length(as.numeric(exceedances.its))
  p.less.thresh <- 1 - n.exceed/n
  if (!is.na(run)) {
    exceedances.its <- evir::decluster(exceedances.its, run, picture)
    n.exceed <- length(exceedances.its)
  }
  intensity <- n.exceed/span
  exceedances <- as.numeric(exceedances.its)
  xbar <- mean(exceedances) - threshold
  s2 <- stats::var(exceedances)
  shape0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  extra <- ((length(exceedances)/span)^(-shape0) - 1)/shape0
  betahat <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  scale0 <- betahat/(1 + shape0 * extra)
  loc0 <- 0
  theta <- c(shape0, scale0, loc0)
  negloglik <- function(theta, exceedances, threshold, span) {
    if ((theta[2] <= 0) || (min(1 + (theta[1] * (exceedances -
                                                 theta[3]))/theta[2]) <= 0))
      f <- 1e+06
    else {
      y <- logb(1 + (theta[1] * (exceedances - theta[3]))/theta[2])
      term3 <- (1/theta[1] + 1) * sum(y)
      term1 <- span * (1 + (theta[1] * (threshold - theta[3]))/theta[2])^(-1/theta[1])
      term2 <- length(y) * logb(theta[2])
      f <- term1 + term2 + term3
    }
    f
  }
  fit <- stats::optim(theta, negloglik, hessian = TRUE, ...,
                      exceedances = exceedances, threshold = threshold,
                      span = span)
  if (fit$convergence)
    warning("optimization may not have succeeded")
  par.ests <- fit$par
  varcov <- solve(fit$hessian)
  par.ses <- sqrt(diag(varcov))
  beta <- par.ests[2] + par.ests[1] * (threshold - par.ests[3])
  par.ests <- c(par.ests, beta)
  out <- list(n = length(data), period = c(start, end), data = exceedances.its,
              span = span, threshold = threshold, p.less.thresh = p.less.thresh,
              n.exceed = n.exceed, run = run, par.ests = par.ests,
              par.ses = par.ses, varcov = varcov, intensity = intensity,
              nllh.final = fit$value, converged = fit$convergence,
              input_data = input_data)
  names(out$par.ests) <- c("xi", "sigma", "mu", "beta")
  names(out$par.ses) <- c("xi", "sigma", "mu")
  class(out) <- "potd"
  out
}
