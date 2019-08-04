#' Maximum-likelihood Fitting for the Point Process Model
#'
#' This is a slightly modified versions of the \code{\link[ismev]{gev.fit}}
#' function in the \code{\link[ismev]{ismev}} package.
#' The main modification is to add to the returned object the arguments
#' \code{xdat, ydat, mulink, siglink, shlink} and matrices
#' \code{mumat, sigmat, shmat} giving the respective regression design matrices
#' for the location, scale and shape parameters of the model.  In addition, a
#' bug in the code that sets initial estimates has been corrected: the bug
#' meant that if \code{threshold} is a vector then the optimization hangs.
#' @inheritParams ismev::pp.fit
#' @references Heffernan, J. E. and Stephenson, A. G. (2018). ismev: An
#'   Introduction to Statistical Modeling of Extreme Values.
#'   R package version 1.42.
#'   \url{https://CRAN.R-project.org/package=ismev}.
#' @examples
#' # We need the evd package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#' if (got_ismev) {
#'   data(rain)
#'   fit1 <- pp.fit(rain, 100, show = FALSE)
#'   ls(fit1)
#'   fit2 <- pp_refit(rain, 10, show = FALSE)
#'   ls(fit2)
#' }
#' @export
pp_refit <- function (xdat, threshold, npy = 365, ydat = NULL, mul = NULL,
                      sigl = NULL, shl = NULL, mulink = identity,
                      siglink = identity, shlink = identity, muinit = NULL,
                      siginit = NULL, shinit = NULL, show = TRUE,
                      method = "Nelder-Mead", maxit = 10000, ...) {
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold))
    stop("`threshold' cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1)
    z$trans <- TRUE
  uInd <- xdat > u
  lrate <- sum(uInd)/n
  xdatu <- xdat[uInd]
  in2 <- sqrt(6 * stats::var(xdatu))/pi
  in1 <- mean(xdatu) - 0.57722 * in2
  if (is.null(shinit))
    in3 <- 1e-08
  else in3 <- shinit
  # Change from ismev::pp.fit(). These initial estimates are only sensible if
  # threshold is a scalar.  Using them means that the wooster.temps demo hangs
  if (length(threshold) == 1) {
    in2 <- exp(log(in2) + in3 * log(lrate))
    in1 <- threshold - (in2/in3) * (lrate^(-in3) - 1)
  }
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit))
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit))
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(siginit))
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
    if (is.null(siginit))
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit))
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit))
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(muinit, siginit, shinit)
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  z$threshold <- threshold
  z$npy <- npy
  z$nexc <- length(xdatu)
  z$data <- xdatu
  pp.lik <- function(a) {
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    if (any(sc^uInd <= 0))
      return(10^6)
    if (min((1 + ((xi * (u - mu))/sc))^uInd) < 0) {
      l <- 10^6
    }
    else {
      y <- (xdat - mu)/sc
      y <- 1 + xi * y
      if (min(y) <= 0)
        l <- 10^6
      else l <- sum(uInd * log(sc)) +
        sum(uInd * log(y) * (1/xi + 1)) +
        n/npy * mean((1 + (xi * (u - mu))/sc)^(-1/xi))
    }
    l
  }
  x <- stats::optim(init, pp.lik, hessian = TRUE, method = method,
                    control = list(maxit = maxit, ...))
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$conv <- x$convergence
  z$nllh <- x$value
  z$vals <- cbind(mu, sc, xi, u)
  z$gpd <- apply(z$vals, 1, ppp, npy)
  if (z$trans) {
    z$data <- as.vector((1 + (xi[uInd] * (xdatu - u[uInd])) /
                           z$gpd[2,uInd])^(-1/xi[uInd]))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  if (show) {
    if (z$trans)
      print(z[c(2, 3)])
    if (length(z[[4]]) == 1)
      print(z[4])
    print(z[c(5, 6, 8)])
    if (!z$conv)
      print(z[c(9, 12, 14)])
  }
  z$xdat <- xdat
  z$ydat <- ydat
  z$mumat <- mumat
  z$sigmat <- sigmat
  z$shmat <- shmat
  z$mulink <- mulink
  z$siglink <- siglink
  z$shlink <- shlink
  class(z) <- "pp.fit"
  invisible(z)
}
