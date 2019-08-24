#' Maximum-likelihood (Re-)Fitting using the ismev package
#'
#' These are a slightly modified versions of the \code{\link[ismev]{gev.fit}},
#' \code{\link[ismev]{gpd.fit}} and \code{\link[ismev]{pp.fit}}
#' functions in the \code{\link[ismev]{ismev}} package.
#' The modification is to add to the returned object regression design matrices
#' for the parameters of the model.  That is,
#' \code{xdat, ydat, mulink, siglink, shlink} and matrices
#' \code{mumat, sigmat, shmat} for the location, scale and shape parameters
#' \code{\link[ismev]{gev.fit}} and \code{\link[ismev]{pp.fit}} and
#' \code{xdat}, \code{ydat, siglink, shlink} and matrices
#' \code{sigmat, shmat} for the scale and shape parameters for
#' \code{\link[ismev]{gpd.fit}}.
#' @inheritParams ismev::gev.fit
#' @inheritParams ismev::gpd.fit
#' @inheritParams ismev::pp.fit
#' @references Heffernan, J. E. and Stephenson, A. G. (2018). ismev: An
#'   Introduction to Statistical Modeling of Extreme Values.
#'   R package version 1.42.
#'   \url{https://CRAN.R-project.org/package=ismev}.
#' @examples
#' # We need the ismev package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#' if (got_ismev) {
#'   library(ismev)
#'   fit1 <- gev.fit(revdbayes::portpirie, show = FALSE)
#'   ls(fit1)
#'   fit2 <- gev_refit(revdbayes::portpirie, show = FALSE)
#'   ls(fit2)
#'
#'   data(rain)
#'   fit1 <- gpd.fit(rain, 10)
#'   ls(fit1)
#'   fit2 <- gpd_refit(rain, 10)
#'   ls(fit2)
#'
#'   fit1 <- pp.fit(rain, 10, show = FALSE)
#'   ls(fit1)
#'   fit2 <- pp_refit(rain, 10, show = FALSE)
#'   ls(fit2)
#' }
#' @name ismev_refits
NULL
## NULL

#' @rdname ismev_refits
#' @export
gev_refit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                       mulink = identity, siglink = identity,
                       shlink = identity, muinit = NULL, siginit = NULL,
                       shinit = NULL, show = TRUE, method = "Nelder-Mead",
                       maxit = 10000, ...) {
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  in2 <- sqrt(6 * stats::var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
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
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  init <- c(muinit, siginit, shinit)
  gev.lik <- function(a) {
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    y <- (xdat - mu)/sc
    y <- 1 + xi * y
    if (any(y <= 0) || any(sc <= 0))
      return(10^6)
    sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
  }
  x <- stats::optim(init, gev.lik, hessian = TRUE, method = method,
                    control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  if (show) {
    if (z$trans)
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv)
      print(z[c(5, 7, 9)])
  }
  z$xdat <- xdat
  z$ydat <- ydat
  z$mumat <- mumat
  z$sigmat <- sigmat
  z$shmat <- shmat
  z$mulink <- mulink
  z$siglink <- siglink
  z$shlink <- shlink
  class(z) <- "gev.fit"
  invisible(z)
}

#' @rdname ismev_refits
#' @export
gpd_refit <- function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL,
                       shl = NULL, siglink = identity, shlink = identity,
                       siginit = NULL, shinit = NULL, show = TRUE,
                       method = "Nelder-Mead", maxit = 10000, ...) {
  z <- list()
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold))
    stop("`threshold' cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1)
    z$trans <- TRUE
  xdatu <- xdat[xdat > u]
  xind <- (1:n)[xdat > u]
  u <- u[xind]
  in2 <- sqrt(6 * stats::var(xdatu))/pi
  in1 <- mean(xdatu, na.rm = TRUE) - 0.57722 * in2
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(siginit))
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
    if (is.null(siginit))
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(shinit))
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
    if (is.null(shinit))
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(siginit, shinit)
  z$model <- list(sigl, shl)
  z$link <- deparse(substitute(c(siglink, shlink)))
  z$threshold <- threshold
  z$nexc <- length(xdatu)
  z$data <- xdatu
  gpd.lik <- function(a) {
    sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0)
      l <- 10^6
    else {
      if (min(y) <= 0)
        l <- 10^6
      else {
        l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
      }
    }
    l
  }
  x <- stats::optim(init, gpd.lik, hessian = TRUE, method = method,
                    control = list(maxit = maxit, ...))
  sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))
  z$conv <- x$convergence
  z$nllh <- x$value
  z$vals <- cbind(sc, xi, u)
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$rate <- length(xdatu)/n
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$n <- n
  z$npy <- npy
  z$xdata <- xdat
  if (show) {
    if (z$trans)
      print(z[c(2, 3)])
    if (length(z[[4]]) == 1)
      print(z[4])
    print(z[c(5, 7)])
    if (!z$conv)
      print(z[c(8, 10, 11, 13)])
  }
  z$xdat <- xdat
  z$ydat <- ydat
  z$sigmat <- sigmat
  z$shmat <- shmat
  z$siglink <- siglink
  z$shlink <- shlink
  class(z) <- "gpd.fit"
  invisible(z)
}

#' @rdname ismev_refits
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
  z$gpd <- apply(z$vals, 1, ismev_ppp, npy)
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
