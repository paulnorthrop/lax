#' Maximum-likelihood Fitting of the GEV Distribution
#'
#' This is a slightly modified version of \code{\link[ismev]{gev.fit}}.
#' The modification is to add to the returned object the arguments
#' \code{xdat, ydat, mulink, siglink, shlink} and matrices
#' \code{mumat, sigmat, shmat} giving the respective regression design matrices
#' for the location, scale and shape parameters of the model.
#'
#' @examples
#' # We need the evd package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#' if (got_ismev) {
#'   fit1 <- gev.fit(revdbayes::portpirie, show = FALSE)
#'   ls(fit1)
#'   fit2 <- oogev.fit(revdbayes::portpirie, show = FALSE)
#'   ls(fit2)
#' }
#' @export
oogev.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                       mulink = identity, siglink = identity,
                       shlink = identity, muinit = NULL, siginit = NULL,
                       shinit = NULL, show = TRUE, method = "Nelder-Mead",
                       maxit = 10000, ...) {
  z <- list()
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
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
  x <- optim(init, gev.lik, hessian = TRUE, method = method,
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
