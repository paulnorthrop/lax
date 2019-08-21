#' Return level inferences for Stationary Extreme Value Models
#'
#' Calculates point estimates and confidence intervals for \code{m}-observation
#' return levels for \strong{stationary} extreme value fitted model objects
#' returned from \code{\link{alogLik}}.  Two types of interval may be returned:
#' (a) intervals based on approximate large-sample normality of the maximum
#' likelihood estimator for return level, which are symmetric about the
#' estimate, and (b) profile likelihood-based intervals based on an (adjusted)
#' loglikelihood.
#'
#' @param x An object inheriting from class \code{"lax"} returned from
#'   \code{\link{alogLik}}.
#' @param m A numeric scalar.  The return period in units of the number
#'   of observations.  See \strong{Details} for information.
#' @param level A numeric scalar in (0, 1).  The (maximum) confidence level
#'   required for the confidence intervals for the \code{m}-year return level.
#' @param npy A numeric scalar.  The
#' @param prof A logical scalar.  Should we calculate intervals based on
#'   profile log-likelihood?
#' @param inc A numeric scalar. Only relevant if \code{prof = TRUE}. The
#'   increment in return level by which we move upwards and downwards from the
#'   MLE for the return level in the search for the lower and upper confidence
#'   limits.  If this is not supplied then \code{inc} is set to one hundredth
#'   of the length of the symmetric confidence interval for return level.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link[chandwich]{adjust_loglik}}, that is, the type of
#'   adjustment made to the independence loglikelihood function in creating
#'   an adjusted loglikelihood function.
#' @details The input value of \code{m} ... GEV, GP, PP
#'
#'   The profile likelihood-based intervals are calculated by
#'   reparameterising in terms of the \code{m}-year return level and estimating
#'   the values at which the (adjusted) profile log-likelihood reaches
#'   the critical value \code{logLik(x) - 0.5 * stats::qchisq(level, 1)}.
#'   This is achieved by calculating the profile loglikelihood for a sequence
#'   of values of this return level as governed by \code{inc}. Once the profile
#'   loglikelhood drops below the critical value the lower and upper limits are
#'   estimated by interpolating linearly between the cases lying either side of
#'   the critical value. The smaller \code{inc} the more accurate (but slower)
#'   the calculation will be.
#' @return A object (a list) of class \code{"retlev", "lax"} with the
#'   components
#'   \item{rl_sym,rl_prof }{Named numeric vectors containing the respective
#'     lower 100\code{level}\% limit, the MLE and the upper
#'     100\code{level}\% limit.  If \code{prof = FALSE} then \code{rl_prof}
#'     will be missing.}
#'   \item{max_loglik,crit,for_plot }{If \code{prof = TRUE} then
#'     these components will be present, containing respectively: the maximised
#'     loglikelihood; the critical value and a matrix with return levels in
#'     the first column (\code{ret_levs}) and the corresponding values of the
#'     (adjusted) profile loglikelihood (\code{prof_loglik}).}
#'   \item{m,level }{The input values of \code{m} and \code{level}.}
#' @examples
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   set.seed(4082019)
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata)
#'   adj_fgev <- alogLik(M1)
#'   rl <- return_level(adj_fgev)
#'   plot(rl)
#' }
#'
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
#'   # An example from the ismev::gev.fit documentation
#'   gev_fit <- gev.fit(revdbayes::portpirie, show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit)
#'   rl <- return_level(adj_gev_fit)
#'   plot(rl)
#'   ismev::gev.prof(gev_fit, m = 100, xlow = 4.45, xup = 5.5)
#' }
#' @export
return_level <- function(x, m = 100, level = 0.95, npy = 1, prof = TRUE,
                         inc = NULL,
                         type = c("vertical", "cholesky", "spectral", "none")) {
  if (!inherits(x, "lax")) {
    stop("use only with \"lax\" objects")
  }
  if (!inherits(x, "stat")) {
    stop("use only with stationary extreme value models")
  }
  type <- match.arg(type)
  if (inherits(x, "gev")) {
    temp <- return_level_gev(x, m, level, npy, prof, inc, type)
  }
  temp$m <- m
  temp$level <- level
  class(temp) <- c("retlev", "lax")
  return(temp)
}

return_level_gev <- function(x, m, level, npy, prof, inc, type) {
  # MLE and symmetric conf% CI for the return level
  rl_sym <- gev_rl_CI(x, m, level, npy, type)
  if (!prof) {
    return(list(rl_sym = rl_sym, rl_prof = NULL))
  }
  temp <- gev_rl_prof(x, m, level, npy, inc, type, rl_sym)
  return(list(rl_sym = rl_sym, rl_prof = temp$rl_prof, max_loglik = logLik(x),
              crit = temp$crit, for_plot = temp$for_plot))
}

gev_rl_prof <- function(x, m, level, npy, inc, type, rl_sym) {
  if (is.null(inc)) {
    inc <- (rl_sym["upper"] - rl_sym["lower"]) / 100
  }
  p <- 1 / (m * npy)
  # Calculates the negated profile log-likelihood of the m-year return level
  gev_neg_prof_loglik <- function(a, xp) {
    if (a[1] <= 0) {
      return(10 ^ 10)
    }
    mu <- xp - revdbayes::qgev(1 - p, loc = 0, scale = a[1], shape = a[2])
    gev_pars <- c(mu, a[1:2])
    return(-x(gev_pars))
  }
  rl_mle <- rl_sym["mle"]
  max_loglik <- attr(x, "max_loglik")
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- rl_mle
  v2[1] <- v1[1] <- max_loglik
  #
  # Starting from the MLE, we search upwards and downwards until we pass the
  # cutoff for the 100level% confidence interval
  #
  ### Upper tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- attr(x, "MLE")[2:3]
  while (my_val > conf_line){
    xp <- xp + inc
    opt <- stats::optim(sol, gev_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- xp
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  sol_up <- sol
  ### Lower tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- attr(x, "MLE")[2:3]
  while (my_val > conf_line){
    xp <- xp - inc
    opt <- stats::optim(sol, gev_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- xp
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  sol_low <- sol
  #
  # Find the limits of the confidence interval
  #
  prof_lik <- c(rev(v1), v2)
  ret_levs <- c(rev(x1), x2)
  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc + 1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc+1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc+1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  rl_prof <- c(lower = low_lim, mle = rl_mle, upper = up_lim)
  names(rl_prof) <- c("lower", "mle", "upper")
  return(list(rl_prof = rl_prof, crit = conf_line,
              for_plot = cbind(ret_levs = ret_levs, prof_loglik = prof_lik)))
}

gev_rl_CI <- function (x, m, level, npy, type){
  mle <- attr(x, "MLE")
  mu <- mle[1]
  sigma <- mle[2]
  xi <- mle[3]
  if (type == "none") {
    mat <- attr(x, "VC")
  } else {
    mat <- attr(x, "adjVC")
  }
  p <- 1 / (m * npy)
  rl_mle <- revdbayes::qgev(p, loc = mu, scale = sigma, shape = xi,
                            lower.tail = FALSE)
  yp <- -log(1 - p)
  delta <- matrix(0, 3, 1)
  delta[1,] <- 1
  delta[2,] <- revdbayes::qgev(p, loc = 0, scale = 1, shape = xi,
                               lower.tail = FALSE)
  delta[3,] <- sigma * box_cox_deriv(yp, lambda = -xi)
  rl_var <- t(delta) %*% mat %*% delta
  rl_se <- sqrt(rl_var)
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  rl_lower <- rl_mle - z_val * rl_se
  rl_upper <- rl_mle + z_val * rl_se
  res <- c(lower = rl_lower, mle = rl_mle, upper = rl_upper)
  return(res)
}

# ------------------------------- plot.retlev ------------------------------- #

#' Plot diagnostics for a ret_lev object
#'
#' \code{plot} method for an objects of class \code{c("ret_lev", "lax")}.
#'
#' @param x an object of class \code{c("ret_lev", "lax")}, a result of
#'   a call to \code{\link{return_level}}.
#' @param y Not used.
#' @param level A numeric scalar in (0, 1).  The confidence level
#'   required for the confidence intervals for the \code{m}-year return level.
#'   This must be no larger than the value stored in \code{x$level}.
#' @param legend A logical scalar.  Should we add a legend (in the top right
#'   of the plot) that gives the approximate values of the MLE and
#'   100\code{level}\% confidence limits?
#' @param digits An integer. Passed to \code{\link[base:Round]{signif}} to
#'   round the values in the legend.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}}.
#' @return Nothing is returned.
#' @seealso \code{\link{return_level}}.
#' @section Examples:
#' See the examples in \code{\link{return_level}}.
#' @export
plot.retlev <- function(x, y = NULL, level = NULL, legend = TRUE, digits = 3,
                        ...) {
  if (!inherits(x, "retlev")) {
    stop("use only with \"retlev\" objects")
  }
  if (!inherits(x, "lax")) {
    stop("use only with \"lax\" objects")
  }
  if (is.null(x$rl_prof)) {
    stop("No prof loglik info: call return_level() using prof = TRUE")
  }
  my_xlab <- paste0(x$m, "-observation return level")
  my_ylab <- "profile loglikelihood"
  my_plot <- function(x, y, ..., xlab = my_xlab, ylab = my_ylab, type = "l") {
    graphics::plot(x, y, ..., xlab = xlab, ylab = ylab, type = type)
  }
  my_plot(x$for_plot[, "ret_levs"], x$for_plot[, "prof_loglik"], ...)
  hline <- function(x, ..., col = "blue", lty = 2) {
    graphics::abline(h = x, ..., col = col, lty = lty)
  }
  hline(x$max_loglik, ...)
  hline(x$crit, ...)
  # Add a legend, if requested
  if (legend) {
    mle_leg <- paste0("     MLE ", signif(x$rl_prof["mle"], digits))
    conf_leg <- paste0(100 * x$level, "% CI (",
                       signif(x$rl_prof["lower"], digits), ",",
                       signif(x$rl_prof["upper"], digits), ")")
    graphics::legend("topright", legend = c(mle_leg, conf_leg))
  }
  return(invisible())
}

