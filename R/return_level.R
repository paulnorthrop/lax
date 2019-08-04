# Can I just reparameterise a chandwich object?
# ... or use the plot.chandwich emthod

# GEV, GP, PP
# npy
# plot.retlevs()
# lax-internal.R

#' Return levels
#'
#' Add description
#'
#' @examples
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata)
#'   adj_fgev <- alogLik(M1)
#'   summary(adj_fgev)
#' }
#'
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
#'   # An example from the ismev::gev.fit documentation
#'   gev_fit <- gev.fit(revdbayes::portpirie, show = FALSE)
#'   adj_gev_fit <- alogLik(gev_fit)
#'   return_level(adj_gev_fit)
#'   ismev::gev.prof(gev_fit, m = 100, xlow = 4.45, xup = 5.5)
#' }
#' @export
return_level <- function(x, m = 100, level = 0.95, npy = 1, inc = NULL,
                         type = c("vertical", "cholesky", "spectral", "none")) {
  if (!inherits(x, "lax")) {
    stop("use only with \"lax\" objects")
  }
  if (!inherits(x, "stat")) {
    stop("use only with stationary extreme value models")
  }
  type <- match.arg(type)
  if (inherits(x, "gev")) {
    temp <- return_level_gev(x, m, level, npy, inc, type)
  }
  return(temp)
}

return_level_gev <- function(x, m, level, npy, inc, type) {
  # MLE and symmetric conf% CI for the return level
  rl_sym <- gev_rl_CI(x, m, level, npy, type)
  # Set inc (if it hasn't been supplied)
  rl_prof <- gev_rl_prof(x, m, level, npy, inc, type, rl_sym)
  return(list(rl_sym = rl_sym, rl_prof = rl_prof))
}

gev_rl_prof <- function(x, m, level, npy, inc, type, rl_sym) {
  if (is.null(inc)) {
    inc <- (rl_sym$upper - rl_sym$lower) / 100
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
  rl_mle <- rl_sym$mle
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
    opt <- optim(sol, gev_neg_prof_loglik, method = "BFGS", xp = xp)
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
    opt <- optim(sol, gev_neg_prof_loglik, method = "BFGS",xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- xp
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  sol_low <- sol
#  plot(x, v, type = "l", xlab = paste(round(m, 0), "year return level"),
#       ylab = "profile Log-likelihood")
#  abline(h = max_loglik, col = 4)
#  abline(h = conf_line, col = 4)
  return(list(ret_levs = c(rev(x1), x2), prof_lik = c(rev(v1), v2)))
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
  z_val <- qnorm(1 - (1 - level) / 2)
  rl_lower <- rl_mle - z_val * rl_se
  rl_upper <- rl_mle + z_val * rl_se
  list(mle = rl_mle, lower = rl_lower, upper = rl_upper)
}
