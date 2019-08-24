#' Internal lax functions
#'
#' Internal lax functions
#' @details
#' These functions are not intended to be called by the user.
#' @name lax-internal
#' @keywords internal
NULL

#' @keywords internal
#' @rdname lax-internal
adj_object <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # Find all available methods for x
  find_methods_fn <- function(i) as.vector(utils::methods(class = class(x)[i]))
  all_methods <- unlist(sapply(1:length(class(x)), find_methods_fn))
  #
  # Set logLikVec (if a method exists) ----------
  #
  has_logLikVec_method <- paste0("logLikVec.", class(x)) %in% all_methods
  if (!any(has_logLikVec_method)) {
    stop("A logLikVec method must be available for x")
  }
  loglik_fn <- function(pars, fitted_object, ...) {
    return(logLikVec(fitted_object, pars = pars))
  }
  #
  # Set H, but not if use_vcov = FALSE or no vcov method exists ----------
  #
  if (!use_vcov) {
    H <- NULL
  } else {
    # Check whether a vcov method exists for object x
    has_vcov_method <- paste0("vcov.", class(x)) %in% all_methods
    if (any(has_vcov_method)) {
      H <- -solve(vcov(x))
    } else {
      H <- NULL
    }
  }
  #
  # Set mle and nobs ----------
  #
  mle <- coef(x)
  n_obs <- nobs(x)
  #
  # Set V, using meat() or meatCL() from the sandwich package ----------
  #
  if (is.null(cluster)) {
    V <- sandwich::meat(x, fitted_object = x, loglik_fn = loglik_fn,
                        ...) * n_obs
  } else {
    V <- sandwich::meatCL(x, cluster = cluster, fitted_object = x,
                          loglik_fn = loglik_fn, ...) * n_obs
  }
  # We don't pass cluster because it would only be used in the estimation of
  # V: we have already estimated V using sandwich::meat() or sandwich::meatCL()
  res <- chandwich::adjust_loglik(loglik = loglik_fn,
                                  fitted_object = x,
                                  p = length(mle),
                                  par_names = names(mle),
                                  name = paste(class(x), collapse = "_"),
                                  mle = mle, H = H, V = V)
  class(res) <- c("lax", "chandwich")
  return(res)
}


#' @keywords internal
#' @rdname lax-internal
ismev_ppp <- function (a, npy) {
  u <- a[4]
  la <- 1 - exp(-(1 + (a[3] * (u - a[1]))/a[2])^(-1/a[3])/npy)
  sc <- a[2] + a[3] * (u - a[1])
  xi <- a[3]
  c(la, sc, xi)
}

#' @keywords internal
#' @rdname lax-internal
return_level_gev <- function(x, m, level, npy, prof, inc, type) {
  # MLE and symmetric conf% CI for the return level
  rl_sym <- gev_rl_CI(x, m, level, npy, type)
  # Extract SE
  rl_se <- rl_sym["se"]
  # Remove SE
  rl_sym <- rl_sym[c("lower", "mle", "upper")]
  if (!prof) {
    return(list(rl_sym = rl_sym, rl_prof = NULL, rl_se = rl_se))
  }
  temp <- gev_rl_prof(x, m, level, npy, inc, type, rl_sym)
  return(list(rl_sym = rl_sym, rl_prof = temp$rl_prof, rl_se = rl_se,
              max_loglik = logLik(x), crit = temp$crit,
              for_plot = temp$for_plot))
}

#' @keywords internal
#' @rdname lax-internal
gev_rl_prof <- function(x, m, level, npy, inc, type, rl_sym) {
  if (is.null(inc)) {
    inc <- (rl_sym["upper"] - rl_sym["lower"]) / 100
  }
  p <- 1 / (m * npy)
  # Calculates the negated profile loglikelihood of the m-year return level
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
  rl_prof <- c(lower = low_lim, rl_mle, upper = up_lim)
  return(list(rl_prof = rl_prof, crit = conf_line,
              for_plot = cbind(ret_levs = ret_levs, prof_loglik = prof_lik)))
}

#' @keywords internal
#' @rdname lax-internal
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
  res <- c(lower = rl_lower, mle = rl_mle, upper = rl_upper, se = rl_se)
  return(res)
}

# ====================== box_cox_deriv ==========================

#' @keywords internal
#' @rdname lax-internal
box_cox_deriv <- function (x, lambda = 1, lambda_tol = 1e-6,
                           poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda.tol
  #
  # Returns:
  #   A numeric vector.  The derivative with respect to lambda of
  #     (x^lambda - 1) / lambda
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * log(x) - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(log(x) ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}
