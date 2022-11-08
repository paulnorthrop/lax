#' Internal lax functions
#'
#' Internal lax functions
#' @details
#' These functions are not intended to be called by the user.
#' @name lax-internal
#' @keywords internal
NULL

# ====================== Log-likelihood adjustment function ================= #

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
  # If cluster was supplied then overwrite the default (1,2, ...) returned by
  # chandwich::adjust_loglik()
  if (!is.null(cluster)) {
    attr(res, "cluster") <- cluster
  }
  # Add the original fitted model object as an attribute
  attr(res, "original_fit") <- x
  class(res) <- c("lax", "chandwich")
  return(res)
}

# ====================== GEV return levels functions ======================== #

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
gev_rl_CI <- function (x, m, level, npy, type) {
  mle <- attr(x, "MLE")
  mu <- mle[1]
  sigma <- mle[2]
  xi <- mle[3]
  if (type == "none") {
    mat <- attr(x, "VC")
  } else {
    mat <- attr(x, "adjVC")
  }
  p <- 1 - (1 - 1 / m) ^ (1 / npy)
  rl_mle <- revdbayes::qgev(p, loc = mu, scale = sigma, shape = xi,
                            lower.tail = FALSE)
  yp <- -log(1 - p)
  delta <- matrix(0, 3, 1)
  delta[1, ] <- 1
  delta[2, ] <- revdbayes::qgev(p, loc = 0, scale = 1, shape = xi,
                               lower.tail = FALSE)
  delta[3, ] <- sigma * box_cox_deriv(yp, lambda = -xi)
  rl_var <- t(delta) %*% mat %*% delta
  rl_se <- sqrt(rl_var)
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  rl_lower <- rl_mle - z_val * rl_se
  rl_upper <- rl_mle + z_val * rl_se
  res <- c(lower = rl_lower, mle = rl_mle, upper = rl_upper, se = rl_se)
  return(res)
}

#' @keywords internal
#' @rdname lax-internal
gev_rl_prof <- function(x, m, level, npy, inc, type, rl_sym) {
  if (is.null(inc)) {
    inc <- (rl_sym["upper"] - rl_sym["lower"]) / 100
  }
  p <- (1 - 1 / m) ^ (1 / npy)
  # Calculates the negated profile loglikelihood of the m-year return level
  gev_neg_prof_loglik <- function(a, xp) {
    # a[1] is sigma, a[2] is xi
    # Check that sigma is positive
    if (a[1] <= 0) {
      return(10 ^ 10)
    }
    mu <- xp - revdbayes::qgev(p, loc = 0, scale = a[1], shape = a[2])
    gev_pars <- c(mu, a[1:2])
    return(-x(gev_pars, type = type))
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

# ==================== Binomial-GP return levels functions ================== #

#' @keywords internal
#' @rdname lax-internal
return_level_bingp <- function(x, m, level, npy, prof, inc, type, npy_given) {
  # Extract the threshold and npy from the original fitted model object
  # The names may vary between packages
  if (inherits(x, "ismev")) {
    u <- attr(x, "original_fit")$threshold
    if (!npy_given) {
      npy <- attr(x, "original_fit")$npy
    }
    if (is.na(npy)) {
      stop("'npy' has not been supplied")
    }
  }
  # MLE and symmetric conf% CI for the return level
  rl_sym <- bingp_rl_CI(x, m, level, npy, type, u)
  # Extract SE
  rl_se <- rl_sym["se"]
  # Remove SE
  rl_sym <- rl_sym[c("lower", "mle", "upper")]
  if (!prof) {
    return(list(rl_sym = rl_sym, rl_prof = NULL, rl_se = rl_se))
  }
  temp <- bingp_rl_prof(x, m, level, npy, inc, type, rl_sym, u)
  return(list(rl_sym = rl_sym, rl_prof = temp$rl_prof, rl_se = rl_se,
              max_loglik = logLik(x), crit = temp$crit,
              for_plot = temp$for_plot))
}

#' @keywords internal
#' @rdname lax-internal
bingp_rl_CI <- function (x, m, level, npy, type, u) {
  # Check that binom = TRUE was used in the call to aloglik()
  bin_object <- attr(x, "pu_aloglik")
  if (is.null(bin_object)) {
    stop("The argument ''binom = TRUE'' must be used when calling alogLik()")
  }
  # Extract information from the binomial inference
  pu <- attr(bin_object, "MLE")
  if (type == "none") {
    bin_mat <- attr(bin_object, "VC")
  } else {
    bin_mat <- attr(bin_object, "adjVC")
  }
  # Extract information from the GP inference
  mle <- attr(x, "MLE")
  sigmau <- mle[1]
  xi <- mle[2]
  if (type == "none") {
    gp_mat <- attr(x, "VC")
  } else {
    gp_mat <- attr(x, "adjVC")
  }
  # If it exists, extract information from the extremal index inference
  theta_exists <- !is.null(attr(x, "theta"))
  if (theta_exists) {
    theta_info <- attr(x, "theta")
    theta <- theta_info$theta
    theta_se <- theta_info$se
    theta_mat <- theta_se ^ 2
  } else {
    theta <- 1
  }
  # Create the covariance matrix for all 3 (or 4) parameters:
  # (pu, sigmau, xi) or (pu, sigmau, xi, theta)
  npars <- ifelse(theta_exists, 4, 3)
  mat <- matrix(0, npars, npars)
  mat[1, 1] <- bin_mat
  mat[2:3, 2:3] <- gp_mat
  # pmnpy is approximately equal to 1 / (m * npy * theta)
  con <- 1 - 1 / m
  nt <- npy * theta
  pmnpy <- 1 - con ^ (1 / nt)
  p <- pmnpy / pu
  rp <- 1 / p
  rl_mle <- revdbayes::qgp(p, loc = u, scale = sigmau, shape = xi,
                           lower.tail = FALSE)
  delta <- matrix(0, npars, 1)
  delta[1, ] <- sigmau * pu ^ (xi - 1) / pmnpy ^ xi
  delta[2, ] <- revdbayes::qgp(p, loc = 0, scale = 1, shape = xi,
                               lower.tail = FALSE)
  delta[3, ] <- sigmau * box_cox_deriv(rp, lambda = xi)
  # Add information about theta, if it is available
  if (theta_exists) {
    delta[4, ] <- -sigmau * rp ^ (xi - 1) * pu * con ^ (1 / nt) * log(con) /
      (npy * theta ^ 2 * pmnpy ^ 2)
    mat[4, 4] <- theta_mat
  }
  rl_var <- t(delta) %*% mat %*% delta
  rl_se <- sqrt(rl_var)
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  rl_lower <- rl_mle - z_val * rl_se
  rl_upper <- rl_mle + z_val * rl_se
  res <- c(lower = rl_lower, mle = rl_mle, upper = rl_upper, se = rl_se)
  return(res)
}

#' @keywords internal
#' @rdname lax-internal
bingp_rl_prof <- function(x, m, level, npy, inc, type, rl_sym, u) {
  if (is.null(inc)) {
    inc <- (rl_sym["upper"] - rl_sym["lower"]) / 100
  }
  # Extract the binomial adjusted loglikelihood object
  bin_object <- attr(x, "pu_aloglik")
  if (is.null(bin_object)) {
    stop("The argument ''binom = TRUE'' must be used when calling alogLik()")
  }
  # If it exists, extract the theta adjusted loglikelihood object
  theta_exists <- !is.null(attr(x, "theta"))
  if (theta_exists) {
    theta_info <- attr(x, "theta")
    theta_mle <- theta_info$theta
    theta_loglik <- function(tval) {
      theta_list <- c(list(theta = tval), theta_info$ss)
      return(do.call(kgaps_loglik, theta_list))
    }
  }
  # Extract the MLEs of (pu, sigmau, xi)
  bin_mle <- attr(bin_object, "MLE")
  gp_mle <- attr(x, "MLE")
  # Function to calculate the negated adjusted binomial-GP loglikelihood
  if (theta_exists) {
    bingp_mle <- c(bin_mle, gp_mle, theta = theta_mle)
    bingp_negloglik <- function(pars, type) {
      gp_negloglik <- -x(pars[2:3], type = type)
      bin_negloglik <- -bin_object(pars[1], type = type)
      theta_negloglik <- -theta_loglik(pars[4])
      return(gp_negloglik + bin_negloglik + theta_negloglik)
    }
    # Calculates the negated profile loglikelihood of the m-year return level
    bingp_neg_prof_loglik <- function(a, xp) {
      # a[1] is pu, a[2] is xi, a[3] is theta
      # Check that pu is in (0, 1) and theta is in (0, 1]
      if (a[1] <= 0 || a[1] >= 1 || a[3] <= 0 || a[3] > 1) {
        return(10 ^ 10)
      }
      pmnpy <- 1 - (1 - 1 / m) ^ (1 / (npy * a[3]))
      p <- pmnpy / a[1]
      sigmau <- (xp - u) / revdbayes::qgp(p, loc = 0, scale = 1, shape = a[2],
                                          lower.tail = FALSE)
      # Check that sigmau is positive
      if (sigmau <= 0) {
        return(10 ^ 10)
      }
      bingp_pars <- c(a[1], sigmau, a[2], a[3])
      return(bingp_negloglik(bingp_pars, type = type))
    }
    max_loglik <- attr(x, "max_loglik") + attr(bin_object, "max_loglik") +
       theta_loglik(theta_mle)
#      This could be used once exdex 1.0.2 is on CRAN
#      theta_info$max_loglik
  } else {
    bingp_mle <- c(bin_mle, gp_mle)
    bingp_negloglik <- function(pars, type) {
      gp_negloglik <- -x(pars[2:3], type = type)
      bin_negloglik <- -bin_object(pars[1], type = type)
      return(gp_negloglik + bin_negloglik)
    }
    pmnpy <- 1 - (1 - 1 / m) ^ (1 / npy)
    # Calculates the negated profile loglikelihood of the m-year return level
    bingp_neg_prof_loglik <- function(a, xp) {
      # a[1] is pu, a[2] is xi
      # Check that pu is in (0, 1)
      if (a[1] <= 0 || a[1] >= 1) {
        return(10 ^ 10)
      }
      p <- pmnpy / a[1]
      sigmau <- (xp - u) / revdbayes::qgp(p, loc = 0, scale = 1, shape = a[2],
                                          lower.tail = FALSE)
      # Check that sigmau is positive
      if (sigmau <= 0) {
        return(10 ^ 10)
      }
      bingp_pars <- c(a[1], sigmau, a[2])
      return(bingp_negloglik(bingp_pars, type = type))
    }
    max_loglik <- attr(x, "max_loglik") + attr(bin_object, "max_loglik")
  }
  rl_mle <- rl_sym["mle"]
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
  sol <- bingp_mle[-2]
  while (my_val > conf_line){
    xp <- xp + inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
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
  sol <- bingp_mle[-2]
  while (my_val > conf_line){
    xp <- xp - inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
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

# ============================== box_cox_deriv ============================== #

#' @keywords internal
#' @rdname lax-internal
box_cox_deriv <- function(x, lambda = 1, lambda_tol = 1 / 50,
                          poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda_tol
  #
  # Returns:
  #   A numeric vector.  The derivative with respect to lambda of
  #     (x^lambda - 1) / lambda
  #
  lnx <- log(x)
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * lnx - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(lnx ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}

# ============== Used in (ismev) pp_refit() for GPD parameters ============== #

#' @keywords internal
#' @rdname lax-internal
ismev_ppp <- function (a, npy) {
  u <- a[4]
  la <- 1 - exp(-(1 + (a[3] * (u - a[1]))/a[2])^(-1/a[3])/npy)
  sc <- a[2] + a[3] * (u - a[1])
  xi <- a[3]
  c(la, sc, xi)
}

# =============================== kgaps_loglik ============================== #
# Included because this is not exported from exdex
# The argument n_kgaps is not used here but it is included because it is
# included in the list returned by exdex::kgaps_stat()

#' @keywords internal
#' @rdname lax-internal
kgaps_loglik <- function(theta, N0, N1, sum_qs, n_kgaps){
  if (theta < 0 || theta > 1) {
    return(-Inf)
  }
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}
