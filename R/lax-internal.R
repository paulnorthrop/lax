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
  # Set H, but not if use_cov = FALSE or no vcov method exists ----------
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
