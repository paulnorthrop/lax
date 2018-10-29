# ================================== adj_object ============================= #

#' Loglikelihood adjustment of fitted model objects
#'
#' Performs adjustments, for model misspecification or the the presence of
#' cluster dependence, of the loglikelihood associated with fitted
#' model objects, following
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' Certain classes of model objects are supported automatically.
#' User-supplied objects can also be supported: the requirements for these
#' objects are explained in \strong{Details}.
#'
#' @param x A fitted model object with certain associated S3 methods.
#'   See \strong{Details}.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   This must have the same length as the vector returned by the
#'   \code{logLikVec} method for object like \code{x}.
#'   If \code{cluster} is not supplied then it is assumed that
#'   each observation forms its own cluster.
#'
#'   If the sandwich package
#'   \href{http://dx.doi.org/10.18637/jss.v016.i09}{(Zeleis, 2006)}
#'   is used to estimate the quantities required to adjust the loglikelihood
#'   (i.e. \code{use_sanswich = TRUE}) then \code{cluster} determines whether
#'   the variance matrix \code{V} of the score vector is estimated using
#'   \code{\link[sandwich]{meat}} (\code{cluster} is \code{NULL}) or
#'   \code{\link[sandwich:vcovCL]{meatCL}} (\code{cluster} is not \code{NULL}).
#'   See \code{use_sandwich} and \strong{Details}.
#' @param use_sandwich A logical scalar.  Should we use the \code{sandwich}
#'   package \href{http://dx.doi.org/10.18637/jss.v016.i09}{(Zeleis, 2006)}
#'   to estimate the variance \eqn{V} of the score function to be passed as
#'   the argument \code{V} to \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \eqn{V} is based on numerical derivatives, estimated using
#'   the \code{\link[numDeriv:numDeriv-package]{numDeriv}} package.
#'   See \strong{Details} for more information.
#'
#'   The main purpose of \code{use_sandwich} is to enable a check that
#'   equivalent results are obtained using the sandwich package
#'   and the numerical derivatives.
#' @param use_vcov A logical scalar.  Should we use the \code{vcov} S3 method
#'   for \code{x} (if this exists) to estimate the Hessian of the independence
#'   loglikelihood to be passed as the argument \code{H} to
#'   \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \code{H} is estimated inside
#'   \code{\link[chandwich]{adjust_loglik}} using
#'   \code{\link[stats:optim]{optimHess}}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}}, if \code{cluster = NULL},
#'   or \code{\link[sandwich:vcovCL]{meatCL}}, otherwise.
#' @details
#'   Supported classes, and the functions from which they are returned are:
#'   \itemize{
#'     \item{"lm": }{\code{\link[stats]{lm}} in the
#'       \code{\link[stats:stats-package]{stats}} package,}
#'     \item{"glm": }{\code{\link[stats]{glm}} in the
#'       \code{\link[stats:stats-package]{stats}} package}
#'   }
#'   Objects of other classes are supported provided that they have
#'   certain S3 methods.
#'   The class of the object \code{x} \emph{must} have the following S3
#'   methods:
#'   \itemize{
#'     \item{\code{logLikVec: }}{returns a vector of the contributions to the
#'       independence loglikelihood from individual observations;}
#'     \item{\code{coef: }}{returns a vector of model coefficients, see
#'       \code{\link[stats]{coef}};}
#'     \item{\code{nobs: }}{returns the number of (non-missing) observations
#'       used in a model fit, see \code{\link[stats]{nobs}}};
#'   }
#'   and \emph{may} have the following S3 methods
#'   \itemize{
#'     \item{\code{vcov: }}{returns the estimated variance-covariance matrix of
#'       the (main) parameters of a fitted model, see
#'       \code{\link[stats]{vcov}};}
#'     \item{\code{estfun: }}{returns an \eqn{n x k} matrix, in which each
#'       column gives the derivative of the loglikelihood at each of \eqn{n}
#'       observation with respect to the \eqn{k} parameters of the model, see
#'       \code{\link[sandwich]{estfun}}.}
#'   }
#'   If a \code{vcov} method is not available then the variance-covariance
#'   matrix is estimated inside \code{\link[chandwich]{adjust_loglik}}.
#'   If an \code{estfun} method is not available then the matrix is estimated
#'   using \code{\link[numDeriv]{jacobian}}.
#'
#'   Loglikelihood adjustment is performed using the
#'   \code{\link[chandwich]{adjust_loglik}} function in the
#'   \code{\link[chandwich]{chandwich}} package.
#'   The relevant arguments to \code{\link[chandwich]{adjust_loglik}}, namely
#'   \code{loglik, mle, H} and \code{V}, are created based on the class of
#'   the object \code{x}. If \code{use_sandwich = TRUE} then
#'   \code{H} is inferred using \code{\link[sandwich]{bread}} in the
#'   sandwich package and, similarly, \code{V} is inferred using
#'   either \code{\link[sandwich]{meat}}, if \code{cluster = NULL}
#'   and \code{\link[sandwich:vcovCL]{meatCL}}, otherwise.
#'
#'   If \code{cluster} is \code{NULL} then arguments of
#'   \code{\link[sandwich:vcovCL]{meatCL}} present in \dots will be ignored.
#'   Similarly, if \code{cluster} is not \code{NULL} then arguments of
#'   \code{\link[sandwich]{meat}} present in \dots will be ignored.
#' @return An object of class inheriting from \code{"oola"}, which inherits
#'   from the class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'   The attribute \code{"name"} of the returned object is the elements of
#'   \code{class(x)} concatenated into a character scalar and separated
#'   by \code{_}.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @references Zeleis (2006) Object-Oriented Computation and Sandwich
#'   Estimators.  \emph{Journal of Statistical Software}, \strong{16}, 1-16.
#'   \url{http://dx.doi.org/10.18637/jss.v016.i09}
#' @seealso \code{\link[chandwich]{summary.chandwich}},
#'   \code{\link[chandwich]{plot.chandwich}},
#'   \code{\link[chandwich]{confint.chandwich}},
#'   \code{\link[chandwich]{anova.chandwich}},
#'   \code{\link[chandwich]{coef.chandwich}},
#'   \code{\link[chandwich]{vcov.chandwich}}
#'   and \code{\link[chandwich]{logLik.chandwich}}
#'   for S3 methods for objects of class \code{"chandwich"}.
#' @seealso \code{\link[chandwich]{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood.
#' @seealso \code{\link[sandwich]{bread}}, \code{\link[sandwich]{meat}},
#'   \code{\link[sandwich:vcovCL]{meatCL}} and
#'   \code{\link[sandwich]{sandwich}} in the sandwich package.
#' @export
adj_object <- function(x, cluster = NULL, use_sandwich = TRUE,
                       use_vcov = TRUE, ...) {
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
  loglik_fn <- function(pars, fitted_object, contrib = TRUE, ...) {
    return(logLikVec(fitted_object, pars = pars, contrib = contrib))
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
  if (use_sandwich) {
    if (is.null(cluster)) {
      V <- sandwich::meat(x, fitted_object = x, loglik_fn = loglik_fn,
                          ...) * n_obs
    } else {
      V <- sandwich::meatCL(x, cluster = cluster, fitted_object = x,
                            loglik_fn = loglik_fn, ...) * n_obs
    }
  } else {
    V <- NULL
  }
  # We don't pass cluster because it would only be used in the estimation of
  # V: we have already estimated V using sandwich::meat() or sandwich::meatCL()
  res <- chandwich::adjust_loglik(loglik = loglik_fn,
                                  fitted_object = x,
                                  p = length(mle),
                                  par_names = names(mle),
                                  name = paste(class(x), collapse = "_"),
                                  mle = mle, H = H, V = V)
  class(res) <- c("oola", "chandwich")
  return(res)
}
