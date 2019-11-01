#' Evaluate loglikelihood contributions from specific observations
#'
#' Generic function for calculating loglikelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#' @export
logLikVec <- function(object, ...) {
  UseMethod("logLikVec")
}

#' Loglikelihood adjustment for model fits
#'
#' This function is generic.  It performs adjustment of the loglikelihood
#' associated with fitted model objects, following
#' \href{http://doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' Certain classes of extreme value model objects are supported automatically.
#' For details see the \code{alogLik} help pages for the packages:
#' \code{\link[lax]{evd}},
#' \code{\link[lax]{evir}},
#' \code{\link[lax]{extRemes}},
#' \code{\link[lax]{fExtremes}},
#' \code{\link[lax]{ismev}},
#' \code{\link[lax]{mev}},
#' \code{\link[lax]{POT}},
#' \code{\link[lax]{texmex}}.
#' User-supplied objects can also be supported: the requirements for these
#' objects are explained in \strong{Details}.
#'
#' @param x A fitted model object with certain associated S3 methods.
#'   See \strong{Details}.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   This must have the same length as the vector returned by the
#'   \code{logLikVec} method for an object like \code{x}.
#'   If \code{cluster} is not supplied (i.e. is \code{NULL}) then it is
#'   assumed that each observation forms its own cluster.
#'   See \strong{Details}.
#' @param use_vcov A logical scalar.  Should we use the \code{vcov} S3 method
#'   for \code{x} (if this exists) to estimate the Hessian of the independence
#'   loglikelihood to be passed as the argument \code{H} to
#'   \code{\link[chandwich]{adjust_loglik}}?
#'   Otherwise, \code{H} is estimated inside
#'   \code{\link[chandwich]{adjust_loglik}} using
#'   \code{\link[stats:optim]{optimHess}}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}} (if \code{cluster = NULL}),
#'   or \code{\link[sandwich:vcovCL]{meatCL}} (if \code{cluster} is not
#'   \code{NULL}).
#' @details Object \code{x} \emph{must} have the following S3
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
#'   Loglikelihood adjustment is performed using the
#'   \code{\link[chandwich]{adjust_loglik}} function in the
#'   \code{\link[chandwich]{chandwich}} package.
#'   The relevant arguments to \code{\link[chandwich]{adjust_loglik}}, namely
#'   \code{loglik, mle, H} and \code{V}, are created based on the class of
#'   the object \code{x}.
#'
#'   If a \code{vcov} method is not available, or if \code{use_vcov = FALSE},
#'   then the variance-covariance matrix of the MLE (from which \code{H} is
#'   calculated) is estimated inside \code{\link[chandwich]{adjust_loglik}}
#'   using \code{\link[stats:optim]{optimHess}}.
#'
#'   The \code{sandwich} package is used to estimate the variance matrix
#'   \code{V} of the score vector: \code{\link[sandwich]{meat}} is used if
#'   \code{cluster = NULL}; \code{\link[sandwich:vcovCL]{meatCL}} is used if
#'   \code{cluster} is not \code{NULL}.
#'   If \code{cluster} is \code{NULL} then any arguments of
#'   \code{\link[sandwich:vcovCL]{meatCL}} present in \dots will be ignored.
#'   Similarly, if \code{cluster} is not \code{NULL} then any arguments of
#'   \code{\link[sandwich]{meat}} present in \dots will be ignored.
#'   \code{\link[sandwich]{meat}} and \code{\link[sandwich:vcovCL]{meatCL}}
#'   require an \code{\link[sandwich]{estfun}} method to be available, which,
#'   in the current context, provides matrix of score contributions.
#'   If a bespoke \code{estfun} method is not provided then this is constructed
#'   by estimating the score contributions using \code{\link[numDeriv]{jacobian}}.
#' @return An object inheriting from class \code{"chandwich"}.  See
#'   \code{\link[chandwich]{adjust_loglik}}.
#'
#'   If \code{x} is one of the supported models then \code{class(x)} is a
#'   vector of length 5. The first 3 components are
#'   \code{c("lax", "chandwich", "name_of_package")}, where
#'   \code{"name_of_package"} is the name of the package from which the input
#'   object \code{x} originated.  The remaining 2 components depend on the
#'   model that was fitted.  See the documentation of the relevant package
#'   for details:
#'   \code{\link[lax]{evd}},
#'   \code{\link[lax]{evir}},
#'   \code{\link[lax]{extRemes}},
#'   \code{\link[lax]{fExtremes}},
#'   \code{\link[lax]{ismev}},
#'   \code{\link[lax]{mev}},
#'   \code{\link[lax]{POT}},
#'   \code{\link[lax]{texmex}}.
#'
#'   Otherwise, \code{class(x)} is \code{c("lax", "chandwich", class(x))}.
#'
#'   Objects returned from `aloglik` have `anova`, `coef`, `confint`, `logLik`,
#'   `nobs`, `plot`, `print`, `summary` and `vcov` methods.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://doi.org/10.1093/biomet/asm015}
#' @references Zeleis (2006) Object-Oriented Computation and Sandwich
#'   Estimators.  \emph{Journal of Statistical Software}, \strong{16}, 1-16.
#'   \url{http://doi.org/10.18637/jss.v016.i09}
#' @seealso \code{\link[chandwich]{summary.chandwich}},
#'   \code{\link[chandwich]{plot.chandwich}},
#'   \code{\link[chandwich]{confint.chandwich}},
#'   \code{\link[chandwich]{anova.chandwich}},
#'   \code{\link[chandwich]{coef.chandwich}},
#'   \code{\link[chandwich]{vcov.chandwich}}
#'   and \code{\link[chandwich]{logLik.chandwich}}
#'   for S3 methods for objects of class \code{"chandwich"}.
#'
#'   \code{\link[chandwich]{conf_region}} for confidence regions for
#'   pairs of parameters.
#' @seealso \code{\link[chandwich]{adjust_loglik}} in the
#'   \code{\link[chandwich]{chandwich}} package to adjust a user-supplied
#'   loglikelihood.
#' @seealso \code{\link[sandwich]{meat}} and
#'   \code{\link[sandwich:vcovCL]{meatCL}} in the sandwich package.
#' @section Examples:
#' See the (package-specific) examples in \code{\link{evd}},
#'   \code{\link{evir}}, \code{\link{extRemes}},\code{\link{fExtremes}},
#'   \code{\link{ismev}}, \code{\link{mev}}, \code{\link{POT}} and
#'   \code{\link{texmex}}.
#' @export
alogLik <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  UseMethod("alogLik")
}

#' @export
alogLik.default <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # Call adj_object() to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  class(res) <- c("lax", "chandwich", class(x))
  return(res)
}

#' Sum loglikelihood contributions from individual observations
#'
#' S3 logLik method for logLikVec objects
#'
#' @param object An object of class \code{"logLikVec"} return from a
#'   \code{logLikVec} method.
#' @param ... Further arguments.
#' @export
logLik.logLikVec <- function(object, ...) {
  save_attributes <- attributes(object)
  object <- sum(object)
  attributes(object) <- save_attributes
  class(object) <- "logLik"
  return(object)
}
