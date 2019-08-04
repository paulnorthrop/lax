#' Evaluate loglikelihood contributions from specific observations
#'
#' Generic function for calculating loglikelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#' @export
logLikVec <- function(x, ...) {
  UseMethod("logLikVec")
}

#' Loglikelihood adjustment of model fits
#'
#' This function is generic.  It performs adjustment of the loglikelihood
#' associated with fited model objects, following
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' Certain classes of extreme value model objects are supported automatically.
#' For details see the \code{alogLik} help pages for the packages:
#' \code{\link[lax]{evd}},
#' \code{\link[lax]{evir}},
#' \code{\link[lax]{extRemes}},
#' \code{\link[lax]{fExtremes}},
#' \code{\link[lax]{ismev}},
#' \code{\link[lax]{POT}},
#' \code{\link[lax]{texmex}}.
#' User-supplied objects can also be supported: the requirements for these
#' objects are explained in \strong{Details}.
#'
#' @inherit adj_object params details return references seealso
#' @export
alogLik <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  UseMethod("alogLik")
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
