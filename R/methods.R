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
#' Description
#'
#' @inherit adj_object params details return references seealso
#' @export
alogLik <- function(x, ...) {
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
