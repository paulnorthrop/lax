logLikVec <- function(x, ...) {
  UseMethod("logLikFn")
}

alogLik <- function(x, ...) {
  UseMethod("alogLik")
}

estfun <- function(x, ...) {
  UseMethod("estfun")
}
