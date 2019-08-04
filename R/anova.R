# =============================== anova.lax =============================== #

#' Comparison of nested models
#'
#' \code{anova} method for objects of class \code{"lax"}.
#' Compares two or more nested models using the adjusted likelihood ratio
#' test statistic (ALRTS) described in Section 3.5 of
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' The nesting must result from the simple constraint that a subset of the
#' parameters of the larger model is held fixed.
#'
#' @param object An object of class \code{"lax"}, returned by
#'   \code{\link{alogLik}}.
#' @param object2 An object of class \code{"chandwich"}, returned by
#'   \code{\link{alogLik}}.
#' @param ... Further objects of class \code{"lax"} and/or arguments
#'   to be passed to \code{\link[chandwich]{anova.chandwich}}.
#'
#' @details The objects of class \code{"lax"} need not be provided in nested
#'   order: they will be ordered inside \code{anova.lax} based on the
#'   values of \code{attr(., "p_current")}.
#' @return An object of class \code{"anova"} inheriting from class
#'  \code{"data.frame"}, with four columns:
#'     \item{Model.Df}{The number of parameters in the model}
#'     \item{Df}{The decrease in the number of parameter compared the model
#'       in the previous row}
#'     \item{ALRTS}{The adjusted likelihood ratio test statistic}
#'     \item{Pr(>ALRTS)}{The p-value associated with the test that the
#'       model is a valid simplification of the model in the previous row.}
#'  The row names are the names of the model objects.
#' @seealso \code{\link[chandwich]{anova.chandwich}}.
#' @seealso \code{\link{alogLik}}.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @examples
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- rep(c(1, -1), each = length(y) / 2)
#'   owfit <- evd::fgev(y, nsloc = x)
#'   year <- rep(1:(length(y) / 2), 2)
#'   small <- alogLik(owfit, cluster = year)
#'
#'   owfit <- evd::fgev(y)
#'   year <- rep(1:(length(y) / 2), 2)
#'   tiny <- alogLik(owfit, cluster = year)
#'
#'   anova(small, tiny)
#' }
#' @export
anova.lax <- function (object, object2, ...) {
  if (missing(object)) {
    stop("model one must be supplied, using object")
  }
  if (missing(object2)) {
    stop("model two must be supplied, using object2")
  }
  # Extract the names of object and object2
  model1 <- deparse(substitute(object))
  model2 <- deparse(substitute(object2))
  # Extract arguments supplied in ... and determine which are named
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep_len(FALSE, length(dotargs))
  else (names(dotargs) != "")
  which_named <- which(named)
  which_not_named <- which(!named)
  # Named objects are intended for compare_models()
  for_compare_models <- dotargs[named]
  # Create list of model objects:  unnamed arguments may be model objects
  model_list <- c(list(object, object2), dotargs[!named])
  # Check for objects that do not have class "lax"
  is_chand <- vapply(model_list, function(x) inherits(x, "lax"), NA)
  if (any(!is_chand)) {
    stop("The following are not 'lax' objects: ",
         paste(names(model_list)[!is_chand], collapse = ", "))
  }
  extra_names <- as.list(substitute(list(...)))[-1][which_not_named]
  extra_names <- sapply(extra_names, function(x) deparse(x))
  model_names <- c(model1, model2, extra_names)
  # Check for duplicate names
  if (anyDuplicated(model_names)) {
    stop("A model name has been supplied more than once")
  }
  # Order the models in order of the number of parameters
  n_pars <- vapply(model_list, function(x) attr(x, "p_current"), 0)
  # Check for models with the same number of parameters
  if (anyDuplicated(n_pars)) {
    stop("At least two models have the same number of parameters")
  }
  m_order <- order(n_pars, decreasing = TRUE)
  model_list <- model_list[m_order]
  n_pars <- n_pars[m_order]
  n_models <- length(model_list)
  # Nested models: the largest model is model_list[[1]]
  alrts <- p_value <- numeric(n_models - 1)
  for (i in 2:n_models) {
    larger <- model_list[[i - 1]]
    smaller <- model_list[[i]]
    all_pars <- attr(larger, "free_pars")
    larger_names <- names(all_pars)
    smaller_names <- names(attr(model_list[[i]], "free_pars"))
    # Check that the names of all the parameters in the smaller model are
    # present in the larger model
    if (!all(smaller_names %in% larger_names)) {
      stop("parameter names are not nested")
    }
    # Which parameters have been fixed in moving from larger to smaller?
    fixed_pars <- which(!(larger_names %in% smaller_names))
    fixed_pars <- all_pars[fixed_pars]
    fixed_at <- rep(0, length(fixed_pars))
    names(fixed_at) <- names(fixed_pars)
    attr(smaller, "fixed_pars") <- fixed_pars
    attr(smaller, "fixed_at") <- fixed_at
    # In comparing smaller to larger, treat larger as the full model,
    # i.e. having no fixed parameters
    attr(larger, "fixed_pars") <- NULL
    # Do the testing
    res <- do.call(chandwich::compare_models,
                   c(list(larger = larger, smaller = smaller),
                     for_compare_models))
    alrts[i - 1] <- res$alrts
    p_value[i - 1] <- res$p_value
  }
  df <- -diff(n_pars)
  my_table <- data.frame(n_pars, c(NA, df), c(NA, alrts), c(NA, p_value))
  dimnames(my_table) <- list(model_names,
                             c("Model.Df", "Df", "ALRTS", "Pr(>ALRTS)"))
  structure(my_table, heading = c("Analysis of (Adjusted) Deviance Table\n"),
            class = c("anova", "data.frame"))
}
