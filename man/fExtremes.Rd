% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fExtremes.R
\name{fExtremes}
\alias{fExtremes}
\alias{alogLik.fGEVFIT}
\alias{alogLik.fGPDFIT}
\title{Loglikelihood adjustment for fExtremes fits}
\usage{
\method{alogLik}{fGEVFIT}(x, cluster = NULL, use_vcov = TRUE, ...)

\method{alogLik}{fGPDFIT}(x, cluster = NULL, use_vcov = TRUE, ...)
}
\arguments{
\item{x}{A fitted model object with certain associated S3 methods.
See \strong{Details}.}

\item{cluster}{A vector or factor indicating from which cluster the
  respective log-likelihood contributions from \code{loglik} originate.
  The length of \code{cluster} must be consistent with the \code{estfun}
  method to be used in the estimation of the 'meat' \code{V} of the sandwich
  estimator of the covariance matrix of the parameters to be passed to
  \code{\link[chandwich]{adjust_loglik}}.  In most cases, \code{cluster}
  must have length equal to the number of observations in data.  The
  exception is the GP (only) model (\code{binom = FALSE}), where the
  \code{cluster} may either contain a value for each observation in the raw
  data, or for each threshold exceedance in the data.

  If \code{cluster} is not supplied (is \code{NULL}) then it is
  assumed that each observation forms its own cluster.
  See \strong{Details} for further details.}

\item{use_vcov}{A logical scalar.  Should we use the \code{vcov} S3 method
for \code{x} (if this exists) to estimate the Hessian of the independence
loglikelihood to be passed as the argument \code{H} to
\code{\link[chandwich]{adjust_loglik}}?
Otherwise, \code{H} is estimated inside
\code{\link[chandwich]{adjust_loglik}} using
\code{\link[stats:optim]{optimHess}}.}

\item{...}{Further arguments to be passed to the functions in the
sandwich package \code{\link[sandwich]{meat}} (if \code{cluster = NULL}),
or \code{\link[sandwich:vcovCL]{meatCL}} (if \code{cluster} is not
\code{NULL}).}
}
\value{
An object inheriting from class \code{"chandwich"}.  See
  \code{\link[chandwich]{adjust_loglik}}.
  \code{class(x)} is a vector of length 5. The first 3 components are
  \code{c("lax", "chandwich", "fExtremes")}.
  The remaining 2 components depend on the model that was fitted.
  If \code{\link[fExtremes:GevModelling]{gevFit}} or
  \code{\link[fExtremes:GevModelling]{gumbelFit}} was used then these
  components are \code{c("gev", "stat")}.
  If \code{\link[fExtremes:GpdModelling]{gpdFit}} was used then these
  components are \code{c("gpd", "stat")}.
}
\description{
S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions
\code{\link[fExtremes:GevModelling]{gevFit}},
\code{\link[fExtremes:GevModelling]{gumbelFit}} and
\code{\link[fExtremes:GpdModelling]{gpdFit}}
in the \code{\link[fExtremes:00Extremes-package]{fExtremes}} package.
The model must have been fitted using maximum likelihood estimation.
}
\details{
See \code{\link{alogLik}} for details.
}
\examples{
# We need the fExtremes package
got_fExtremes <- requireNamespace("fExtremes", quietly = TRUE)
if (got_fExtremes) {
  library(fExtremes)

  # GEV
  # An example from the fExtremes::gevFit documentation
  set.seed(4082019)
  x <- fExtremes::gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
  # Fit GEV distribution by maximum likelihood estimation
  fit <- fExtremes::gevFit(x)
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # GP
  # An example from the fExtremes::gpdFit documentation
  # Simulate GP data
  x <- fExtremes::gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
  # Fit GP distribution by maximum likelihood estimation
  fit <- fExtremes::gpdFit(x, u = min(x))
  adj_fit <- alogLik(fit)
  summary(adj_fit)
}
}
\references{
Chandler, R. E. and Bate, S. (2007). Inference for clustered
  data using the independence loglikelihood. \emph{Biometrika},
  \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}

Suveges, M. and Davison, A. C. (2010) Model
  misspecification in peaks over threshold analysis, \emph{The Annals of
  Applied Statistics}, \strong{4}(1), 203-221.
  \doi{10.1214/09-AOAS292}

Zeileis (2006) Object-Oriented Computation and Sandwich
  Estimators.  \emph{Journal of Statistical Software}, \strong{16}, 1-16.
  \doi{10.18637/jss.v016.i09}
}
\seealso{
\code{\link{alogLik}}: loglikelihood adjustment for model fits.
}
