% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ismev.R
\name{ismev}
\alias{ismev}
\alias{alogLik.gev.fit}
\alias{alogLik.pp.fit}
\alias{alogLik.gpd.fit}
\alias{alogLik.rlarg.fit}
\title{Loglikelihood adjustment for ismev fits}
\usage{
\method{alogLik}{gev.fit}(x, cluster = NULL, use_vcov = TRUE, ...)

\method{alogLik}{pp.fit}(x, cluster = NULL, use_vcov = TRUE, ...)

\method{alogLik}{gpd.fit}(
  x,
  cluster = NULL,
  use_vcov = TRUE,
  binom = FALSE,
  k,
  inc_cens = TRUE,
  ...
)

\method{alogLik}{rlarg.fit}(x, cluster = NULL, use_vcov = TRUE, ...)
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

\item{binom}{A logical scalar.  This option is only relevant to
\strong{GP models} and is only available in the \strong{stationary}
(no covariates) case. If \code{binom = FALSE} then loglikelihood
adjustment is only performed using the GP model. If \code{binom = TRUE}
then loglikelihood adjustment is also performed for inferences about the
probability of threshold exceedance, using a Bernoulli model for the
instances of threshold exceedance.}

\item{k}{A non-negative integer scalar.  This option is only relevant to
\strong{GP models} and is only available in the \strong{stationary}
(no covariates) case.  If \code{k} is supplied then it is passed as the
run parameter \eqn{K} to \code{\link[exdex]{kgaps}} for making inferences
about the extremal index \eqn{\theta} using the \eqn{K}-gaps model of
Suveges and Davison (2010).}

\item{inc_cens}{A logical scalar.  This argument is only relevant if
\code{k} is supplied.  Passed to \code{\link[exdex]{kgaps}} to indicate
whether or not to include censored inter-exceedance times, relating to
the first and last observations.}
}
\value{
An object inheriting from class \code{"chandwich"}.  See
  \code{\link[chandwich]{adjust_loglik}}.

  \code{class(x)} is a vector of length 5. The first 3 components are
  \code{c("lax", "chandwich", "ismev")}.
  The remaining 2 components depend on the model that was fitted.
  The 4th component is:
  \code{"gev"} if \code{\link[ismev]{gev.fit}}
  (or \code{\link{gev_refit}}) was used;
  \code{"gpd"} if \code{\link[ismev]{gpd.fit}}
  (or \code{\link{gpd_refit}}) was used;
  \code{"pp"} \code{\link[ismev]{pp.fit}}
  (or \code{\link{pp_refit}}) was used;
  \code{"rlarg"} \code{\link[ismev]{rlarg.fit}}
  (or \code{\link{rlarg_refit}}) was used.
  The 5th component is
  \code{"stat"} if \code{x$trans = FALSE} and
  \code{"nonstat"} if \code{x$trans = TRUE}.
}
\description{
S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions
\code{\link[ismev]{gev.fit}}, \code{\link[ismev]{gpd.fit}},
\code{\link[ismev]{pp.fit}} and \code{\link[ismev]{rlarg.fit}} in the
\code{\link[ismev]{ismev}} package.  If regression modelling is used then
the model will need to be re-fitted, see \code{\link{ismev_refits}}.
}
\details{
See \code{\link{alogLik}} for details.

If regression modelling is used then the ismev functions
\code{\link[ismev]{gev.fit}}, \code{\link[ismev]{gpd.fit}},
\code{\link[ismev]{pp.fit}} and \code{\link[ismev]{rlarg.fit}}
return residuals but \code{\link{alogLik}} needs the raw data.
The model will need to be re-fitted, using one of the functions in
\code{\link{ismev_refits}}, and the user will be prompted to do this
by an error message produced by \code{\link{alogLik}}.
}
\examples{
# We need the ismev package
got_ismev <- requireNamespace("ismev", quietly = TRUE)

if (got_ismev) {
  library(ismev)

  # GEV model -----

  # An example from the ismev::gev.fit documentation
  gev_fit <- gev.fit(revdbayes::portpirie, show = FALSE)
  adj_gev_fit <- alogLik(gev_fit)
  summary(adj_gev_fit)

  # An example from chapter 6 of Coles (2001)
  data(fremantle)
  xdat <- fremantle[, "SeaLevel"]
  # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
  ydat <- cbind(fremantle[, "Year"] - 1896, fremantle[, "SOI"])
  gev_fit <- gev_refit(xdat, ydat, mul = 1:2, show = FALSE)
  adj_gev_fit <- alogLik(gev_fit)
  summary(adj_gev_fit)

  # An example from Chandler and Bate (2007)
  gev_fit <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4,
                       show = FALSE)
  adj_gev_fit <- alogLik(gev_fit, cluster = ow$year)
  summary(adj_gev_fit)
  # Get closer to the values reported in Table 2 of Chandler and Bate (2007)
  gev_fit <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4,
                       show = FALSE, method = "BFGS")
  # Call sandwich::meatCL() with cadjust = FALSE
  adj_gev_fit <- alogLik(gev_fit, cluster = ow$year, cadjust = FALSE)
  summary(adj_gev_fit)

  # GP model -----

  # An example from the ismev::gpd.fit documentation
  \donttest{
  data(rain)
  rain_fit <- gpd.fit(rain, 10, show = FALSE)
  adj_rain_fit <- alogLik(rain_fit)
  summary(adj_rain_fit)
  # Continuing to the regression example on page 119 of Coles (2001)
  ydat <- as.matrix((1:length(rain)) / length(rain))
  reg_rain_fit <- gpd_refit(rain, 30, ydat = ydat, sigl = 1, siglink = exp,
                            show = FALSE)
  adj_reg_rain_fit <- alogLik(reg_rain_fit)
  summary(adj_reg_rain_fit)
  }
  # Binomial-GP model -----

  # Use Newlyn seas surges data from the exdex package
  surges <- exdex::newlyn
  u <- quantile(surges, probs = 0.9)
  newlyn_fit <- gpd.fit(surges, u, show = FALSE)
  # Create 5 clusters each corresponding approximately to 1 year of data
  cluster <- rep(1:5, each = 579)[-1]
  adj_newlyn_fit <- alogLik(newlyn_fit, cluster = cluster, binom = TRUE,
                            cadjust = FALSE)
  summary(adj_newlyn_fit)
  summary(attr(adj_newlyn_fit, "pu_aloglik"))

  # Add inference about the extremal index theta, using K = 1
  adj_newlyn_theta <- alogLik(newlyn_fit, cluster = cluster, binom = TRUE,
                              k = 1, cadjust = FALSE)
  summary(attr(adj_newlyn_theta, "theta"))

  # PP model -----

  # An example from the ismev::pp.fit documentation
  data(rain)
  # Start from the mle to save time
  init <- c(40.55755732, 8.99195409, 0.05088103)
  muinit <- init[1]
  siginit <- init[2]
  shinit <- init[3]
  rain_fit <- pp_refit(rain, 10, muinit = muinit, siginit = siginit,
                       shinit = shinit, show = FALSE)
  adj_rain_fit <- alogLik(rain_fit)
  summary(adj_rain_fit)

  # An example from chapter 7 of Coles (2001).
  # Code from demo ismev::wooster.temps
  data(wooster)
  x <- seq(along = wooster)
  usin <- function(x, a, b, d) {
    return(a + b * sin(((x - d) * 2 * pi) / 365.25))
  }
  wu <- usin(x, -30, 25, -75)
  ydat <- cbind(sin(2 * pi * x / 365.25), cos(2 * pi *x / 365.25))
  # Start from the mle to save time
  init <- c(-15.3454188, 9.6001844, 28.5493828, 0.5067104, 0.1023488,
            0.5129783, -0.3504231)
  muinit <- init[1:3]
  siginit <- init[4:6]
  shinit <- init[7]
  wooster.pp <- pp_refit(-wooster, threshold = wu, ydat = ydat, mul = 1:2,
                         sigl = 1:2, siglink = exp, method = "BFGS",
                         muinit = muinit, siginit = siginit, shinit = shinit,
                         show = FALSE)
  adj_pp_fit <- alogLik(wooster.pp)
  summary(adj_pp_fit)

  # r-largest order statistics model -----

  # An example based on the ismev::rlarg.fit() documentation
  vdata <- revdbayes::venice
  rfit <- rlarg.fit(vdata, muinit = 120.54, siginit = 12.78,
                    shinit = -0.1129, show = FALSE)
  adj_rfit <- alogLik(rfit)
  summary(adj_rfit)

  \donttest{
  # Adapt this example to add a covariate
  set.seed(30102019)
  ydat <- matrix(runif(nrow(vdata)), nrow(vdata), 1)
  rfit2 <- rlarg_refit(vdata, ydat = ydat, mul = 1,
                       muinit = c(120.54, 0), siginit = 12.78,
                       shinit = -0.1129, show = FALSE)
  adj_rfit2 <- alogLik(rfit2)
  summary(adj_rfit2)
  }
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
