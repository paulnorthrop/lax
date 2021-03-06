% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extRemes.R
\name{extRemes}
\alias{extRemes}
\alias{alogLik.fevd}
\title{Loglikelihood adjustment for extRemes fits}
\usage{
\method{alogLik}{fevd}(x, cluster = NULL, use_vcov = TRUE, ...)
}
\arguments{
\item{x}{A fitted model object with certain associated S3 methods.
See \strong{Details}.}

\item{cluster}{A vector or factor indicating from which cluster the
respective loglikelihood contributions from \code{loglik} originate.
This must have the same length as the vector returned by the
\code{logLikVec} method for an object like \code{x}.
If \code{cluster} is not supplied (i.e. is \code{NULL}) then it is
assumed that each observation forms its own cluster.
See \strong{Details}.}

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
  \code{c("lax", "chandwich", "extRemes")}.
  The remaining 2 components depend on the model that was fitted.
  The 4th component is: \code{"gev"} if \code{x$type = "GEV"} or
  \code{x$type = "Gumbel"}; \code{"gp"} if \code{x$type = "GP"} or
  \code{x$type = "Exponential"}; \code{"pp"} if \code{x$type = "PP"}.
  The 5th component is
  \code{"stat"} if \code{\link[extRemes]{is.fixedfevd} = TRUE} and
  \code{"nonstat"} if \code{\link[extRemes]{is.fixedfevd} = FALSE}.
}
\description{
S3 \code{alogLik} method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the function
\code{\link[extRemes]{fevd}} in the
\code{\link[extRemes:extRemes-package]{extRemes}} package.
The model must have been fitted using maximum likelihood estimation.
}
\details{
See \code{\link{alogLik}} for details.
}
\examples{
# We need the extRemes and distillery packages
got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
got_distillery <- requireNamespace("distillery", quietly = TRUE)

if (got_extRemes & got_distillery) {
  library(extRemes)
  library(distillery)
  # Examples from the extRemes::fevd documentation
  data(PORTw)

  # GEV
  fit0 <- fevd(TMX1, PORTw, units = "deg C", use.phi = TRUE)
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GEV regression
  fitPORTstdmax <- fevd(TMX1, PORTw, scale.fun = ~STDTMAX, use.phi = TRUE)
  adj_fit1 <- alogLik(fitPORTstdmax)
  summary(adj_fit1)
  fitPORTstdmax2 <- fevd(TMX1, PORTw, location.fun = ~STDTMAX,
                         scale.fun = ~STDTMAX, use.phi = TRUE)
  adj_fit2 <- alogLik(fitPORTstdmax2)
  summary(adj_fit2)
  anova(adj_fit0, adj_fit1)
  anova(adj_fit1, adj_fit2)
  anova(adj_fit0, adj_fit2)
  anova(adj_fit0, adj_fit1, adj_fit2)

  # Gumbel
  fit0 <- fevd(TMX1, PORTw, type = "Gumbel", units = "deg C")
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GP
  data(damage)
  fit1 <- fevd(Dam, damage, threshold = 6, type = "GP",
               time.units = "2.05/year")
  adj_fit1 <- alogLik(fit1)
  summary(adj_fit1)

  # Exponential
  fit0 <- fevd(Dam, damage, threshold = 6, type="Exponential",
               time.units = "2.05/year")
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GP non-constant threshold
  data(Fort)
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
              type = "GP")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # Exponential non-constant threshold
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
              type = "Exponential")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # PP model
  fit <- fevd(Prec, Fort, threshold = 0.475, type = "PP", units = "inches")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # PP non-constant threshold
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
              type = "PP")
  adj_fit <- alogLik(fit)
  summary(adj_fit)
}
}
\references{
Chandler, R. E. and Bate, S. (2007). Inference for clustered
  data using the independence loglikelihood. \emph{Biometrika},
  \strong{94}(1), 167-183. \url{http://doi.org/10.1093/biomet/asm015}

Zeleis (2006) Object-Oriented Computation and Sandwich
  Estimators.  \emph{Journal of Statistical Software}, \strong{16}, 1-16.
  \url{http://doi.org/10.18637/jss.v016.i09}
}
\seealso{
\code{\link{alogLik}}: loglikelihood adjustment for model fits.
}
