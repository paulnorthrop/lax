% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lax.R
\docType{package}
\name{lax}
\alias{lax}
\title{lax: Loglikelihood Adjustment for Extreme Value Models}
\description{
Performs adjusted inferences based on model objects fitted, using maximum
likelihood estimation, by the extreme value analysis packages
\href{https://cran.r-project.org/package=evd}{evd},
\href{https://cran.r-project.org/package=evir}{evir},
\href{https://cran.r-project.org/package=extRemes}{extRemes},
\href{https://cran.r-project.org/package=fExtremes}{fExtremes},
\href{https://cran.r-project.org/package=ismev}{ismev},
\href{https://cran.r-project.org/package=mev}{mev},
\href{https://cran.r-project.org/package=POT}{POT} and
\href{https://cran.r-project.org/package=texmex}{texmex}.
Univariate extreme value models, including regression models, are supported.
Adjusted standard errors and an adjusted loglikelihood are provided,
using the
\href{https://cran.r-project.org/package=chandwich}{chandwich package}
 and the object-oriented features of the
 \href{https://cran.r-project.org/package=sandwich}{sandwich package}.
The adjustment is based on a robust sandwich estimator of the parameter
covariance matrix, based on the methodology in Chandler and Bate (2007).
This can be used for cluster correlated data
when interest lies in the parameters of the marginal distributions, or for
performing inferences that are robust to certain types of model
misspecification.  Univariate extreme value models, including regression
models, are supported.
}
\details{
Main function is \code{\link[lax]{alogLik}}, which works in an
object-oriented way, operating on fitted model objects.
This function performs the loglikelihood adjustments using
\code{\link[chandwich]{adjust_loglik}}.
See the following package-specific help pages for details and examples:
\href{https://cran.r-project.org/package=evd}{evd},
\href{https://cran.r-project.org/package=evir}{evir},
\href{https://cran.r-project.org/package=extRemes}{extRemes},
\href{https://cran.r-project.org/package=fExtremes}{fExtremes},
\href{https://cran.r-project.org/package=ismev}{ismev},
\href{https://cran.r-project.org/package=mev}{mev},
\href{https://cran.r-project.org/package=POT}{POT},
\href{https://cran.r-project.org/package=texmex}{texmex}.

See \code{vignette("lax-vignette", package = "lax")} for an overview of the
package.
}
\references{
Belzile, L., Wadsworth, J. L., Northrop, P. J., Grimshaw, S. D.
  and Huser, R. (2019). mev: Multivariate Extreme Value Distributions.
  R package version 1.12.2. \url{https://github.com/lbelzile/mev/}

Berger S., Graham N., Zeileis A. (2017). Various Versatile
  Variances: An Object-Oriented Implementation of Clustered Covariances in R.
  Technical Report 2017-12, Working Papers in Economics and Statistics,
  Research Platform Empirical and Experimental Economics, Universitat
  Innsbruck. \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2017-12}.

Chandler, R. E. and Bate, S. (2007). Inference for clustered
  data using the independence loglikelihood. \emph{Biometrika},
  \strong{94}(1), 167-183. \url{http://doi.org/10.1093/biomet/asm015}

Gilleland, E. and Katz, R. W. (2016). extRemes 2.0: An Extreme
  Value Analysis Package in R. \emph{Journal of Statistical Software},
  \strong{72}(8), 1-39. \url{http://doi.org/10.18637/jss.v072.i08}

Northrop, P. J. and Chandler, R. E. (2018).
  chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
  version 1.1. \url{https://CRAN.R-project.org/package=chandwich}.

Pfaff, B. and McNeil, A. (2018). evir: Extreme Values in R.
  R package version 1.7-4. \url{https://CRAN.R-project.org/package=evir}

Ribatet, M. and Dutang, C. (2019). POT: Generalized Pareto
  Distribution and Peaks Over Threshold. R package version 1.1-7.
  \url{https://CRAN.R-project.org/package=POT}

Southworth, H., Heffernan, J. E. and Metcalfe, P. D. (2017).
  texmex: Statistical modelling of extreme values. R package version 2.4.
  \url{https://CRAN.R-project.org/package=texmex}.

Stephenson, A. G. evd: Extreme Value Distributions.
  \emph{R News}, \strong{2}(2):31-32, June 2002.
  \url{https://CRAN.R-project.org/doc/Rnews/}

Stephenson, A. G., Heffernan, J. E. and Gilleland, E. (2018).
  ismev: An Introduction to Statistical Modeling of Extreme Values.
  R package version 1.42. \url{https://CRAN.R-project.org/package=ismev}.

Wuertz, D., Setz, T. and Chalabi, Y. (2017). fExtremes:
  Rmetrics - Modelling Extreme Events in Finance. R package version
  3042.82. \url{https://CRAN.R-project.org/package=fExtremes}

Zeileis A. (2004). Econometric Computing with HC and HAC
  Covariance Matrix Estimators. \emph{Journal of Statistical Software},
  \strong{11}(10), 1-17. \url{http://doi.org/10.18637/jss.v011.i10}.

Zeileis A. (2006). Object-Oriented Computation of Sandwich
  Estimators. \emph{Journal of Statistical Software}, \strong{16}(9),
  1-16. \url{http://doi.org/10.18637/jss.v016.i09}.
}
