#' lax: Loglikelihood Adjustment for Extreme Value Models
#'
#' Performs adjusted inferences based on model objects fitted, using maximum
#' likelihood estimation, by the extreme value analysis packages
#' \href{https://cran.r-project.org/package=evd}{evd},
#' \href{https://cran.r-project.org/package=evir}{evir},
#' \href{https://cran.r-project.org/package=extRemes}{extRemes},
#' \href{https://cran.r-project.org/package=fExtremes}{fExtremes},
#' \href{https://cran.r-project.org/package=ismev}{ismev},
#' \href{https://cran.r-project.org/package=POT}{POT} and
#' \href{https://cran.r-project.org/package=texmex}{texmex}.
#' Univariate extreme value models, including regression models, are supported.
#' Adjusted standard errors and an adjusted loglikelihood are provided,
#' using the
#' \href{https://cran.r-project.org/package=chandwich}{chandwich package}
#'  and the object-oriented features of the
#'  \href{https://cran.r-project.org/package=sandwich}{sandwich package}.
#' The adjustment is based on a robust sandwich estimator of the parameter
#' covariance matrix, based on the methodology in Chandler and Bate (2007).
#' This can be used for cluster correlated data
#' when interest lies in the parameters of the marginal distributions, or for
#' performing inferences that are robust to certain types of model
#' misspecification.  Univariate extreme value models, including regression
#' models, are supported.
#' @details Main function is \code{\link[lax]{alogLik}}, which works in an
#' object-oriented way, operating on fitted model objects.
#' This function performs the loglikelihood adjustments using
#' \code{\link[chandwich]{adjust_loglik}}.
#' See the following package-specific help pages for details and examples:
#' \href{https://cran.r-project.org/package=evd}{evd},
#' \href{https://cran.r-project.org/package=evir}{evir},
#' \href{https://cran.r-project.org/package=extRemes}{extRemes},
#' \href{https://cran.r-project.org/package=fExtremes}{fExtremes},
#' \href{https://cran.r-project.org/package=ismev}{ismev},
#' \href{https://cran.r-project.org/package=POT}{POT},
#' \href{https://cran.r-project.org/package=texmex}{texmex}.
#'
#' See \code{vignette("lax-vignette", package = "lax")} for an overview of the
#' package.
#' @references Berger S., Graham N., Zeileis A. (2017). Various Versatile
#'   Variances: An Object-Oriented Implementation of Clustered Covariances in R.
#'   Technical Report 2017-12, Working Papers in Economics and Statistics,
#'   Research Platform Empirical and Experimental Economics, Universitat
#'   Innsbruck. \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2017-12}.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://doi.org/10.1093/biomet/asm015}
#' @references Gilleland, E. and Katz, R. W. (2016). extRemes 2.0: An Extreme
#'   Value Analysis Package in R. \emph{Journal of Statistical Software},
#'   \strong{72}(8), 1-39. \url{http://doi.org/10.18637/jss.v072.i08}
#' @references Northrop, P. J. and Chandler, R. E. (2018).
#'   chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
#'   version 1.1. \url{https://CRAN.R-project.org/package=chandwich}.
#' @references Pfaff, B. and McNeil, A. (2018). evir: Extreme Values in R.
#'   R package version 1.7-4. \url{https://CRAN.R-project.org/package=evir}
#' @references Ribatet, M. and Dutang, C. (2019). POT: Generalized Pareto
#'   Distribution and Peaks Over Threshold. R package version 1.1-7.
#'   \url{https://CRAN.R-project.org/package=POT}
#' @references Southworth, H., Heffernan, J. E. and Metcalfe, P. D. (2017).
#'   texmex: Statistical modelling of extreme values. R package version 2.4.
#'   \url{https://CRAN.R-project.org/package=texmex}.
#' @references Stephenson, A. G. evd: Extreme Value Distributions.
#'   \emph{R News}, \strong{2}(2):31-32, June 2002.
#'   \url{https://CRAN.R-project.org/doc/Rnews/}
#' @references Stephenson, A. G., Heffernan, J. E. and Gilleland, E. (2018).
#'   ismev: An Introduction to Statistical Modeling of Extreme Values.
#'   R package version 1.42. \url{https://CRAN.R-project.org/package=ismev}.
#' @references Wuertz, D., Setz, T. and Chalabi, Y. (2017). fExtremes:
#'   Rmetrics - Modelling Extreme Events in Finance. R package version
#'   3042.82. \url{https://CRAN.R-project.org/package=fExtremes}
#' @references Zeileis A. (2004). Econometric Computing with HC and HAC
#'   Covariance Matrix Estimators. \emph{Journal of Statistical Software},
#'   \strong{11}(10), 1-17. \url{http://doi.org/10.18637/jss.v011.i10}.
#' @references Zeileis A. (2006). Object-Oriented Computation of Sandwich
#'   Estimators. \emph{Journal of Statistical Software}, \strong{16}(9),
#'   1-16. \url{http://doi.org/10.18637/jss.v016.i09}.
#' @docType package
#' @name lax
#' @import methods
#' @import sandwich
#' @importFrom stats nobs vcov coef logLik
#' @importFrom graphics plot
NULL

#' Oxford and Worthing annual maximum temperatures
#'
#' Annual maximum temperatures at Oxford and Worthing (England), for the
#' period 1901 to 1980.
#'
#' @format A dataframe with 80 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{temp}: }{annual maximum temperatures in degrees
#'       Fahrenheit.}
#'     \item{Column 2, \code{year}: }{year in which the maximum was recorded.}
#'     \item{Column 3, \code{name}: }{name of location, "oxford" or "worthing"}
#'     \item{Column 4, \code{loc}: }{location: 1 for "oxford", -1 for
#'       "worthing"}
#'  }
#' @source Tabony, R. C. (1983) Extreme value analysis in meteorology.
#'  \emph{The Meteorological Magazine}, \strong{112}, 77-98.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
"ow"
