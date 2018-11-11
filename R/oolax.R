#' oolax: Objected-orientated Loglikelihood Adjustment for Extreme Value Models
#'
#' Performs adjusted inferences using fitted model objects returned from the
#' extreme value analysis packages
#' \href{https://cran.r-project.org/package=evd}{evd}
#' \href{https://cran.r-project.org/package=evir}{evir}
#' \href{https://cran.r-project.org/package=ismev}{ismev}
#' \href{https://cran.r-project.org/package=texmex}{texmex}.
#' Adjusted standard errors and an adjusted loglikelihood are provided,
#' using the
#' \href{https://cran.r-project.org/package=chandwich}{chandwich package}
#'  and the object-orientated features of the
#'  \href{https://cran.r-project.org/package=sandwich}{sandwich package}.
#' Adjusted standard errors and an adjusted loglikelihood are provided.
#' The adjustment is based on a robust sandwich estimator of the parameter
#' covariance matrix, based on the methodology in Chandler and Bate (2007)
#' <doi:10.1093/biomet/asm015>. This can be used for cluster correlated data
#' when interest lies in the parameters of the marginal distributions or for
#' performing inferences that are robust to certain types of model
#' misspecification.  Currently only univariate extreme value models are
#' supported.  As a byproduct, generic S3 logLik, coef, vcov and nobs methods
#' are provided, where these are not already provided in the original package.
#' @details
#' \code{\link[oolax]{evd}},
#' \code{\link[oolax]{evir}}.
#' \code{\link[oolax]{ismev}},
#' \code{\link[oolax]{texmex}}.
#'
#' See Chandler and Bate (2007) for full details and
#' \code{vignette("oola-vignette", package = "oola")} for an
#' overview of the package.
#' @references Berger S., Graham N., Zeileis A. (2017). Various Versatile
#'   Variances: An Object-Oriented Implementation of Clustered Covariances in R.
#'   Technical Report 2017-12, Working Papers in Economics and Statistics,
#'   Research Platform Empirical and Experimental Economics, Universitat
#'   Innsbruck. \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2017-12}.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @references Northrop, P. J. and Chandler, R. E. (2018).
#'   chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
#'   version 1.1. \url{https://CRAN.R-project.org/package=chandwich}.
#' @references Pfaff, B. and McNeil, A. (2018). evir: Extreme Values in R.
#'   R package version 1.7-4. \url{https://CRAN.R-project.org/package=evir}
#' @references Southworth, H., Heffernan, J. E. and Metcalfe, P. D. (2017).
#'   texmex: Statistical modelling of extreme values. R package version 2.4.
#'   \url{https://CRAN.R-project.org/package=texmex}.
#' @references Stephenson, A. G. evd: Extreme Value Distributions.
#'   \emph{R News}, \strong{2}(2):31-32, June 2002.
#'   \url{ttps://CRAN.R-project.org/doc/Rnews/}
#' @references Stephenson, A. G., Heffernan, J. E. and Gilleland, E. (2018).
#'   ismev: An Introduction to Statistical Modeling of Extreme Values.
#'   R package version 1.42. \url{https://CRAN.R-project.org/package=ismev}.
#' @references Zeileis A. (2004). Econometric Computing with HC and HAC
#'   Covariance Matrix Estimators. \emph{Journal of Statistical Software},
#'   \strong{11}(10), 1-17. \url{http://doi.org/10.18637/jss.v011.i10}.
#' @references Zeileis A. (2006). Object-Oriented Computation of Sandwich
#'   Estimators. \emph{Journal of Statistical Software}, \strong{16}(9),
#'   1-16. \url{http://doi.org/10.18637/jss.v016.i09}.
#' @docType package
#' @name oolax
#' @import methods
#' @import sandwich
NULL
