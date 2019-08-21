
<!-- README.md is generated from README.Rmd. Please edit that file -->
lax
===

[![Travis-CI Build Status](https://travis-ci.org/paulnorthrop/lax.svg?branch=master)](https://travis-ci.org/paulnorthrop/lax) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/lax?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/lax) [![Coverage Status](https://codecov.io/github/paulnorthrop/lax/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/lax?branch=master)

Loglikelihood Adjustment for Extreme Value Models
-------------------------------------------------

### What does lax do?

The [CRAN Task View on Extreme Value Analysis](https://CRAN.R-project.org/view=ExtremeValue) provides information about R packages that perform various extreme value analyses. The *lax* package supplements the functionality of 7 of these packages, namely [evd](https://cran.r-project.org/package=evd), [evir](https://cran.r-project.org/package=evir), [extRemes](https://cran.r-project.org/package=extRemes), [fExtremes](https://cran.r-project.org/package=fExtremes), [ismev](https://cran.r-project.org/package=ismev), [POT](https://cran.r-project.org/package=POT) and [texmex](https://cran.r-project.org/package=texmex).

The [chandwich](https://cran.r-project.org/package=chandwich) package is used to provide robust sandwich estimation of parameter covariance matrix and loglikelihood adjustment for models fitted by maximum likelihood estimation. This may be useful for cluster correlated data when interest lies in the parameters of the marginal distributions, or for performing inferences that are robust to certain types of model misspecification.

*lax* works in an object-oriented way, operating on R objects returned from functions in other packages that summarise the fit of an extreme value model. Univariate extreme value models, including regression models, are supported. Loglikelihood adjustment and sandwich estimation is performed by an `alogLik` S3 method, illustrated by the following example.

### A simple example

This example is based on the analysis presented in Section 5.2 of [Chandler and Bate (2007)](https://doi.org/10.1093/biomet/asm015). The data are a bivariate time series of annual maximum temperatures, recorded in degress Fahrenheit, at Oxford and Worthing in England, for the period 1901 to 1980. If interest is only in the marginal distributions of high temperatures in Oxford and Worthing, then we might fit a GEV regression model in which some or all of the parameters may vary between Oxford and Worthing. However, we should adjust for the cluster dependence between temperatures recorded during the same year.

The following code fits such a model, allowing only the GEV location parameter to vary. The model is fitted using the `fgev` function from the [evd](https://cran.r-project.org/package=evd) package. Then `alogLik` is used to provide adjusted standard errors and an adjusted loglikelihood. Finally, a `confint` method is used to calculate approximate 95% confidence intervals for the parameters, based on the adjusted loglikelihood.

``` r
  library(lax)
  y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
  x <- rep(c(1, -1), each = length(y) / 2)
  # Fit a GEV model with separate location parameters for Oxford and Worthing
  owfit <- evd::fgev(y, nsloc = x)
  # Create a vector to 
  year <- rep(rownames(chandwich::owtemps), 2)
  adj_owfit <- alogLik(owfit, cluster = year)
  summary(adj_owfit)
#>              MLE      SE adj. SE
#> loc      81.1600 0.32980 0.40830
#> loctrend  2.5040 0.31080 0.19910
#> scale     3.7900 0.22810 0.25230
#> shape    -0.2097 0.04765 0.04063
  confint(adj_owfit)
#> Waiting for profiling to be done...
#>               2.5 %     97.5 %
#> loc      80.3518858 81.9584587
#> loctrend  2.0989464  2.8841321
#> scale     3.3449036  4.3451673
#> shape    -0.2870718 -0.1231754
```

The analysis in Chandler and Bate (2007) allowed all three GEV parameters to vary between Oxford and Worthing. The `fgev` function cannot fit this model but the function `gev.fit` in the [ismev](https://cran.r-project.org/package=ismev) package can. However, an object returned from `gev.fit` does not provide all the information about a fitted regression model that `alogLik` needs. Therefore, *lax* provides the function `gev_refit`, which is a version of `gev.fit` that adds this information. The following code reproduces the values in rows 1, 3 and 4 in Table 2 of Chandler and Bate (2007).

``` r
x <- as.matrix(rep(c(1, -1), each = length(y) / 2))
gev_fit <- gev_refit(y, x, mul = 1, sigl = 1, shl = 1, show = FALSE,
                     method = "BFGS")
adj_gev_fit <- alogLik(gev_fit, cluster = year, cadjust = FALSE)
summary(adj_gev_fit)
#>             MLE      SE adj. SE
#> loc    81.17000 0.32820 0.40360
#> loc1    2.66900 0.32820 0.21280
#> scale   3.72900 0.22920 0.24250
#> scale1  0.53120 0.22920 0.19110
#> shape  -0.19890 0.04937 0.03944
#> shape1 -0.08835 0.04937 0.03625
```

The estimation of the 'meat' of the sandwich adjustment is performed using the [sandwich](https://cran.r-project.org/package=sandwich) package. In this example, we need to pass `cadjust = FALSE` to `sandwich::meatCL` in order that the adjustment is the same as that used in Chandler and Bate (2007).

### Installation

To get the current released version from CRAN:

``` r
install.packages("lax")
```

### Vignettes
