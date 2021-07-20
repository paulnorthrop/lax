
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lax

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/lax?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/lax)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/lax/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/lax?branch=master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/lax)](https://cran.r-project.org/package=lax)

## Loglikelihood Adjustment for Extreme Value Models

### What does lax do?

The [CRAN Task View on Extreme Value
Analysis](https://CRAN.R-project.org/view=ExtremeValue) provides
information about R packages that perform various extreme value
analyses. The *lax* package supplements the univariate extreme value
modelling, including regression modelling, provided by 9 of these
packages, namely [eva](https://cran.r-project.org/package=eva),
[evd](https://cran.r-project.org/package=evd),
[evir](https://cran.r-project.org/package=evir),
[extRemes](https://cran.r-project.org/package=extRemes),
[fExtremes](https://cran.r-project.org/package=fExtremes),
[ismev](https://cran.r-project.org/package=ismev),
[mev](https://cran.r-project.org/package=mev),
[POT](https://cran.r-project.org/package=POT) and
[texmex](https://cran.r-project.org/package=texmex). *lax* works in an
object-oriented way, operating on R objects returned from functions in
other packages that summarise the fit of an extreme value model. It uses
the [chandwich](https://cran.r-project.org/package=chandwich) package to
provide robust sandwich estimation of parameter covariance matrix and
loglikelihood adjustment for models fitted by maximum likelihood
estimation. This is performed by an `alogLik` S3 method, illustrated by
the following example.

### An example

This example is based on the analysis presented in Section 5.2 of
[Chandler and Bate (2007)](https://doi.org/10.1093/biomet/asm015). The
data, which are available in the data frame `ow`, are a bivariate time
series of annual maximum temperatures, recorded in degrees Fahrenheit,
at Oxford and Worthing in England, for the period 1901 to 1980. If
interest is only in the marginal distributions of high temperatures in
Oxford and Worthing, then we might fit a GEV regression model in which
some or all of the parameters may vary between Oxford and Worthing.
However, we should adjust for the cluster dependence between
temperatures recorded during the same year.

The following code fits such a model using the `evm` function in the
[texmex](https://cran.r-project.org/package=texmex) package and the uses
`alogLik` to perform adjusted inferences.

``` r
library(lax)
library(texmex, quietly = TRUE)
# Fit a GEV model with separate location, scale and shape for Oxford and Worthing
# Note: phi = log(scale)
evm_fit <- evm(temp, ow, gev, mu = ~ loc, phi = ~ loc, xi = ~loc)
# Adjust the loglikelihood and standard errors
adj_evm_fit <- alogLik(evm_fit, cluster = ow$year, cadjust = FALSE)
# MLEs, SEs and adjusted SEs
summary(adj_evm_fit)
#>                       MLE      SE adj. SE
#> mu: (Intercept)  81.17000 0.32820 0.40360
#> mu: loc           2.66800 0.32820 0.21280
#> phi: (Intercept)  1.30600 0.06091 0.06490
#> phi: loc          0.14330 0.06091 0.05074
#> xi: (Intercept)  -0.19900 0.04937 0.03943
#> xi: loc          -0.08821 0.04937 0.03624
```

An object returned from `aloglik` is a function to evaluate the adjusted
loglikelihood, with `anova`, `coef`, `confint`, `logLik`, `nobs`,
`plot`, `print`, `summary` and `vcov` methods.

### Installation

To get the current released version from CRAN:

``` r
install.packages("lax")
```

### Vignette

See `vignette("lax-vignette", package = "lax")` for an overview of the
package.
