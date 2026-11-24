# Loglikelihood adjustment of texmex fits

S3 `alogLik` method to perform loglikelihood adjustment of fitted
extreme value model objects returned from the `evm` function in the
`texmex` package. The model must have been fitted using maximum
likelihood estimation.

## Usage

``` r
# S3 method for class 'evmOpt'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)
```

## Arguments

- x:

  A fitted model object with certain associated S3 methods. See
  **Details**.

- cluster:

  A vector or factor indicating from which cluster the respective
  log-likelihood contributions from `loglik` originate. The length of
  `cluster` must be consistent with the `estfun` method to be used in
  the estimation of the 'meat' `V` of the sandwich estimator of the
  covariance matrix of the parameters to be passed to
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).
  In most cases, `cluster` must have length equal to the number of
  observations in data. The exception is the GP (only) model
  (`binom = FALSE`), where the `cluster` may either contain a value for
  each observation in the raw data, or for each threshold exceedance in
  the data.

  If `cluster` is not supplied (is `NULL`) then it is assumed that each
  observation forms its own cluster. See **Details** for further
  details.

- use_vcov:

  A logical scalar. Should we use the `vcov` S3 method for `x` (if this
  exists) to estimate the Hessian of the independence loglikelihood to
  be passed as the argument `H` to
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)?
  Otherwise, `H` is estimated inside
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
  using [`optimHess`](https://rdrr.io/r/stats/optim.html).

- ...:

  Further arguments to be passed to the functions in the sandwich
  package
  [`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html)
  (if `cluster = NULL`), or
  [`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  (if `cluster` is not `NULL`).

## Value

An object inheriting from class `"chandwich"`. See
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).
`class(x)` is a vector of length 5. The first 3 components are
`c("lax", "chandwich", "texmex")`. The remaining 2 components depend on
the model that was fitted. The 4th component is: `"gev"` if
`x$family$name = "GEV"`; `"gpd"` if `x$family$name = "GPD"`; `"egp3"` if
`x$family$name = "EGP3"`. The 5th component is `"stat"` if there are no
covariates in the mode and `"nonstat"` otherwise.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
for details.

## References

Chandler, R. E. and Bate, S. (2007). Inference for clustered data using
the independence loglikelihood. *Biometrika*, **94**(1), 167-183.
[doi:10.1093/biomet/asm015](https://doi.org/10.1093/biomet/asm015)

Suveges, M. and Davison, A. C. (2010) Model misspecification in peaks
over threshold analysis, *The Annals of Applied Statistics*, **4**(1),
203-221. [doi:10.1214/09-AOAS292](https://doi.org/10.1214/09-AOAS292)

Zeileis (2006) Object-Oriented Computation and Sandwich Estimators.
*Journal of Statistical Software*, **16**, 1-16.
[doi:10.18637/jss.v016.i09](https://doi.org/10.18637/jss.v016.i09)

## See also

[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md):
loglikelihood adjustment for model fits.

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run to avoid a CRAN check error inherited from the texmex package
# We need the texmex package, and ismev for the fremantle dataset
got_texmex <- requireNamespace("texmex", quietly = TRUE)
got_ismev <- requireNamespace("ismev", quietly = TRUE)
if (got_texmex) {
  library(texmex)
  # Examples from the texmex::evm documentation

  # GEV
  mod <- evm(SeaLevel, data = texmex::portpirie, family = gev)
  adj_mod <- alogLik(mod)
  summary(adj_mod)

  # GP
  mod <- evm(rain, th = 30)
  adj_mod <- alogLik(mod)
  summary(adj_mod)
  mod <- evm(rain, th = 30, cov = "sandwich")
  mod$se
  vcov(adj_mod)
  vcov(mod)

  # EGP3
  mod <- evm(rain, th = 30, family = egp3)
  adj_mod <- alogLik(mod)
  summary(adj_mod)

  # GP regression
  # An example from page 119 of Coles (2001)
  n_rain <- length(rain)
  rain_df <- data.frame(rain = rain, time = 1:n_rain / n_rain)
  evm_fit <- evm(y = rain, data = rain_df, family = gpd, th = 30,
                 phi = ~ time)
  adj_evm_fit <- alogLik(evm_fit)
  summary(adj_evm_fit)
  evm_fit <- evm(y = rain, data = rain_df, family = gpd, th = 30,
                 phi = ~ time, cov = "sandwich")
  evm_fit$se
  vcov(adj_evm_fit)
  vcov(evm_fit)

  # GEV regression
  # An example from page 113 of Coles (2001)
  if (got_ismev) {
    library(ismev)
    data(fremantle)
    new_fremantle <- fremantle
    # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
    new_fremantle[, "Year"] <- new_fremantle[, "Year"] - 1896
    evm_fit <- evm(y = SeaLevel, data = new_fremantle, family = gev,
                   mu = ~ Year + SOI)
    adj_evm_fit <- alogLik(evm_fit)
    summary(adj_evm_fit)
  }

  # An example from Chandler and Bate (2007)
  # Note: evm uses phi = log(sigma)
  evm_fit <- evm(temp, ow, gev, mu = ~ loc, phi = ~ loc, xi = ~loc)
  adj_evm_fit <- alogLik(evm_fit, cluster = ow$year, cadjust = FALSE)
  summary(adj_evm_fit)
}
} # }
```
