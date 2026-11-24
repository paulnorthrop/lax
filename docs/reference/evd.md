# Loglikelihood adjustment for evd fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions `fgev` and
`fpot` in the evd package. If `x` is returned from `fgev` then the call
must have used `prob = NULL`.

## Usage

``` r
# S3 method for class 'evd'
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
`c("lax", "chandwich", "evd")`. The remaining 2 components depend on the
model that was fitted. If `fgev` was used then these components are
`c("gev", "stat")` if `nsloc` was `NULL` and `c("gev", "nonstat")` if
`nsloc` was not `NULL`. If `fpot` was used then these components are
`c("pot", "gpd")` if `model` was `"gpd"` and `c("pot", "pp")` if `model`
was `"pp"`.

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
# We need the evd package
got_evd <- requireNamespace("evd", quietly = TRUE)

if (got_evd) {
  library(evd)
  # An example from the evd::fgev documentation
  set.seed(3082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  adj_fgev <- alogLik(M1)
  summary(adj_fgev)

  # An example from Chandler and Bate (2007)
  owfit <- fgev(ow$temp, nsloc = ow$loc)
  adj_owfit <- alogLik(owfit, cluster = ow$year)
  summary(adj_owfit)

  # An example from the evd::fpot documentation
  set.seed(3082019)
  uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
  M1 <- fpot(uvdata, 1)
  adj_fpot <- alogLik(M1)
  summary(adj_fpot)
  # Fit using the pp model, rather than the gpd
  M1 <- fpot(uvdata, 1, model = "pp", npp = 365)
  adj_fpot <- alogLik(M1)
  summary(adj_fpot)
}
```
