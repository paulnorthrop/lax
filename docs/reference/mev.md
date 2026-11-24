# Loglikelihood adjustment for mev fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions
[`fit.gev`](https://lbelzile.github.io/mev/reference/fit.gev.html),
[`fit.gpd`](https://lbelzile.github.io/mev/reference/fit.gpd.html), and
[`fit.pp`](https://lbelzile.github.io/mev/reference/fit.pp.html) and
[`fit.rlarg`](https://lbelzile.github.io/mev/reference/fit.rlarg.html)
in the mev package.

## Usage

``` r
# S3 method for class 'mev_gev'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'mev_pp'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'mev_gpd'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'mev_egp'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'mev_rlarg'
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
`c("lax", "chandwich", "mev")`. The 4th component depends on which model
was fitted. `"gev"` if
[`fit.gev`](https://lbelzile.github.io/mev/reference/fit.gev.html) was
used; `"gpd"` if
[`fit.gpd`](https://lbelzile.github.io/mev/reference/fit.gpd.html) was
used; `"pp"`
[`fit.pp`](https://lbelzile.github.io/mev/reference/fit.pp.html) was
used; `"egp"`
[`fit.egp`](https://lbelzile.github.io/mev/reference/fit.egp.html) was
used; `"rlarg"`
[`fit.rlarg`](https://lbelzile.github.io/mev/reference/fit.rlarg.html)
was used; The 5th component is `"stat"` (for stationary).

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
for details.

If `x` was returned from
[`fit.pp`](https://lbelzile.github.io/mev/reference/fit.pp.html) then
the data `xdat` supplied to
[`fit.pp`](https://lbelzile.github.io/mev/reference/fit.pp.html) must
contain *all* the data, both threshold exceedances and non-exceedances.

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
# We need the mev package
got_mev <- requireNamespace("mev", quietly = TRUE)

if (got_mev) {
  library(mev)
  # An example from the mev::gev.fit documentation
  gev_mev <- fit.gev(revdbayes::portpirie)
  adj_gev_mev <- alogLik(gev_mev)
  summary(adj_gev_mev)

  # Use simulated data
  set.seed(1112019)
  x <- revdbayes::rgp(365 * 10, loc = 0, scale = 1, shape = 0.1)
  pfit <- fit.pp(x, threshold = 1, npp = 365)
  adj_pfit <- alogLik(pfit)
  summary(adj_pfit)

  # An example from the mev::fit.gpd documentation
  gpd_mev <- fit.gpd(eskrain, threshold = 35, method = 'Grimshaw')
  adj_gpd_mev <- alogLik(gpd_mev)
  summary(adj_gpd_mev)

  # An example from the mev::fit.egp documentation
  # (model = "egp1" and model = "egp3" also work)
  xdat <- evd::rgpd(n = 100, loc = 0, scale = 1, shape = 0.5)
  fitted <- fit.egp(xdat = xdat, thresh = 1, model = "egp2", show = FALSE)
  adj_fitted <- alogLik(fitted)
  summary(adj_fitted)

  # An example from the mev::fit.rlarg documentation
  set.seed(31102019)
  xdat <- rrlarg(n = 10, loc = 0, scale = 1, shape = 0.1, r = 4)
  fitr <- fit.rlarg(xdat)
  adj_fitr <- alogLik(fitr)
  summary(adj_fitr)
}
#> 
#> Attaching package: 'mev'
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     venice
#> Error in loadNamespace(x): there is no package called 'evd'
```
