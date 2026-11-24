# Loglikelihood adjustment for model fits.

This function is generic. It performs adjustment of the loglikelihood
associated with fitted model objects, following Chandler and Bate
(2007). Certain classes of extreme value model objects are supported
automatically. For details see the `alogLik` help pages for the
packages: [`evd`](https://github.com/paulnorthrop/lax/reference/evd.md),
[`evir`](https://github.com/paulnorthrop/lax/reference/evir.md),
[`extRemes`](https://github.com/paulnorthrop/lax/reference/extRemes.md),
[`fExtremes`](https://github.com/paulnorthrop/lax/reference/fExtremes.md),
[`ismev`](https://github.com/paulnorthrop/lax/reference/ismev.md),
[`mev`](https://github.com/paulnorthrop/lax/reference/mev.md),
[`POT`](https://github.com/paulnorthrop/lax/reference/POT.md),
[`texmex`](https://github.com/paulnorthrop/lax/reference/texmex.md).
User-supplied objects can also be supported: the requirements for these
objects are explained in **Details**.

## Usage

``` r
alogLik(
  x,
  cluster = NULL,
  use_vcov = TRUE,
  binom = FALSE,
  k,
  inc_cens = TRUE,
  ...
)
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

- binom:

  A logical scalar. This option is only relevant to **GP models** and is
  only available in the **stationary** (no covariates) case. If
  `binom = FALSE` then loglikelihood adjustment is only performed using
  the GP model. If `binom = TRUE` then loglikelihood adjustment is also
  performed for inferences about the probability of threshold
  exceedance, using a Bernoulli model for the instances of threshold
  exceedance.

- k:

  A non-negative integer scalar. This option is only relevant to **GP
  models** and is only available in the **stationary** (no covariates)
  case. If `k` is supplied then it is passed as the run parameter \\K\\
  to
  [`kgaps`](https://github.com/paulnorthrop/exdex/reference/kgaps.html)
  for making inferences about the extremal index \\\theta\\ using the
  \\K\\-gaps model of Suveges and Davison (2010).

- inc_cens:

  A logical scalar. This argument is only relevant if `k` is supplied.
  Passed to
  [`kgaps`](https://github.com/paulnorthrop/exdex/reference/kgaps.html)
  to indicate whether or not to include censored inter-exceedance times,
  relating to the first and last observations.

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

The original fitted model object is available as an attribute named
`"original_fit"`, accessible using `attr(name, "original_fit")`, where
`name` is the name of the object to which the object returned from
`alogLik` is assigned.

If `binom = TRUE` then the returned object has an extra attribute named
`pu_aloglik` that contains an object inheriting from class `"chandwich"`
relating specifically to inferences about the probability of threshold
exceedance. Also, the 4th component of the class of the returned object
becomes `"bin-gpd"`.

If `k` is supplied then the returned object has an extra attribute named
`theta` that contains an object inheriting from class
`c("kgaps", "exdex")` relating specifically to inferences about the
extremal index \\\theta\\. See the **Value** section in
[`kgaps`](https://github.com/paulnorthrop/exdex/reference/kgaps.html).

If `x` is one of the supported models then the class of the returned
object is a vector of length 5. The first 3 components are
`c("lax", "chandwich", "name_of_package")`, where `"name_of_package"` is
the name of the package from which the input object `x` originated. The
remaining 2 components depend on the model that was fitted. See the
documentation of the relevant package for details:
[`evd`](https://github.com/paulnorthrop/lax/reference/evd.md),
[`evir`](https://github.com/paulnorthrop/lax/reference/evir.md),
[`extRemes`](https://github.com/paulnorthrop/lax/reference/extRemes.md),
[`fExtremes`](https://github.com/paulnorthrop/lax/reference/fExtremes.md),
[`ismev`](https://github.com/paulnorthrop/lax/reference/ismev.md),
[`mev`](https://github.com/paulnorthrop/lax/reference/mev.md),
[`POT`](https://github.com/paulnorthrop/lax/reference/POT.md),
[`texmex`](https://github.com/paulnorthrop/lax/reference/texmex.md).

Otherwise, the class of the returned object is
`c("lax", "chandwich", class(x))`.

Objects returned from \`aloglik\` have \`anova\`, \`coef\`, \`confint\`,
\`logLik\`, \`nobs\`, \`plot\`, \`print\`, \`summary\` and \`vcov\`
methods.

## Details

Object `x` *must* have the following S3 methods:

- `logLikVec`: returns a vector of the contributions to the independence
  loglikelihood from individual observations;

- `coef`: returns a vector of model coefficients, see
  [`coef`](https://rdrr.io/r/stats/coef.html);

- `nobs`: returns the number of (non-missing) observations used in a
  model fit, see [`nobs`](https://rdrr.io/r/stats/nobs.html);

and *may* have the following S3 methods

- `vcov`: returns the estimated variance-covariance matrix of the (main)
  parameters of a fitted model, see
  [`vcov`](https://rdrr.io/r/stats/vcov.html);

- `estfun`: returns an \\n\\ by \\k\\ matrix, in which each column gives
  the derivative of the loglikelihood at each of \\n\\ observation with
  respect to the \\k\\ parameters of the model, see
  [`estfun`](https://zeileis.codeberg.page/sandwich/reference/estfun.html).

Loglikelihood adjustment is performed using the
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
function in the
[`chandwich`](https://github.com/paulnorthrop/chandwich/reference/chandwich-package.html)
package. The relevant arguments to
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html),
namely `loglik, mle, H` and `V`, are created based on the class of the
object `x`.

If a `vcov` method is not available, or if `use_vcov = FALSE`, then the
variance-covariance matrix of the MLE (from which `H` is calculated) is
estimated inside
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
using [`optimHess`](https://rdrr.io/r/stats/optim.html).

The `sandwich` package is used to estimate the variance matrix `V` of
the score vector:
[`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html) is
used if `cluster = NULL`;
[`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
is used if `cluster` is not `NULL`. If `cluster` is `NULL` then any
arguments of
[`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
present in ... will be ignored. Similarly, if `cluster` is not `NULL`
then any arguments of
[`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html)
present in ... will be ignored.
[`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html) and
[`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
require an
[`estfun`](https://zeileis.codeberg.page/sandwich/reference/estfun.html)
method to be available, which, in the current context, provides matrix
of score contributions. If a bespoke `estfun` method is not provided
then this is constructed by estimating the score contributions using
[`jacobian`](https://rdrr.io/pkg/numDeriv/man/jacobian.html).

## Examples

See the (package-specific) examples in
[`evd`](https://github.com/paulnorthrop/lax/reference/evd.md),
[`evir`](https://github.com/paulnorthrop/lax/reference/evir.md),
[`extRemes`](https://github.com/paulnorthrop/lax/reference/extRemes.md),[`fExtremes`](https://github.com/paulnorthrop/lax/reference/fExtremes.md),
[`ismev`](https://github.com/paulnorthrop/lax/reference/ismev.md),
[`mev`](https://github.com/paulnorthrop/lax/reference/mev.md),
[`POT`](https://github.com/paulnorthrop/lax/reference/POT.md) and
[`texmex`](https://github.com/paulnorthrop/lax/reference/texmex.md).

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

[`summary.chandwich`](https://github.com/paulnorthrop/chandwich/reference/summary.chandwich.html),
[`plot.chandwich`](https://github.com/paulnorthrop/chandwich/reference/plot.chandwich.html),
[`confint.chandwich`](https://github.com/paulnorthrop/chandwich/reference/confint.chandwich.html),
[`anova.chandwich`](https://github.com/paulnorthrop/chandwich/reference/anova.chandwich.html),
[`coef.chandwich`](https://github.com/paulnorthrop/chandwich/reference/coef.chandwich.html),
[`vcov.chandwich`](https://github.com/paulnorthrop/chandwich/reference/vcov.chandwich.html)
and
[`logLik.chandwich`](https://github.com/paulnorthrop/chandwich/reference/logLik.chandwich.html)
for S3 methods for objects of class `"chandwich"`.

[`conf_region`](https://github.com/paulnorthrop/chandwich/reference/conf_region.html)
for confidence regions for pairs of parameters.

[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
in the `chandwich` package to adjust a user-supplied loglikelihood.

[`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html) and
[`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
in the sandwich package.
