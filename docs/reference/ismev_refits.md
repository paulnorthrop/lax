# Maximum-likelihood (Re-)Fitting using the ismev package

These are a slightly modified versions of the
[`gev.fit`](https://rdrr.io/pkg/ismev/man/gev.fit.html),
[`gpd.fit`](https://rdrr.io/pkg/ismev/man/gpd.fit.html),
[`pp.fit`](https://rdrr.io/pkg/ismev/man/pp.fit.html) and
[`rlarg.fit`](https://rdrr.io/pkg/ismev/man/rlarg.fit.html) functions in
the [`ismev`](https://rdrr.io/pkg/ismev/man/ismev.html) package. The
modification is to add to the returned object regression design matrices
for the parameters of the model. That is,
`xdat, ydat, mulink, siglink, shlink` and matrices
`mumat, sigmat, shmat` for the location, scale and shape parameters
[`gev.fit`](https://rdrr.io/pkg/ismev/man/gev.fit.html),
[`pp.fit`](https://rdrr.io/pkg/ismev/man/pp.fit.html) and
[`rlarg.fit`](https://rdrr.io/pkg/ismev/man/rlarg.fit.html), and `xdat`,
`ydat, siglink, shlink` and matrices `sigmat, shmat` for the scale and
shape parameters for
[`gpd.fit`](https://rdrr.io/pkg/ismev/man/gpd.fit.html).

## Usage

``` r
gev_refit(
  xdat,
  ydat = NULL,
  mul = NULL,
  sigl = NULL,
  shl = NULL,
  mulink = identity,
  siglink = identity,
  shlink = identity,
  muinit = NULL,
  siginit = NULL,
  shinit = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  maxit = 10000,
  ...
)

gpd_refit(
  xdat,
  threshold,
  npy = 365,
  ydat = NULL,
  sigl = NULL,
  shl = NULL,
  siglink = identity,
  shlink = identity,
  siginit = NULL,
  shinit = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  maxit = 10000,
  ...
)

pp_refit(
  xdat,
  threshold,
  npy = 365,
  ydat = NULL,
  mul = NULL,
  sigl = NULL,
  shl = NULL,
  mulink = identity,
  siglink = identity,
  shlink = identity,
  muinit = NULL,
  siginit = NULL,
  shinit = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  maxit = 10000,
  ...
)

rlarg_refit(
  xdat,
  r = dim(xdat)[2],
  ydat = NULL,
  mul = NULL,
  sigl = NULL,
  shl = NULL,
  mulink = identity,
  siglink = identity,
  shlink = identity,
  muinit = NULL,
  siginit = NULL,
  shinit = NULL,
  show = TRUE,
  method = "Nelder-Mead",
  maxit = 10000,
  ...
)
```

## Arguments

- xdat:

  A numeric vector of data to be fitted.

- ydat:

  A matrix of covariates for generalized linear modelling of the
  parameters (or `NULL` (the default) for stationary fitting). The
  number of rows should be the same as the length of `xdat`.

- mul, sigl, shl:

  Numeric vectors of integers, giving the columns of `ydat` that contain
  covariates for generalized linear modelling of the location, scale and
  shape parameters repectively (or `NULL` (the default) if the
  corresponding parameter is stationary).

- mulink, siglink, shlink:

  Inverse link functions for generalized linear modelling of the
  location, scale and shape parameters repectively.

- muinit, siginit, shinit:

  numeric of length equal to total number of parameters used to model
  the location, scale or shape parameter(s), resp. See Details section
  for default (NULL) initial values.

- show:

  Logical; if `TRUE` (the default), print details of the fit.

- method:

  The optimization method (see
  [`optim`](https://rdrr.io/r/stats/optim.html) for details).

- maxit:

  The maximum number of iterations.

- ...:

  Other control parameters for the optimization. These are passed to
  components of the `control` argument of `optim`.

- threshold:

  The threshold; a single number or a numeric vector of the same length
  as `xdat`.

- npy:

  The number of observations per year/block.

- r:

  The largest `r` order statistics are used for the fitted model.

## References

Heffernan, J. E. and Stephenson, A. G. (2018). ismev: An Introduction to
Statistical Modeling of Extreme Values. R package version 1.42.
<https://CRAN.R-project.org/package=ismev>.

## Examples

``` r
# We need the ismev package
got_ismev <- requireNamespace("ismev", quietly = TRUE)
if (got_ismev) {
  library(ismev)
  fit1 <- gev.fit(revdbayes::portpirie, show = FALSE)
  ls(fit1)
  fit2 <- gev_refit(revdbayes::portpirie, show = FALSE)
  ls(fit2)

  data(rain)
  fit1 <- gpd.fit(rain, 10)
  ls(fit1)
  fit2 <- gpd_refit(rain, 10)
  ls(fit2)

  fit1 <- pp.fit(rain, 10, show = FALSE)
  ls(fit1)
  fit2 <- pp_refit(rain, 10, show = FALSE)
  ls(fit2)

  data(venice)
  fit1 <- rlarg.fit(venice[, -1], muinit = 120.54, siginit = 12.78,
                   shinit = -0.1129, show = FALSE)
  ls(fit1)
  fit2 <- rlarg_refit(venice[, -1], muinit = 120.54, siginit = 12.78,
                   shinit = -0.1129, show = FALSE)
  ls(fit2)
}
#> $threshold
#> [1] 10
#> 
#> $nexc
#> [1] 2003
#> 
#> $conv
#> [1] 0
#> 
#> $nllh
#> [1] 6123.465
#> 
#> $mle
#> [1] 7.43768624 0.05045225
#> 
#> $rate
#> [1] 0.1142547
#> 
#> $se
#> [1] 0.23606472 0.02256649
#> 
#> $threshold
#> [1] 10
#> 
#> $nexc
#> [1] 2003
#> 
#> $conv
#> [1] 0
#> 
#> $nllh
#> [1] 6123.465
#> 
#> $mle
#> [1] 7.43768624 0.05045225
#> 
#> $rate
#> [1] 0.1142547
#> 
#> $se
#> [1] 0.23606472 0.02256649
#> 
#>  [1] "conv"    "cov"     "data"    "link"    "mle"     "model"   "mulink" 
#>  [8] "mumat"   "nllh"    "r"       "se"      "shlink"  "shmat"   "siglink"
#> [15] "sigmat"  "trans"   "vals"    "xdat"   
```
