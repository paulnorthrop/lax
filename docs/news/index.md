# Changelog

## lax 1.2.5

### Documentation correction

- Corrected the number of rows in the dataset `ow` from 80 to 160.

## lax 1.2.4

CRAN release: 2025-10-21

### Bug fixes

- Updates to methods for
  [`mev::fit.egp()`](https://lbelzile.github.io/mev/reference/fit.egp.html)
  to comply with versions 2.0 and add new models.
- Add Details, Value and Examples sections in `logLik.logLikVec.Rd` and
  `logLikVec.Rd` to avoid CRAN HTML validation NOTE on r-devel Debian.

## lax 1.2.3

CRAN release: 2024-02-25

### Bug fixes

- Calls to `texmex::evm()`, which resulted in CRAN package check ERRORs
  on some platforms, have been avoided.

## lax 1.2.2

CRAN release: 2023-12-02

### Bug fixes

- Fixed issues with the incorrect use of in some Rd files.

## lax 1.2.1

CRAN release: 2023-09-02

### Bug fixes and minor improvements

- The original model object `x` is added as an attribute
  `"original_fit"` to the object returned from `alogLik(x)`.

- In the documentation of
  [`return_level()`](https://github.com/paulnorthrop/lax/reference/return_level.md)
  the role of `npy` has been explained and a more accurate calculation
  is used for the estimation of return levels in the case where `npy` is
  not equal to 1.

- If the argument `cluster` was supplied an
  [`alogLik()`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
  method then this is now returned as the attribute `cluster` in the
  returned object, rather than the default returned by
  [`chandwich::adjust_loglik()`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).

- Create the help file for the package correctly, with alias
  lax-package.

- README.md: Used app.codecov.io as base for codecov link.

- Activated 3rd edition of the `testthat` package

## lax 1.2.0

CRAN release: 2021-07-20

### New features

- The eva package is now supported: functions `gpdFit` and `gevrFit`.

### Bug fixes and minor improvements

- The links at the end of the Details section of the main lax package
  help page have been corrected.

- Depreciated function
  [`testthat::context`](https://testthat.r-lib.org/reference/context.html)
  is no longer used.

- Some obsolete code has been deleted from the lax help file for mev.

## lax 1.1.0

CRAN release: 2019-12-05

### New features

- The mev package is now supported: functions `fit.gev`, `fit.gpd`,
  `fit.egp`, `fit.pp` and `fit.rlarg`.

- The function `rlarg.fit` in the ismev package is now supported.

### Bug fixes and minor improvements

- Unnecessary generic information concerning the availability of S3
  methods has been removed from the Details sections of the
  package-specific loglikelihood adjustment documentation.

- More tests of internal function `box_cox_deriv`.
