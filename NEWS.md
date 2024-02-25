# lax 1.2.3

## Bug fixes

# lax 1.2.2

## Bug fixes

* Fixed issues with the incorrect use of \itemize in some Rd files.

# lax 1.2.1

## Bug fixes and minor improvements

* The original model object `x` is added as an attribute `"original_fit"` to the object returned from `alogLik(x)`.

* In the documentation of `return_level()` the role of `npy` has been explained and a more accurate calculation is used for the estimation of return levels in the case where `npy` is not equal to 1. 

* If the argument `cluster` was supplied an `alogLik()` method then this is now returned as the attribute `cluster` in the returned object, rather than the default returned by `chandwich::adjust_loglik()`.

* Create the help file for the package correctly, with alias lax-package.

* README.md: Used app.codecov.io as base for codecov link.

* Activated 3rd edition of the `testthat` package

# lax 1.2.0

## New features

* The eva package is now supported: functions `gpdFit` and `gevrFit`.

## Bug fixes and minor improvements

* The links at the end of the Details section of the main lax package help page have been corrected.

* Depreciated function `testthat::context` is no longer used.

* Some obsolete code has been deleted from the lax help file for mev.

# lax 1.1.0

## New features

* The mev package is now supported: functions `fit.gev`, `fit.gpd`, `fit.egp`, `fit.pp` and `fit.rlarg`.

* The function `rlarg.fit` in the ismev package is now supported.

## Bug fixes and minor improvements

* Unnecessary generic information concerning the availability of S3 methods has been removed from the Details sections of the package-specific loglikelihood adjustment documentation. 

* More tests of internal function `box_cox_deriv`.
