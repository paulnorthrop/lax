# Package index

## Package

- [`lax`](https://github.com/paulnorthrop/lax/reference/lax-package.md)
  [`lax-package`](https://github.com/paulnorthrop/lax/reference/lax-package.md)
  : lax: Loglikelihood Adjustment for Extreme Value Models

## Loglikelihood adjustment

- [`alogLik()`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
  : Loglikelihood adjustment for model fits.
- [`alogLik(`*`<gevrFit>`*`)`](https://github.com/paulnorthrop/lax/reference/eva.md)
  [`alogLik(`*`<gpdFit>`*`)`](https://github.com/paulnorthrop/lax/reference/eva.md)
  : Loglikelihood adjustment for eva fits
- [`alogLik(`*`<evd>`*`)`](https://github.com/paulnorthrop/lax/reference/evd.md)
  : Loglikelihood adjustment for evd fits
- [`alogLik(`*`<gev>`*`)`](https://github.com/paulnorthrop/lax/reference/evir.md)
  [`alogLik(`*`<gpd>`*`)`](https://github.com/paulnorthrop/lax/reference/evir.md)
  [`alogLik(`*`<potd>`*`)`](https://github.com/paulnorthrop/lax/reference/evir.md)
  : Loglikelihood adjustment for evir fits
- [`alogLik(`*`<fevd>`*`)`](https://github.com/paulnorthrop/lax/reference/extRemes.md)
  : Loglikelihood adjustment for extRemes fits
- [`alogLik(`*`<fGEVFIT>`*`)`](https://github.com/paulnorthrop/lax/reference/fExtremes.md)
  [`alogLik(`*`<fGPDFIT>`*`)`](https://github.com/paulnorthrop/lax/reference/fExtremes.md)
  : Loglikelihood adjustment for fExtremes fits
- [`alogLik(`*`<gev.fit>`*`)`](https://github.com/paulnorthrop/lax/reference/ismev.md)
  [`alogLik(`*`<pp.fit>`*`)`](https://github.com/paulnorthrop/lax/reference/ismev.md)
  [`alogLik(`*`<gpd.fit>`*`)`](https://github.com/paulnorthrop/lax/reference/ismev.md)
  [`alogLik(`*`<rlarg.fit>`*`)`](https://github.com/paulnorthrop/lax/reference/ismev.md)
  : Loglikelihood adjustment for ismev fits
- [`alogLik(`*`<mev_gev>`*`)`](https://github.com/paulnorthrop/lax/reference/mev.md)
  [`alogLik(`*`<mev_pp>`*`)`](https://github.com/paulnorthrop/lax/reference/mev.md)
  [`alogLik(`*`<mev_gpd>`*`)`](https://github.com/paulnorthrop/lax/reference/mev.md)
  [`alogLik(`*`<mev_egp>`*`)`](https://github.com/paulnorthrop/lax/reference/mev.md)
  [`alogLik(`*`<mev_rlarg>`*`)`](https://github.com/paulnorthrop/lax/reference/mev.md)
  : Loglikelihood adjustment for mev fits
- [`alogLik(`*`<uvpot>`*`)`](https://github.com/paulnorthrop/lax/reference/POT.md)
  : Loglikelihood adjustment for POT fits
- [`alogLik(`*`<evmOpt>`*`)`](https://github.com/paulnorthrop/lax/reference/texmex.md)
  : Loglikelihood adjustment of texmex fits

## Refitting models

- [`gev_refit()`](https://github.com/paulnorthrop/lax/reference/ismev_refits.md)
  [`gpd_refit()`](https://github.com/paulnorthrop/lax/reference/ismev_refits.md)
  [`pp_refit()`](https://github.com/paulnorthrop/lax/reference/ismev_refits.md)
  [`rlarg_refit()`](https://github.com/paulnorthrop/lax/reference/ismev_refits.md)
  : Maximum-likelihood (Re-)Fitting using the ismev package
- [`pot_refit()`](https://github.com/paulnorthrop/lax/reference/pot_refit.md)
  : Fits a Poisson point process to the data, an approach sometimes
  known as peaks over thresholds (POT), and returns an object of class
  "potd".

## Statistical inference

- [`anova(`*`<lax>`*`)`](https://github.com/paulnorthrop/lax/reference/anova.lax.md)
  : Comparison of nested models

- [`fit_bernoulli()`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`logLikVec(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`nobs(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`coef(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`vcov(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`logLik(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  [`alogLik(`*`<bernoulli>`*`)`](https://github.com/paulnorthrop/lax/reference/bernoulli.md)
  : Inference for the Bernoulli distribution

- [`return_level()`](https://github.com/paulnorthrop/lax/reference/return_level.md)
  : Return Level Inferences for Stationary Extreme Value Models

- [`plot(`*`<retlev>`*`)`](https://github.com/paulnorthrop/lax/reference/plot.retlev.md)
  : Plot diagnostics for a retlev object

- [`print(`*`<retlev>`*`)`](https://github.com/paulnorthrop/lax/reference/print.retlev.md)
  : Print method for retlev object

- [`summary(`*`<retlev>`*`)`](https://github.com/paulnorthrop/lax/reference/summary.retlev.md)
  :

  Summary method for a `"retlev"` object

- [`print(`*`<summary.retlev>`*`)`](https://github.com/paulnorthrop/lax/reference/print.summary.retlev.md)
  :

  Print method for objects of class `"summary.retlev"`

## Example data

- [`ow`](https://github.com/paulnorthrop/lax/reference/ow.md) : Oxford
  and Worthing annual maximum temperatures

## Generics

- [`logLikVec()`](https://github.com/paulnorthrop/lax/reference/logLikVec.md)
  : Evaluate loglikelihood contributions from specific observations
- [`logLik(`*`<logLikVec>`*`)`](https://github.com/paulnorthrop/lax/reference/logLik.logLikVec.md)
  : Sum loglikelihood contributions from individual observations
