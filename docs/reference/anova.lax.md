# Comparison of nested models

`anova` method for objects of class `"lax"`. Compares two or more nested
models using the adjusted likelihood ratio test statistic (ALRTS)
described in Section 3.5 of Chandler and Bate (2007). The nesting must
result from the simple constraint that a subset of the parameters of the
larger model is held fixed.

## Usage

``` r
# S3 method for class 'lax'
anova(object, object2, ...)
```

## Arguments

- object:

  An object of class `"lax"`, inheriting from class `"chandwich"`,
  returned by
  [`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md).

- object2:

  An object of class `"lax"`, inheriting from class `"chandwich"`,
  returned by
  [`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md).

- ...:

  Further objects of class `"lax"` and/or arguments to be passed to
  [`anova.chandwich`](https://github.com/paulnorthrop/chandwich/reference/anova.chandwich.html),
  and then on to
  [`compare_models`](https://github.com/paulnorthrop/chandwich/reference/compare_models.html),
  in particular `type`, which chooses the type of adjustment.

## Value

An object of class `"anova"` inheriting from class `"data.frame"`, with
four columns:

- Model.Df:

  The number of parameters in the model

- Df:

  The decrease in the number of parameter compared the model in the
  previous row

- ALRTS:

  The adjusted likelihood ratio test statistic

- Pr(\>ALRTS):

  The p-value associated with the test that the model is a valid
  simplification of the model in the previous row.

The row names are the names of the model objects.

## Details

The objects of class `"lax"` need not be provided in nested order: they
will be ordered inside `anova.lax` based on the values of
`attr(., "p_current")`.

## References

Chandler, R. E. and Bate, S. (2007). Inference for clustered data using
the independence loglikelihood. *Biometrika*, **94**(1), 167-183.
[doi:10.1093/biomet/asm015](https://doi.org/10.1093/biomet/asm015)

## See also

[`anova.chandwich`](https://github.com/paulnorthrop/chandwich/reference/anova.chandwich.html):
the `anova` method on which `anova.lax` is based.

[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md):
loglikelihood adjustment for model fits.

## Examples

``` r
got_evd <- requireNamespace("evd", quietly = TRUE)
if (got_evd) {
  library(evd)
  small <- fgev(ow$temp, nsloc = ow[, "loc"])
  adj_small <- alogLik(small, cluster = ow$year)
  tiny <- fgev(ow$temp)
  adj_tiny <- alogLik(tiny, cluster = ow$year)
  anova(adj_small, adj_tiny)

  set.seed(4082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M0 <- fgev(uvdata)
  M1 <- fgev(uvdata, nsloc = (-49:50)/100)
  adj0 <- alogLik(M0)
  adj1 <- alogLik(M1)
  anova(adj1, adj0)
}

got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
if (got_extRemes) {
  library(extRemes)
  large <- fevd(temp, ow, location.fun = ~ loc, scale.fun = ~ loc,
                shape.fun = ~ loc)
  medium <- fevd(temp, ow, location.fun = ~ loc, scale.fun = ~ loc)
  small <- fevd(temp, ow, location.fun = ~ loc)
  tiny <- fevd(temp, ow)
  adj_large <- alogLik(large, cluster = ow$year)
  adj_medium <- alogLik(medium, cluster = ow$year)
  adj_small <- alogLik(small, cluster = ow$year)
  adj_tiny <- alogLik(tiny, cluster = ow$year)
  anova(adj_large, adj_medium, adj_small, adj_tiny)
}
```
