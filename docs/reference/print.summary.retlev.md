# Print method for objects of class `"summary.retlev"`

`print` method for an object `x` of class `"summary.retlev"`.

## Usage

``` r
# S3 method for class 'summary.retlev'
print(x, ...)
```

## Arguments

- x:

  An object of class "summary.retlev", a result of a call to
  [`summary.retlev`](https://github.com/paulnorthrop/lax/reference/summary.retlev.md).

- ...:

  Additional arguments passed on to
  [`print.default`](https://rdrr.io/r/base/print.default.html).

## Value

The argument `x`, invisibly, as for all
[`print`](https://rdrr.io/r/base/print.html) methods.

## Details

Prints the call and the numeric matrix `x$matrix` returned from
[`summary.retlev`](https://github.com/paulnorthrop/lax/reference/summary.retlev.md).

## Examples

See the examples in
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

## See also

[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md)
to perform inferences about return levels.
