# Fits a Poisson point process to the data, an approach sometimes known as peaks over thresholds (POT), and returns an object of class "potd".

This is a slightly modified versions of the `pot` function in the `evir`
package. The main modification is to add to the returned object the
argument `data` supplied by the user. This is added to the returned
(list) object with the name `input_data`.

## Usage

``` r
pot_refit(data, threshold = NA, nextremes = NA, run = NA, picture = TRUE, ...)
```

## Arguments

- data:

  numeric vector of data, which may have a times attribute containing
  (in an object of class `"POSIXct"`, or an object that can be converted
  to that class; see
  [`as.POSIXct`](https://rdrr.io/r/base/as.POSIXlt.html)) the
  times/dates of each observation. If no times attribute exists, the
  data are assumed to be equally spaced.

- threshold:

  a threshold value (either this or `nextremes` must be given but not
  both).

- nextremes:

  the number of upper extremes to be used (either this or `threshold`
  must be given but not both).

- run:

  if the data are to be declustered the run length parameter for the
  runs method (see `decluster`) should be entered here.

- picture:

  whether or not a picture should be drawn if declustering is performed.

- ...:

  arguments passed to [`optim`](https://rdrr.io/r/stats/optim.html).

## References

Bernhard Pfaff and Alexander McNeil (2018). evir: Extreme Values in R. R
package version 1.7-4. <https://CRAN.R-project.org/package=evir>.

## Examples

``` r
# We need the evir package
got_evir <- requireNamespace("evir", quietly = TRUE)
if (got_evir) {
  library(evir)
  data(danish)
  out <- pot(danish, 10)
  ls(out)
  out <- pot_refit(danish, 10)
  ls(out)
}
```
