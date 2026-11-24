# lax: Loglikelihood Adjustment for Extreme Value Models

Performs adjusted inferences based on model objects fitted, using
maximum likelihood estimation, by the extreme value analysis packages
[eva](https://cran.r-project.org/package=evd),
[evd](https://cran.r-project.org/package=evd),
[evir](https://cran.r-project.org/package=evir),
[extRemes](https://cran.r-project.org/package=extRemes),
[fExtremes](https://cran.r-project.org/package=fExtremes),
[ismev](https://cran.r-project.org/package=ismev),
[mev](https://cran.r-project.org/package=mev),
[POT](https://cran.r-project.org/package=POT) and
[texmex](https://cran.r-project.org/package=texmex). Univariate extreme
value models, including regression models, are supported. Adjusted
standard errors and an adjusted loglikelihood are provided, using the
[chandwich package](https://cran.r-project.org/package=chandwich) and
the object-oriented features of the [sandwich
package](https://cran.r-project.org/package=sandwich).

## Details

The adjustment is based on a robust sandwich estimator of the parameter
covariance matrix, based on the methodology in Chandler and Bate (2007).
This can be used for cluster correlated data when interest lies in the
parameters of the marginal distributions, or for performing inferences
that are robust to certain types of model misspecification.

The main function is
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md),
which works in an object-oriented way, operating on fitted model
objects. This function performs the loglikelihood adjustments using
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).
See the following package-specific help pages for details and examples:
[`eva`](https://github.com/paulnorthrop/lax/reference/eva.md),
[`evd`](https://github.com/paulnorthrop/lax/reference/evd.md),
[`evir`](https://github.com/paulnorthrop/lax/reference/evir.md),
[`extRemes`](https://github.com/paulnorthrop/lax/reference/extRemes.md),
[`fExtremes`](https://github.com/paulnorthrop/lax/reference/fExtremes.md),
[`ismev`](https://github.com/paulnorthrop/lax/reference/ismev.md),
[`mev`](https://github.com/paulnorthrop/lax/reference/mev.md),
[`POT`](https://github.com/paulnorthrop/lax/reference/POT.md),
[`texmex`](https://github.com/paulnorthrop/lax/reference/texmex.md).

See
[`vignette("lax-vignette", package = "lax")`](https://github.com/paulnorthrop/lax/articles/lax-vignette.md)
for an overview of the package.

## References

Bader, B. and Yan, J. (2020). eva: Extreme Value Analysis with
Goodness-of-Fit Testing. R package version 0.2.6.
<https://CRAN.R-project.org/package=eva>

Belzile, L., Wadsworth, J. L., Northrop, P. J., Grimshaw, S. D. and
Huser, R. (2019). mev: Multivariate Extreme Value Distributions. R
package version 1.12.2. <https://github.com/lbelzile/mev/>

Berger S., Graham N., Zeileis A. (2017). Various Versatile Variances: An
Object-Oriented Implementation of Clustered Covariances in R. Technical
Report 2017-12, Working Papers in Economics and Statistics, Research
Platform Empirical and Experimental Economics, Universitat Innsbruck.
<https://EconPapers.RePEc.org/RePEc:inn:wpaper:2017-12>.

Chandler, R. E. and Bate, S. (2007). Inference for clustered data using
the independence loglikelihood. *Biometrika*, **94**(1), 167-183.
[doi:10.1093/biomet/asm015](https://doi.org/10.1093/biomet/asm015)

Gilleland, E. and Katz, R. W. (2016). extRemes 2.0: An Extreme Value
Analysis Package in R. *Journal of Statistical Software*, **72**(8),
1-39. [doi:10.18637/jss.v072.i08](https://doi.org/10.18637/jss.v072.i08)

Northrop, P. J. and Chandler, R. E. (2018). chandwich: Chandler-Bate
Sandwich Loglikelihood Adjustment. R package version 1.1.
<https://CRAN.R-project.org/package=chandwich>.

Pfaff, B. and McNeil, A. (2018). evir: Extreme Values in R. R package
version 1.7-4. <https://CRAN.R-project.org/package=evir>

Ribatet, M. and Dutang, C. (2019). POT: Generalized Pareto Distribution
and Peaks Over Threshold. R package version 1.1-7.
<https://CRAN.R-project.org/package=POT>

Southworth, H., Heffernan, J. E. and Metcalfe, P. D. (2017). texmex:
Statistical modelling of extreme values. R package version 2.4.
<https://CRAN.R-project.org/package=texmex>.

Stephenson, A. G. evd: Extreme Value Distributions. *R News*,
**2**(2):31-32, June 2002. <https://CRAN.R-project.org/doc/Rnews/>

Stephenson, A. G., Heffernan, J. E. and Gilleland, E. (2018). ismev: An
Introduction to Statistical Modeling of Extreme Values. R package
version 1.42. <https://CRAN.R-project.org/package=ismev>.

Wuertz, D., Setz, T. and Chalabi, Y. (2017). fExtremes: Rmetrics -
Modelling Extreme Events in Finance. R package version 3042.82.
<https://CRAN.R-project.org/package=fExtremes>

Zeileis A. (2004). Econometric Computing with HC and HAC Covariance
Matrix Estimators. *Journal of Statistical Software*, **11**(10), 1-17.
[doi:10.18637/jss.v011.i10](https://doi.org/10.18637/jss.v011.i10) .

Zeileis A. (2006). Object-Oriented Computation of Sandwich Estimators.
*Journal of Statistical Software*, **16**(9), 1-16.
[doi:10.18637/jss.v016.i09](https://doi.org/10.18637/jss.v016.i09) .

## See also

Useful links:

- <https://paulnorthrop.github.io/lax/>

- <https://github.com/paulnorthrop/lax>

- Report bugs at <https://github.com/paulnorthrop/lax/issues>

## Author

**Maintainer**: Paul J. Northrop <p.northrop@ucl.ac.uk> \[copyright
holder\]

Authors:

- Paul J. Northrop <p.northrop@ucl.ac.uk> \[copyright holder\]

- Camellia Yin \[copyright holder\]
