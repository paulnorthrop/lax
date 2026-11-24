# Internal lax functions

Internal lax functions.

## Usage

``` r
adj_object(x, cluster = NULL, use_vcov = TRUE, ...)

return_level_gev(x, m, level, npy, prof, inc, type)

gev_rl_CI(x, m, level, npy, type)

gev_rl_prof(x, m, level, npy, inc, type, rl_sym)

return_level_bingp(x, m, level, npy, prof, inc, type, npy_given)

bingp_rl_CI(x, m, level, npy, type, u)

bingp_rl_prof(x, m, level, npy, inc, type, rl_sym, u)

box_cox_deriv(x, lambda = 1, lambda_tol = 1/50, poly_order = 3)

ismev_ppp(a, npy)

kgaps_loglik(theta, N0, N1, sum_qs, n_kgaps)
```

## Details

These functions are not intended to be called by the user.
