This is a patch to fix the ERROR at https://cran.r-project.org/web/checks/check_results_lax.html. These problems are inherited from the texmex package, so I have suppressed examples and tests that involve texmex

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Debian Linux, GCC (R-patched and R-devel) on R-hub
- Fedora Linux, GCC (R-devel) on R-hub
- macOS (R-release), ubuntu (R-oldrel, R-release, R-devel), windows (R-release) using the rcmdcheck package
- win-builder (R-devel, R-release and R-oldrelease)

## Downstream dependencies

chantrics passed R CMD check, with the same NOTE as on its CRAN check results page
