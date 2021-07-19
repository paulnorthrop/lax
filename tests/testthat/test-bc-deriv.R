# Check that box_cox_deriv is correct for lambda = 0

x <- 1:10
for (i in 1:length(x)) {
  test_that(paste("check box_cox_deriv for lambda = 0, x = ", x = x[i]), {
    testthat::expect_equal(box_cox_deriv(x = x[i], lambda = 0),
                           log(x[i]) ^ 2 / 2)
  })
}

# Check that box_cox_deriv is correct for lambda very slightly smaller in
# magnitude than lambda_tol = 1 / 50 and poly_order is large

eps <- 1e-10
poly_order <- 10
lambda_tol <- 1 / 50
x <- 1:10

# lambda very slightly less than lambda_tol

lambda <- lambda_tol - eps
check_val <- (lambda * x ^ lambda * log(x) - x ^ lambda + 1) / lambda ^ 2
for (i in 1:length(x)) {
  test_that(paste("box_cox_deriv, 0 < lambda < lambda_tol, x = ", x = x[i]), {
    testthat::expect_equal(box_cox_deriv(x = x[i], lambda = lambda,
                                         lambda_tol = lambda_tol,
                                         poly_order = poly_order),
                           check_val[i])
  })
}

# lambda very slightly greater than -lambda_tol

lambda <- -lambda_tol + eps
check_val <- (lambda * x ^ lambda * log(x) - x ^ lambda + 1) / lambda ^ 2
for (i in 1:length(x)) {
  test_that(paste("box_cox_deriv, -lambda_tol < lambda < 0, x = ", x = x[i]), {
    testthat::expect_equal(box_cox_deriv(x = x[i], lambda = lambda,
                                         lambda_tol = lambda_tol,
                                         poly_order = poly_order),
                           check_val[i])
  })
}
