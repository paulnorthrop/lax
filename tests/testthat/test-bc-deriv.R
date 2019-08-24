context("Box-Cox deriv")

# Check that box_cox_deriv is correct for lamnda = 0

x <- 1:10
for (i in length(x)) {
  test_that(paste("check box_cox_deriv for lambda = 0, x = ", x = x[i]), {
    testthat::expect_equal(box_cox_deriv(x = x[i], lambda = 0),
                           log(x[i]) ^ 2 / 2)
  })
}

