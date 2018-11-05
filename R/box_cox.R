# ====================== box_cox_deriv ==========================

box_cox_deriv <- function (x, lambda = 1, lambda_tol = 1e-6,
                           poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda.tol
  #
  # Returns:
  #   A numeric vector.  The derivative with respect to lambda of
  #     (x^lambda - 1) / lambda
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * log(x) - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(log(x) ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}
