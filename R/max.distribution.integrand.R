#' Function for the max distribution integrand
#'
#' This function is not generally used. This function is used by seamlessTrials to calculate the max distribution integrand
#' @param y value of maximum test statistic
#' @param integrand.arguments The integrand arguments passed to the function that contain K, incv, and x


max.distribution.integrand <- function(y, integrand.arguments)
{

  k = integrand.arguments[1]
  incv = integrand.arguments[2]
  x = integrand.arguments[3]

  root=sqrt(incv/2)
  integrand = k*pnorm((x-y)/root) * dnorm(y/root)/root * (pnorm(y/root))^(k-1)
  integrand
}
