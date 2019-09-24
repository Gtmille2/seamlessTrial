#' A function to get the integral from an integrand, upper limit, lower limit, and integrand arguments
#'
#' This function isn't generally used, but is used in getting density and distribution
#' @param integrand integrand given from max distribution and max density functions
#' @param lower lower limit in the integral
#' @param upper upper limit in the integral
#' @param integrand.arguments arguments used in the integral function


get.integral <- function(integrand, lower, upper, integrand.arguments)
{
  # simple trapezium rule integration routine

  n.grid = 1001
  grid = seq(lower,upper,(upper-lower)/(n.grid-1))
  values = integrand(grid, integrand.arguments)
  integral = (values[1] + 2*sum(values[2:(n.grid-1)]) + values[n.grid])*(upper-lower)/(2*(n.grid-1))
  integral
}
