#' Function to calculate the max desnity integrand
#'
#' This function is not generally used. It returns the integrand used in the max.density calculation
#' @param y The value of the maxmium test statistic in this density
#' @param integrand.arguments The integrand arguments passed to the function including k, incv, and x




max.density.integrand <- function(y, integrand.arguments)
{

  k = integrand.arguments[1]
  incv = integrand.arguments[2]
  x = integrand.arguments[3]

  root=sqrt(incv/2)
  integrand = k*dnorm((x-y)/root)/root * dnorm(y/root)/root * pnorm(y/root)^(k-1)
  integrand
}
