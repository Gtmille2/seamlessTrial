#' Get Max Distribution Function
#'
#' This function isn't generally used. Used to perform integration to return distribution fucntion for K correlated normal distributions
#' Returns the normal distribution with no integration if k = 1
#' @param k The number of treaments being tested in tehe trial
#' @param incv The number of integration points
#' @param x The points in the integration


get.max.distribution <- function(k, incv, x)
{
  # performs integration to return distribution function for maximum of k correlated normals
  # returns normal distribution without integration if k = 1

  if (k==1) distribution = pnorm(x/sqrt(incv))
  else distribution = get.integral(max.distribution.integrand, -20*sqrt(incv), 20*sqrt(incv),integrand.arguments = c(k, incv, x) )
  distribution
}
