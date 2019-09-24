#' Function to get max density
#'
#' This function performs integration to return density function for maximum of k correlated normals
#' @param k The number of different treatment groups
#' @param incv The number of steps in the integration
#' @param x The mean in the dnorm density

get.max.density <- function(k, incv, x)
{
  # performs integration to return density function for maximum of k correlated normals
  # returns normal density without integration if k = 1

  if (k==1) dens = dnorm(x/sqrt(incv))/sqrt(incv)
  else dens = get.integral(max.density.integrand, -20*sqrt(incv), 20*sqrt(incv),integrand.arguments = c(k, incv, x) )
  dens
}
