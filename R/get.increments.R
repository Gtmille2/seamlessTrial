#' Function to generate increments
#'
#' This function is not generally used, but is used in the creation of increments
#' @param x An array is passed from get.boundaries


get.increments <- function(x)
{
  inc.x = rep(0,length(x))
  inc.x[1] = x[1]
  if (length(x) > 1) for (i in seq(2,length(x))) inc.x[i]=x[i]-x[i-1]
  inc.x
}
