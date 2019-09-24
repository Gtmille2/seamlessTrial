#' A distribution search function
#'
#' This function isn't generally used, it is used to get the distribution of the test statistic at each look
#' @param x integrating over x
#' @param k the number of treaments in this trial
#' @param incv Incv is the number of increments in the integration
#' @param look the current interim analysis
#' @param upper.prev previous upper bound
#' @param lower.prev previous lower bound
#' @param n.grid number of points to integrate over
#' @param density.prev previous density of test statistic
#' @param target target value in the spending function


distribution.search.function <- function(x, k, incv, look, upper.prev, lower.prev, n.grid, density.prev, target)
{
  # function to return distribution function for test statistic at each look - this is given by integral for look > 1
  if (look == 1) get.max.distribution(k, incv, x) - target
  else
  {
    z.grid.prev = seq(lower.prev,upper.prev,(upper.prev-lower.prev)/(n.grid-1))
    values = rep(0,n.grid)
    for (j in seq(1,n.grid)) values[j] = get.max.distribution(k,incv,x-z.grid.prev[j])*density.prev[j]
    integral = (values[1] + 2*sum(values[2:(n.grid-1)]) + values[n.grid])*(upper.prev-lower.prev)/(2*(n.grid-1))
    integral - target
  }
}

