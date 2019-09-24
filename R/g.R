#' A function to calculate G in the Pocock Design and return a vector of p
#'
#' This function is not generally used. It is used in the calculation of the treatment assignments in psd
#' @param x x is the vector of deltas in Pocock's design
#' @param p1 = 3/4

g=function(x,p1=3/4,n.trt)
{
  p = NULL
  best = which(x == min(x))
  if (length(best) > 1) {
    p = rep(1/n.trt,n.trt)
  } else {
    p = rep(1-p1/(n.trt-1),n.trt)
    p[best] = 3/4

  }
  p

}
