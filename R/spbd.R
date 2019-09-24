#' spbd Function
#'
#' Function to determine treatment assignments in SPB Design
#' @param x X vector of covariate values, ranging from 0 to total number of factor levels
#' @param n Total sample size
#' @param m The number of blocks
#' @export
#' @examples
#' spbd()

spbd=function(x,n,m=4)
{tr=rep(NA,n)
i=1
while(i<=m)
{
  tr[x==i]=pbr(length(x[x==i]))
  i=i+1
}
return(tr)
}
