#' Function to create the boudnaries that satisfy the spending function
#'
#' calculate and store boundary values
#' if no stopping, single critical value/sqrt(v2) just depend on v1/v2
#' @export

save.boundary.values <- function()
{
  # calculate and store boundary values
  # if no stopping, single critical value/sqrt(v2) just depend on v1/v2
  save.boundary = rep(0,99)
  for (ratio in seq(0.01,0.99,0.01))
  {
    save.boundary[round(ratio*100)]=get.boundaries(v=c(ratio,1))$upper[2]
  }
  save.boundary
}
