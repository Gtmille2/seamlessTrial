#' Function to create the boudnaries that satisfy the spending function
#'
#' @param n.looks The number of interim analyses occuring in the trial
#' @param v Ratio used in the calculation of boundaries
#' @param k Number of treatments used in the trial
#' @param alpha.star.u the upper limit alpha on boundaries for the spending function
#' @param alpha.star.l the lower limit alpha on boundaries for the spending function
#' @export


get.boundaries <- function(n.looks=2, v, k=c(3,1), alpha.star.u=c(0,0.025), alpha.star.l=c(0,0.975))
{

  inc.v = get.increments(v)
  inc.alpha.star.u = get.increments(alpha.star.u)
  inc.alpha.star.l = get.increments(alpha.star.l)

  n.grid = 1001
  upper = rep(0,n.looks)
  lower = rep(0,n.looks)

  for (look in seq(1,n.looks))
  {

    # find boundaries for this look using spending functions unless both spending function increments are 0


    if (look==1) alpha.spent=0
    else alpha.spent=alpha.star.u[look-1]+alpha.star.l[look-1]

    if (look ==1)
    {
      # these values are not used but need to be set to get function call correct
      upper.prev=0
      lower.prev=0
      density.prev=0
    }
    else
    {
      upper.prev = upper[look-1]
      lower.prev = lower[look-1]
    }

    if (inc.alpha.star.u[look] == 0) upper[look] = 10*sqrt(v[look])
    else upper[look] = uniroot(distribution.search.function, lower=-10*sqrt(v[look]), upper=10*sqrt(v[look]), look=look, k=k[look], incv=inc.v[look], upper.prev=upper.prev, lower.prev=lower.prev, n.grid=n.grid, density.prev=density.prev, target=(1-alpha.spent-inc.alpha.star.u[look]))$root

    if (inc.alpha.star.l[look] ==0) lower[look] = -10*sqrt(v[look])
    else lower[look] = uniroot(distribution.search.function, lower=-10*sqrt(v[look]), upper=10*sqrt(v[look]), look=look, k=k[look], incv=inc.v[look], upper.prev=upper.prev, lower.prev=lower.prev, n.grid=n.grid, density.prev=density.prev, target=inc.alpha.star.l[look])$root



    # set up grid and density for next look

    if (look < n.looks)
    {
      z.grid = seq(lower[look],upper[look],(upper[look]-lower[look])/(n.grid-1))

      dens = rep(0,n.grid)
      if (look==1)
        for (i in seq(1,n.grid))
          dens[i] = get.max.density(k[look], inc.v[look], z.grid[i])
      else
      {
        #  integrate to give convolution of current density function with previous density
        z.grid.prev = seq(lower[look-1],upper[look-1],(upper[look-1]-lower[look-1])/(n.grid-1))
        values = rep(0,n.grid)
        for (i in seq(1,n.grid))
        {
          for (j in seq(1,n.grid))
            values[j] = get.max.density(k[look],inc.v[look],z.grid[i]-z.grid.prev[j])*density.prev[j]
          integral = (values[1] + 2*sum(values[2:(n.grid-1)]) + values[n.grid])*(upper[look-1]-lower[look-1])/(2*(n.grid-1))
          dens[i] = integral
        }
      }
      # store density from previous look to use in integration
      density.prev = dens
    }

  }

  # stopping boundaries
  boundaries <- data.frame(lower,upper)

  boundaries
}
