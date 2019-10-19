#' Simulate trials function
#'
#' Simulates 2 look trial and calculates the overall type I error rate, and power.
#' @param n1 Number of patients with primary endpoint available at first analysis
#' @param N1 Number of patients with secondary endpoint available at first analysis
#' @param N The total number of patients in the trial
#' @param n.trt The number of treatments in the trial
#' @param mean.s The mean for short term endpoint sample groups
#' @param mean.t Mean for the long term endpoint sample groups
#' @param p1 The covariate value for the first covariate
#' @param p2 THe covariate value for the second covariate
#' @param sigma0 sigma0 in the bivariate normal distribution
#' @param sigma is the known sigma for the population
#' @param rho is the known correlation between endpoints.
#' @param nsim The number of simulation runs. The default is 1000
#' @param save.boundary The boundaries that satisfy the spending funciton
#' @param design The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @export

simulatetrials <- function(n1=20, N1=200, N=500, n.trt=2, mean.s=rep(0,3), mean.t=rep(0,3), p1 = ps[i,1], p2 = ps[i,2], sigma0=1, sigma=1, rho=0.0, nsim=1000, save.boundary, design = "Pocock",tau1,tau2)
{

  selected=rep(0,nsim)
  reject = rep(0,nsim)
  selected.t=rep(0,nsim)
  reject.t = rep(0,nsim)
  selected.t2=rep(0,nsim)
  reject.t2 = rep(0,nsim)
  selected.ks=rep(0,nsim)
  reject.ks = rep(0,nsim)

  for (sim in seq(1,nsim))
  {

    data = simulatedata.car(mean.s,mean.t, N,sigma0,sigma,p1,p2,rho,design = "Pocock",tau1 = tau1,tau2 = tau2)

    # using short-term and long-term endpoints
    z.v = get.z.v(data,n1,N1,N)
    selected[sim] = seq(1,2)[z.v$v2>0]
    t1percent = min(99,round(100*z.v$v1[1]/z.v$v2[selected[sim]]))
    boundary.value <- sqrt(z.v$v2[selected[sim]])*save.boundary[t1percent]
    if (z.v$z2[selected[sim]]>boundary.value) reject[sim] = 1


  }
  return((data.frame(power=sum(reject)/nsim,power2=sum(reject[selected==2])/nsim)))
}
