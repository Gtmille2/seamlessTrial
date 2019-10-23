#' Simulate trials function
#'
#' Simulates 2 look trial and calculates the overall type I error rate.
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
#' @param design The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @param tau1 The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @param tau2 The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @export
simulatetrials <- function(n1=20, N1=200, N=500, n.trt=2, mean.s=rep(0,3), mean.t=rep(0,3), p1 = ps[i,1], p2 = ps[i,2], sigma0=1, sigma=1, rho=0.0, nsim=1000, save.boundary, design = "Pocock",tau1,tau2)
{


  library(seamlessTrials)
  set.seed(101010)
  n.looks = 2 #There are 5 planned analyses in this trial
  alpha.star.u=c(0.001,0.010,0.019,0.024,0.025) # Upper alpha levels based on spending function
  alpha.star.l=c(0.185,0.733,0.933,0.973,0.975) # Lower alpha levels based on spending function
  mean.s = c(0,0,0) #The mean value for the short-term endpoint for the two treatments and control
  mean.t = c(0,0,0) #The mean value for the long-term ednpoint for the two treatments and control
  sigma0 = 1
  sigma = 1
  rho = 0.5
  trialprogress = seq(1,n.looks)/n.looks
  for (sim in seq(1,nsim)) {

  look = 1
  covValues = genCovValues(p = c(p1,p2),N=N1)
  #Getting treatment value assignments
  if (design == "Pocock" ) treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = n.trt) else treat = spbd(covValues = covValues, m = 4, best=0, tr = NULL, n.trt = n.trt)
  # table(treat)
  #Simulating data for these treatment assignments
  data = simulatedata.car(mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), sigma = 1, sigma0 = 1, rho = 0.5, tau1 = 1, tau2 = 1, treat, covValues,inspection = look)

  #Calculating test statistics z & v for this data, and selecting the best treatment
  z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
  z = z.v[1:look,1]
  v = z.v[1:look,2]
  # z.v
  best = z.v[1,3] #The best treatment was found in this function
  # boundaries = get.boundaries(n.looks = look,v = v, k = c(n.trt,rep(1,look-1)), alpha.star.u, alpha.star.l) # Getting stopping boundaries at this point

  #Generating new covariate values
  look = 2
  covValuesNew = genCovValues(p=c(0.5,0.5),N=N-N1)
  covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
  if (design == "Pocock" ) treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) else treat = spbd(covValues = covValues, m = 4, best=0, tr = NULL, n.trt = n.trt) # Assigning new treatment values
  # print(TRUE)
  data = simulatedata.car(mean.s=rep(0,n.trt+1),rep(0,n.trt+1),sigma = 1,sigma0=1,rho = 0.5, tau1 = 1,tau2 = 1,treat,covValues,data,inspection = look) #Simulating the new patients data

  #Calculating test statistics z & v for this data at the second look.
  z.v = get.z.v.Current(data,n.looks,look,z.v)
  z = z.v[1:look,1] # Retrieving the Z statistic
  v = z.v[1:look,2] # Retrieving the V statistic
  boundaries = get.boundaries(n.looks = look,v = v, k = c(n.trt,rep(1,look-1)),alpha.star.u = alpha.star.u, alpha.star.l = alpha.star.l) # Calculating stopping boundaries at this look

  selected[sim] = best
  if (z.v$z[look] > boundaries$upper[look]) reject[sim] = 1


  }
  reject
}
