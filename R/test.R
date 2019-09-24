#' Simulate Data For CAR Function
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme
#' @param mean.s Mean of the short term response
#' @param mean.t Mean long-term response
#' @param N Total Sample Size
#' @param sigma0 sigma0 in the bivariate normal distribution
#' @param sigma sigma in bivariate normal distribution
#' @param p1 Probability of success in first covariate
#' @param p2 Probability of success in second covariate
#' @param rho rho in bivariate normal distribution
#' @param design Covariate Adapative Randomization scheme. Default is Pocock
#' @param tau1 Tau1 is the weight on the first covariate
#' @param tau2 Tau2 is the weight on the second covariate
#' @export


simulateData  <- function(mean.s,mean.t,N=500,sigma0,sigma,p1,p2,rho, design = "Pocock", tau1, tau2)
  # simulate multiple treatment groups
{
  p = 3/4
  n.trt = length(mean.s)-1
  z1=sample(2,N,p1)
  z2=sample(2,N,p2)
  all = matrix(c(z1,z2),ncol=2)
  if (design == "Pocock") treat = psd(all,N,p=3/4,n.trt) else treat = spbd(x,N)
  n.s = table(treat)
  n.t = ceiling(n.s/2)
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)

  for (trt in seq(0,n.trt))
  {

    mean.full.t[treat==trt] = c(rep(mean.t[trt+1], n.t[trt+1]),rep(NA,diff[trt+1]))
    mean.full.s[treat==trt] = rep(mean.s[trt+1], n.s[trt+1])
    # mean.full.s[treat==trt] = rep(mean.s[trt+1], sum(treat == trt))
    # mean.full.t[treat==trt] = rep(mean.t[trt+1], sum(treat==trt))
  }

  mean.full = cbind(mean.full.s,mean.full.t)

  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm((n.trt)/2*N,mean,var)
  full.data  = full.error + mean.full + tau1*z1 + tau2*z2
  data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
}
