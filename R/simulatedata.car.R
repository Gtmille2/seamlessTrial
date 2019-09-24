#' Simulate Data For CAR Function
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme
#' @param mean.s Mean of the short term response
#' @param mean.t Mean long-term response
#' @param sigma0 sigma0 in the bivariate normal distribution
#' @param sigma sigma in bivariate normal distribution
#' @param rho rho in bivariate normal distribution
#' @param tau1 Tau1 is the weight on the first covariate
#' @param tau2 Tau2 is the weight on the second covariate
#' @param treat List containing the treatment assignments and the base matrix
#' @param covValues covValues is the dummy variable matrix from the covariate values
#' @param data Is the data from the previous analysis if applicable. The default is NULL
#' @export


simulatedata.car  <- function(mean.s,mean.t,sigma0 = 1,sigma=1,rho = .5, tau1 = 1, tau2 = 1,treat,covValues,data=NULL)
  # simulate multiple treatment groups
{

  trts = unique(treat)
  trts = trts[order(trts)]
  if (is.null(data)) keeps = rep(TRUE, nrow(covValues)) else keeps = data$treat %in% trts
  if (is.null(data)) covValues = covValues else covValues = covValues[(length(data$treat)+1):nrow(covValues),]
  n.trt = length(unique(treat))-1
  if (is.null(data)) n.s = table(treat) else n.s = table(treat)- table(data[data$treat %in% trts,]$treat)
  n.t = ceiling(n.s/2)
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)
  treat = tail(treat,sum(n.s))
  for (trt in trts)
  {
    mean.full.t[treat==trt][1:n.s[which(trts==trt)]] = c(rep(mean.t[which(trts==trt)], n.t[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])
 }

  mean.full = cbind(mean.full.s,mean.full.t)

  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  covmatrix = matrix(c(tau1*covValues[,1],tau2*covValues[,2]),ncol=2)
  full.data  = full.error + mean.full + covmatrix
  new.data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
  rbind(data[keeps,],new.data)
}

