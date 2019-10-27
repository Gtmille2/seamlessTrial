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
simulatedata.car  <- function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data=NULL,inspection)
  # simulate multiple treatment groups
{
  # print(head(covValues))
  if (is.null(data)) test = simnew(mean.s = mean.s,mean.t = mean.t,sigma0 = sigma0,sigma=sigma,rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data=NULL,inspection) else simold(mean.s,mean.t,sigma0 = sigma0,sigma=sigma,rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data,inspection)

}
#' Simulate Data For CAR at first look
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme at first look
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
simnew = function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data= NULL,inspection) {
  prop = trialprogress[inspection]
  trts = unique(treat)
  trts = trts[order(trts)]
  keeps = rep(TRUE, nrow(covValues))
  keepdata = data[keeps,]
  nas = is.na(data[keeps,]$t)
  covkeep = NULL
  oldcov = NULL
  covnew = covValues
  covValues <<- rbind(covkeep,covnew)
  term1 = sum(!is.na(data[keeps,]$t))
  term2 = prop*length(treat)
  p = 0
  oldupdate = 0
  # oldupdate = old
  newupdate = length(treat)*prop
  # oldupdate = round(p*(term2-term1))
  # newupdate = round(q*(term2-term1))

  n.trt = length(unique(treat))-1
  n.s = table(treat)
  #Distributing the new available long-term endpoints between the previous group and the incoming group
  # p = (length(treat) - sum(!is.na(data$t)))*(length(data$s)/length(treat)) # The number of patients from the old data that will have their long-term endpoint updated
  # q = (length(treat) - sum(!is.na(data$t)))*(length(data$s)/length(treat)) # The number of patients from the new data that will have long-term endpoints
  n.t = round(n.s*(newupdate/sum(n.s)))
  # n.update = table(data[keeps,]$treat) - table(na.omit(data[keeps,])$treat)
  # n.update = round(n.update*(oldupdate/sum(n.update)))
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)
  treat = tail(treat,sum(n.s))
  # oldtreat = keepdata[is.na(keepdata$t),]$treat
  # olddiff = table(oldtreat)-n.update
  # if (sum(olddiff) < 0) olddiff[1:length(olddiff)] = 0
  # mean.new.t = rep(0,length(oldtreat))

  for (trt in trts)
  {
    mean.full.t[treat==trt][1:n.s[which(trts==trt)]] = c(rep(mean.t[which(trts==trt)], n.t[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])
  }
  ### Have to add on
  mean.full = cbind(mean.full.s,mean.full.t)
  # update.full = cbind(mean.new.t, rep(0,length(mean.new.t)))
  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  # update.error = NULL
  covsum = rowSums(matrix(c(tau1*covnew[,1],tau2*covnew[,2]),ncol=2))
  covmatrix = cbind(covsum, covsum)

  full.data  = full.error + mean.full + covmatrix
  new.data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
  test = new.data
  test
}
#' Simulate Data For CAR after first look
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme after first look
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
simold = function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data,inspection) {
  prop = trialprogress[inspection]
  trts = unique(treat)
  trts = trts[order(trts)]

  keeps = data$treat %in% trts
  keepdata = data[keeps,]
  nas = is.na(data[keeps,]$t)
  covkeep  = covValues[1:length(data$s),][keeps,]
  oldcov = covValues[1:length(data$s),][keeps,][nas, ]
  covnew =  covValues[(length(data$treat)+1):nrow(covValues),]
  covValues <<- rbind(covkeep,covnew)
  term1 = sum(!is.na(data[keeps,]$t))
  term2 = prop*length(treat)
  p = nrow(data[keeps,])/length(treat)
  oldupdate = round(prop*nrow(data[keeps,])) - term1
  # oldupdate = old -
  newupdate = (prop*(length(treat)-nrow(data[keeps,])))
  # oldupdate = round(p*(term2-term1))
  # newupdate = round(q*(term2-term1))

  n.trt = length(unique(treat))-1
  n.s = table(treat)- table(data[data$treat %in% trts,]$treat)
  n.t = round(n.s*(newupdate/sum(n.s)))
  n.update = table(data[keeps,]$treat) - table(na.omit(data[keeps,])$treat)
  n.update = round(n.update*(oldupdate/sum(n.update)))
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)
  treat = tail(treat,sum(n.s))
  oldtreat = keepdata[is.na(keepdata$t),]$treat
  olddiff = table(oldtreat)-n.update
  if (sum(olddiff) < 0) olddiff[1:length(olddiff)] = 0
  mean.new.t = rep(0,length(oldtreat))

  for (trt in trts)
  {
    mean.full.t[treat==trt][1:n.s[which(trts==trt)]] = c(rep(mean.t[which(trts==trt)], n.t[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])
    mean.new.t[oldtreat==trt] = c(rep(mean.t[which(trts==trt)],n.update[which(trts==trt)]) , rep(NA,olddiff[which(trts==trt)]))
  }
  ### Have to add on
  mean.full = cbind(mean.full.s,mean.full.t)
  update.full = cbind(mean.new.t, rep(0,length(mean.new.t)))
  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  update.error = mvtnorm::rmvnorm(nrow(update.full),mean,var)
  covsum = rowSums(matrix(c(tau1*covnew[,1],tau2*covnew[,2]),ncol=2))
  covmatrix = cbind(covsum, covsum)
  updatecovsum = rowSums(matrix(c(tau1*oldcov[,1],tau2*oldcov[,2]),ncol=2))
  updatecovmat = cbind(updatecovsum, updatecovsum)
  updatedata = update.error + mean.new.t + updatecovmat
  olds = keepdata[is.na(keepdata$t),]$s
  updatedata = data.frame(s = olds, t = updatedata[,2],treat=oldtreat)
  olddata = na.omit(keepdata)
  updatefull = rbind(olddata, updatedata)

  full.data  = full.error + mean.full + covmatrix
  new.data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
  test = rbind(updatefull,new.data)
  test
}
#' Simulate Data For CAR at first look
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme at first look
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
simnewnocar = function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data= NULL,inspection) {
  prop = trialprogress[inspection]
  trts = unique(treat)
  trts = trts[order(trts)]

  keeps = rep(TRUE, nrow(covValues))
  keepdata = data[keeps,]
  nas = is.na(data[keeps,]$t)
  covkeep = NULL
  oldcov = NULL
  covnew = covValues
  covValues <<- rbind(covkeep,covnew)
  term1 = sum(!is.na(data[keeps,]$t))
  term2 = prop*length(treat)
  p = 0
  oldupdate = 0
  # oldupdate = old
  newupdate = length(treat)*prop
  # oldupdate = round(p*(term2-term1))
  # newupdate = round(q*(term2-term1))

  n.trt = length(unique(treat))-1
  n.s = table(treat)
  #Distributing the new available long-term endpoints between the previous group and the incoming group
  # p = (length(treat) - sum(!is.na(data$t)))*(length(data$s)/length(treat)) # The number of patients from the old data that will have their long-term endpoint updated
  # q = (length(treat) - sum(!is.na(data$t)))*(length(data$s)/length(treat)) # The number of patients from the new data that will have long-term endpoints
  n.t = round(n.s*(newupdate/sum(n.s)))
  # n.update = table(data[keeps,]$treat) - table(na.omit(data[keeps,])$treat)
  # n.update = round(n.update*(oldupdate/sum(n.update)))
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)
  treat = tail(treat,sum(n.s))
  # oldtreat = keepdata[is.na(keepdata$t),]$treat
  # olddiff = table(oldtreat)-n.update
  # if (sum(olddiff) < 0) olddiff[1:length(olddiff)] = 0
  # mean.new.t = rep(0,length(oldtreat))

  for (trt in trts)
  {
    mean.full.t[treat==trt][1:n.s[which(trts==trt)]] = c(rep(mean.t[which(trts==trt)], n.t[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])
  }
  ### Have to add on
  mean.full = cbind(mean.full.s,mean.full.t)
  # update.full = cbind(mean.new.t, rep(0,length(mean.new.t)))
  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  # update.error = NULL
  # covsum = rowSums(matrix(c(tau1*covnew[,1],tau2*covnew[,2]),ncol=2))
  # covmatrix = cbind(covsum, covsum)

  full.data  =  mean.full + full.error
  new.data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
  test = new.data
  test
}
#' Simulate Data For CAR after first look
#'
#' Function to simulate data with treatment assignments from Covariate Adaptive Randomization Scheme after first look
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
simoldnocar = function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data,inspection) {
  prop = trialprogress[inspection]
  trts = unique(treat)
  trts = trts[order(trts)]

  keeps = data$treat %in% trts
  keepdata = data[keeps,]
  nas = is.na(data[keeps,]$t)
  covkeep  = covValues[1:length(data$s),][keeps,]
  oldcov = covValues[1:length(data$s),][keeps,][nas, ]
  covnew =  covValues[(length(data$treat)+1):nrow(covValues),]
  covValues <<- rbind(covkeep,covnew)
  term1 = sum(!is.na(data[keeps,]$t))
  term2 = prop*length(treat)
  p = nrow(data[keeps,])/length(treat)
  oldupdate = round(prop*nrow(data[keeps,])) - term1
  # oldupdate = old -
  newupdate = (prop*(length(treat)-nrow(data[keeps,])))
  # oldupdate = round(p*(term2-term1))
  # newupdate = round(q*(term2-term1))

  n.trt = length(unique(treat))-1
  n.s = table(treat)- table(data[data$treat %in% trts,]$treat)
  n.t = round(n.s*(newupdate/sum(n.s)))
  n.update = table(data[keeps,]$treat) - table(na.omit(data[keeps,])$treat)
  n.update = round(n.update*(oldupdate/sum(n.update)))
  diff = n.s-n.t
  mean.full.s = rep(rep(0,n.trt+1),n.s)
  mean.full.t = rep(rep(0,n.trt+1),n.s)
  treat = tail(treat,sum(n.s))
  oldtreat = keepdata[is.na(keepdata$t),]$treat
  olddiff = table(oldtreat)-n.update
  if (sum(olddiff) < 0) olddiff[1:length(olddiff)] = 0
  mean.new.t = rep(0,length(oldtreat))

  for (trt in trts)
  {
    mean.full.t[treat==trt][1:n.s[which(trts==trt)]] = c(rep(mean.t[which(trts==trt)], n.t[which(trts==trt)]),rep(NA,diff[which(trts==trt)]))
    mean.full.s[treat==trt] = rep(mean.s[which(trts==trt)], n.s[which(trts==trt)])
    mean.new.t[oldtreat==trt] = c(rep(mean.t[which(trts==trt)],n.update[which(trts==trt)]) , rep(NA,olddiff[which(trts==trt)]))
  }
  ### Have to add on
  mean.full = cbind(mean.full.s,mean.full.t)
  update.full = cbind(mean.new.t, rep(0,length(mean.new.t)))
  # simulate error about mean zero then add on treatment effects
  mean = c(0,0)
  var = c(sigma0^2,rho*sigma0*sigma,rho*sigma0*sigma,sigma^2)
  dim(var) = c(2,2)
  full.error = mvtnorm::rmvnorm(nrow(mean.full),mean,var)
  update.error = mvtnorm::rmvnorm(nrow(update.full),mean,var)

  updatedata = update.error + mean.new.t
  olds = keepdata[is.na(keepdata$t),]$s
  updatedata = data.frame(s = olds, t = updatedata[,2],treat=oldtreat)
  olddata = na.omit(keepdata)
  updatefull = rbind(olddata, updatedata)

  full.data  = full.error + mean.full
  new.data = data.frame(s=full.data[,1],t=full.data[,2],treat=treat)
  test = rbind(updatefull,new.data)
  test
}
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
simulatedata.nocar  <- function(mean.s,mean.t,sigma0,sigma,rho, tau1, tau2,treat,covValues,data=NULL,inspection)
  # simulate multiple treatment groups
{
  if (is.null(data)) test = simnewnocar(mean.s = mean.s,mean.t = mean.t,sigma0 = sigma0,sigma=sigma,rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data=NULL,inspection) else simoldnocar(mean.s,mean.t,sigma0 = 1,sigma=1,rho = .5, tau1 = 1, tau2 = 1,treat,covValues,data,inspection)

}
