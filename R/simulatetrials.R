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
simulatetrials <- function(N1=200, N=500, n.trt=3, mean.s=NULL, mean.t=NULL, p1 = .5, p2 = .5, sigma0=1, sigma=1, rho=0.5, nsim=1000, design = "Pocock",tau1 = 1,tau2 = 1)
{
  Norig =N
  N1orig = N1

  data = NULL
  z.v = NULL
  treat = NULL
  covValues = NULL
  n.looks = 2 #There are 2 planned analyses in this trial
  alpha.star.u <<- c(0.0,0.025)
  alpha.star.l <<- c(0,0.975)
  if (is.null(mean.s)) mean.s = rep(0,n.trt+1)
  if (is.null(mean.t)) mean.t = rep(0,n.trt+1)
  selected=rep(0,nsim)
  reject = rep(0,nsim)
  trialprogress <<- seq(1,n.looks)/n.looks
  for (sim in seq(1,nsim))
    {
    N1 = N1orig*(n.trt+1)
    N = Norig
    look = 1
    covValues = genCovValues(p = c(p1,p2),N=N1)
    #Getting treatment value assignments
    if (design == "Pocock" ) treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = n.trt) else treat = spbd(covValues = covValues, m = 4, best=0, tr = NULL, n.trt = n.trt)
    # table(treat)
    #Simulating data for these treatment assignments
    data = simulatedata.car(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2, treat, covValues,inspection = look)

    #Calculating test statistics z & v for this data, and selecting the best treatment
    z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
    z = z.v[1:look,1]
    v = z.v[1:look,2]
    z.v
    best = z.v[1,3] #The best treatment was found in this function
    # boundaries = get.boundaries(n.looks = look,v = v, k = c(n.trt,rep(1,look-1)), alpha.star.u, alpha.star.l) # Getting stopping boundaries at this point

    #Generating new covariate values
    look = 2
    left = sum(data$treat %in% c(0,best))
    N = N*2 - left
    covValuesNew = genCovValues(p=c(0.5,0.5),N=N)
    covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
    if (design == "Pocock" ) treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) else treat = spbd(covValues = covValues, m = 4, best=best, tr = treat, n.trt = n.trt) # Assigning new treatment values
    # print(TRUE)
    data = simulatedata.car(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data,inspection = look) #Simulating the new patients data

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


#' Testing other simulate trials function
#'
#' This function doesn't use an actual second point, just the projected based on the initial point
#' simulate.trials <- function(n1=20, N1=200, N=200, n.trt=3, mean.s=rep(0,3), mean.t=rep(0,3), p1, sigma0=1, sigma=1, rho=0.0, nsim=1000, save.boundary, simonly=0)
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
simulateest = function(N1=200, N=500, n.trt=3, mean.s=NULL, mean.t=NULL,p1, p2 , sigma0, sigma, rho, nsim, design = "Pocock",tau1,tau2,save.boundary)
{
  Norig =N
  N1orig = N1
  data = NULL
  z.v = NULL
  treat = NULL
  covValues = NULL
  n.looks = 2 #There are 2 planned analyses in this trial
  alpha.star.u <<- c(0.0,0.025)
  alpha.star.l <<- c(0.0,0.975)
  if (is.null(mean.s)) mean.s = rep(0,n.trt+1)
  if (is.null(mean.t)) mean.t = rep(0,n.trt+1)
  selected=rep(0,nsim)
  reject = rep(0,nsim)
  trialprogress <<- seq(1,n.looks)/n.looks
  for (sim in seq(1,nsim))
  {
    N1 = N1orig*(n.trt+1)
    N = Norig
    # data = simulate.data(mean.s,mean.t,N,sigma0,sigma,rho)
    look = 1
    covValues = genCovValues(p = c(p1,p2),N=N1)
    #Getting treatment value assignments
    if (design == "Pocock" ) treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = n.trt) else treat = spbd(covValues = covValues, m = 4, best=0, tr = NULL, n.trt = n.trt)
    table(treat)
    #Simulating data for these treatment assignments
    data = simulatedata.car(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2, treat, covValues,inspection = look)

    #Calculating test statistics z & v for this data, and selecting the best treatment
    z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
    best = z.v[1,3] #The best treatment was found in this function
    selected[sim] = best
    z = z.v[1:look,1]
    v = z.v[1:look,2]
    look = 2
    left = sum(data$treat %in% c(0,best))
    N = N*2 - left
    covValuesNew = genCovValues(p=c(p1,p2),N=N)
    covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
    if (design == "Pocock" ) treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) else treat = spbd(covValues = covValues, m = 4, best=best, tr = treat, n.trt = n.trt) # Assigning new treatment values
    # print(TRUE)
    data = simulatedata.car(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data,inspection = look) #Simulating the new patients data

    #Calculating test statistics z & v for this data at the second look.
    z.v = get.z.v.Current(data,n.looks,look,z.v)
    z = z.v[1:look,1] # Retrieving the Z statistic
    v = z.v[1:look,2] # Retrieving the V statistic

    t1percent = min(99,round(100*v[1]/v[2]))
    boundary.value = sqrt(v[2])*save.boundary[t1percent]
    if (z[2] > boundary.value) reject[sim] = 1
      # # using short-term and long-term endpoints
      # z.v = get.z.v(data,n1,N1,N)
      # print(z.v)
      # selected[sim] = seq(1,3)[z.v$v2>0]
      # t1percent = min(99,round(100*z.v$v1[1]/z.v$v2[selected[sim]]))
      # boundary.value <- sqrt(z.v$v2[selected[sim]])*save.boundary[t1percent]
      # if (z.v$z2[selected[sim]]>boundary.value) reject[sim] = 1


  }
  data.frame(power=sum(reject)/nsim,power3=sum(reject[selected==3])/nsim)
}

#' Simualte for no CAR
#'
#' This function doesn't use an actual second point, just the projected based on the initial point
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
#' @param nsim The number of simulation runs. The default is 10000
#' @param tau1 The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @param tau2 The chosen covariate adaptive randomization procedure. Default is Pocock's design
#' @export
simulateNoCar = function(N1=250, N=500, n.trt=3, mean.s=NULL, mean.t=NULL, p1 = .5, p2 = .5, sigma0=1, sigma=1, rho=0.5, nsim=10000,tau1 = 1,tau2 = 1,save.boundary)
{
  Norig =N
  N1orig = N1

  set.seed(10101)
  data = NULL
  z.v = NULL
  treat = NULL
  covValues = NULL
  n.looks = 2 #There are 2 planned analyses in this trial
  alpha.star.u <<- c(0.0,0.025)
  alpha.star.l <<- c(0.0,0.975)
  if (is.null(mean.s)) mean.s = rep(0,n.trt+1)
  if (is.null(mean.t)) mean.t = rep(0,n.trt+1)
  selected=rep(0,nsim)
  reject = rep(0,nsim)
  trialprogress <<- seq(1,n.looks)/n.looks
  for (sim in seq(1,nsim))
  {
    N1 = N1orig*(n.trt+1)
    N = Norig
    # data = simulate.data(mean.s,mean.t,N,sigma0,sigma,rho)
    look = 1
    covValues = genCovValues(p = c(p1,p2),N=N1)
    #Getting treatment value assignments
    treat = norandom(covValues = covValues, best =  0 , tr = NULL, n.trt = n.trt)
    # table(treat)
    #Simulating data for these treatment assignments
    data = simulatedata.nocar(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2, treat=treat, covValues=covValues,inspection = look,data = NULL)

    #Calculating test statistics z & v for this data, and selecting the best treatment
    z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
    best = z.v[1,3] #The best treatment was found in this function
    selected[sim] = best
    z = z.v[1:look,1]
    v = z.v[1:look,2]
    look = 2
    left = sum(data$treat %in% c(0,best))
    N = N*2 - left
    covValuesNew = genCovValues(p=c(p1,p2),N=N)
    covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
    treat = norandom(covValues = covValues, best =  best , tr = treat, n.trt = 1)
    # print(TRUE)
    data = simulatedata.nocar(mean.s = mean.s, mean.t = mean.t, sigma = sigma, sigma0 = sigma0, rho = rho, tau1 = tau1, tau2 = tau2,treat,covValues,data,inspection = look) #Simulating the new patients data

    #Calculating test statistics z & v for this data at the second look.
    z.v = get.z.v.Current(data,n.looks,look,z.v)
    z = z.v[1:look,1] # Retrieving the Z statistic
    v = z.v[1:look,2] # Retrieving the V statistic

    t1percent = min(99,round(100*v[1]/v[2]))
    boundary.value = sqrt(v[look])*save.boundary[t1percent]
    if (z[look] > boundary.value) reject[sim] = 1
    # # using short-term and long-term endpoints
    # z.v = get.z.v(data,n1,N1,N)
    # print(z.v)
    # selected[sim] = seq(1,3)[z.v$v2>0]
    # t1percent = min(99,round(100*z.v$v1[1]/z.v$v2[selected[sim]]))
    # boundary.value <- sqrt(z.v$v2[selected[sim]])*save.boundary[t1percent]
    # if (z.v$z2[selected[sim]]>boundary.value) reject[sim] = 1


  }
  data.frame(power=sum(reject)/nsim,power3=sum(reject[selected==3])/nsim)
}


