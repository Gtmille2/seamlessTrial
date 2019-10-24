##### This script is to test different parameter settings for the third paper to determine how they affect the outcome of the type I error #####


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
#Generating sample covariate values
simulatetrials <- function(N1=200, N=500, n.trt=3, mean.s=rep(0,n.trt+1), mean.t=rep(0,n.trt+1), p1 = .5, p2 = .5, sigma0=1, sigma=1, rho=0.5, nsim=1000, design = "Pocock",tau1 = 1,tau2 = 1)
{

  selected=rep(0,nsim)
  reject = rep(0,nsim)

  for (sim in seq(1,nsim))
  {
    #Incorporating the Pocock Design and SPBD design to calculate the number of trials that exceed the boundaries
    #Will do so for 2 looks
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
    # # using short-term and long-term endpoints
    # z.v = get.z.v(data,n1,N1,N)
    # selected[sim] = seq(1,2)[z.v$v2>0]
    # t1percent = min(99,round(100*z.v$v1[1]/z.v$v2[selected[sim]]))
    # boundary.value <- sqrt(z.v$v2[selected[sim]])*save.boundary[t1percent]
    # if (z.v$z2[selected[sim]]>boundary.value) reject[sim] = 1


  }
  reject
}
#Running multiple simulations:

n.trt = 3
nsim = 100
allsims = NULL
###Under null hypothesis that all treatment effects are equal
ptm = proc.time()
reject = simulatetrials(N1 = 200, N = 500, n.trt = n.trt, mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 100, design = "Pocock",tau1= 1,tau2 = 1)
proc.time() - ptm
ptm = proc.time()
reject = simulatetrials(nsim=100)
proc.time() - ptm
allsims = cbind(allsims, reject)
reject = simulatetrials(N1 = 300, N = 500, n.trt = n.trt, mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock", tau1= 1,tau2 = 1)
allsims = cbind(allsims, reject)
proc.time() - ptm
###Under alternative hypothesis that at least 1 treatment is not equal to zero
ptm = proc.time()
reject = simulatetrials(N1 = 200, N = 500, n.trt = n.trt, mean.s = rep(0,(n.trt+1)), mean.t = rep(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 100, design = "SPBD", tau1= 1,tau2 = 1)
proc.time() - ptm

reject = simulatetrials(N1 = 300, N = 500, n.trt = n.trt, mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "SPBD", tau1= 1,tau2 = 1)

