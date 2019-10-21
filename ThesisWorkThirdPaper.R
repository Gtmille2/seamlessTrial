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
covValues = genCovValues(p = c(0.5,0.5),N=100)
#Getting treatment value assignments
treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = 2)
table(treat)
#Simulating data for these treatment assignments
data = simulatedata.car(mean.s = c(0,0,0), mean.t = c(0,0,0), sigma = 1, sigma0 = 1, rho = 0.5, tau1 = 1, tau2 = 1, treat, covValues)

#Calculating test statistics z & v for this data, and selecting the best treatment
look = 1 # This is the first look
z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
z = z.v[1:look,1]
v = z.v[1:look,2]
z.v
best = z.v[1,3] #The best treatment was found in this function

get.boundaries(n.looks = look,v = v, k = c(2,rep(1,look-1)), alpha.star.u, alpha.star.l) # Getting stopping boundaries at this point

simulatetrials <- function(n1=20, N1=200, N=100, n.trt=2, mean.s=rep(0,3), mean.t=rep(0,3), p1 = ps[i,1], p2 = ps[i,2], sigma0=1, sigma=1, rho=0.0, nsim=1000, save.boundary, design = "Pocock",tau1,tau2)
{

  selected=rep(0,nsim)
  reject = rep(0,nsim)

  for (sim in seq(1,nsim))
  {
    #Incorporating the Pocock Design and SPBD design to calculate the number of trials that exceed the boundaries
    #Will do so for 2 looks
    n.trt = 2
    covValues = genCovValues(p = c(0.5,0.5),N=100)
    #Getting treatment value assignments
    treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = n.trt )
    #Simulating data for these treatment assignments
    data = simulatedata.car(mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), sigma = 1, sigma0 = 1, rho = 0.5, tau1 = 1, tau2 = 1, treat, covValues)

    #Calculating test statistics z & v for this data, and selecting the best treatment
    look = 1 # This is the first look
    z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
    z = z.v[1:look,1]
    v = z.v[1:look,2]
    best = z.v[1,3]
    selected[sim] = best
    #The best treatment was found in this function
    boundaries = get.boundaries(n.looks = look,v = v, k = c(n.trt,rep(1,look-1)), alpha.star.u, alpha.star.l) # Getting stopping boundaries at this point
    #No stopping at first look
    covValuesNew = genCovValues(p=c(0.5,0.5),N=100)
    covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
    treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) # Assigning new treatment values
    table(treat)
    data = simulatedata.car(mean.s=rep(0,n.trt+1),mean.t=rep(0,n.trt+1),sigma = 1,sigma0=1,rho = 0.5, tau1 = 1,tau2 = 1,treat,covValues,data) #Simulating the new patients data
    table(treat)

    #Calculating test statistics z & v for this data at the second look.
    look = 2
    z.v = get.z.v.Current(data,n.looks,look,z.v)

    z = z.v[1:look,1] # Retrieving the Z statistic
    v = z.v[1:look,2] # Retrieving the V statistic
    boundaries = get.boundaries(n.looks = look,v = v, k = c(n.trt,rep(1,look-1)),alpha.star.u = alpha.star.u, alpha.star.l = alpha.star.l) # Calculating stopping boundaries at this look

    if (z.v$z[look] > boundaries$upper[look]) reject[sim] = 1
    # # using short-term and long-term endpoints
    # z.v = get.z.v(data,n1,N1,N)
    # selected[sim] = seq(1,2)[z.v$v2>0]
    # t1percent = min(99,round(100*z.v$v1[1]/z.v$v2[selected[sim]]))
    # boundary.value <- sqrt(z.v$v2[selected[sim]])*save.boundary[t1percent]
    # if (z.v$z2[selected[sim]]>boundary.value) reject[sim] = 1


  }
  # return((data.frame(power=sum(reject)/nsim,power2=sum(reject[selected==2])/nsim)))
}
