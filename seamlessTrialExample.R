### Example file for seamlessTrials package
### This package corresponds to paper II
### A Confirmatory seamless phase II/phase III clinical trial design incorporating short term endpoint information

library(seamlessTrials)
set.seed(101010)

n.looks = 3 #There are 5 planned analyses in this trial
trialprogress = seq(1,n.looks)/n.looks

alpha.star.u=c(0.001,0.010,0.019,0.024,0.025) # Upper alpha levels based on spending function
alpha.star.l=c(0.185,0.733,0.933,0.973,0.975) # Lower alpha levels based on spending function
mean.s = c(0,0,0) #The mean value for the short-term endpoint for the two treatments and control
mean.t = c(0,0,0) #The mean value for the long-term ednpoint for the two treatments and control
sigma0 = 1
sigma = 1
rho = 0.5
#Generating sample covariate values
covValues = genCovValues(p = c(0.5,0.5),N=100)
#Getting treatment value assignments
treat = psd(covValues, p1 = 3/4, best = 0, tr = NULL, n.trt = 2)
table(treat)
#Simulating data for these treatment assignments
look = 1
data = simulatedata.car(mean.s = c(0,0,0), mean.t = c(0,0,0), sigma = 1, sigma0 = 1, rho = 0.5, tau1 = 1, tau2 = 1, treat, covValues,inspection = look)



#Calculating test statistics z & v for this data, and selecting the best treatment
z.v = get.z.v.Current(data,n.looks,look,z.v.prev=NULL)
z = z.v[1:look,1]
v = z.v[1:look,2]
z.v
best = z.v[1,3] #The best treatment was found in this function

get.boundaries(n.looks = look,v = v, k = c(2,rep(1,look-1)), alpha.star.u, alpha.star.l) # Getting stopping boundaries at this point

#Generating new covariate values
look = 2
covValuesNew = genCovValues(p=c(0.5,0.5),N=100)
covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) # Assigning new treatment values
table(treat)
data = simulatedata.car(mean.s=c(0,0,0),mean.t=c(0,0,0),sigma = 1,sigma0=1,rho = 0.5, tau1 = 1,tau2 = 1,treat,covValues,data,inspection = look) #Simulating the new patients data

#Calculating test statistics z & v for this data at the second look.
z.v = get.z.v.Current(data,n.looks,look,z.v)

z = z.v[1:look,1] # Retrieving the Z statistic
v = z.v[1:look,2] # Retrieving the V statistic
get.boundaries(n.looks = look,v = v, k = c(2,rep(1,look-1)),alpha.star.u = alpha.star.u, alpha.star.l = alpha.star.l) # Calculating stopping boundaries at this look
z.v

look = 3
covValuesNew = genCovValues(p=c(0.5,0.5),N=100)
covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) # Assigning new treatment values
data = simulatedata.car(mean.s=c(0,0,0),mean.t=c(0,0,0),sigma = 1,sigma0=1,rho = 0.5, tau1 = 1,tau2 = 1,treat,covValues,data = data,inspection = look) #Simulating the new patients data

for ( look in 3:n.looks) {

  covValuesNew = genCovValues(p=c(0.5,0.5),N=100)
  covValues = rbind(covValues, covValuesNew) # Combining new covariate values with old covariate values
  treat = psd(covValues,p1=3/4,best = best,tr = treat, n.trt = 1) # Assigning new treatment values
  data = simulatedata.car(mean.s=c(0,0,0),mean.t=c(0,0,0),sigma = 1,sigma0=1,rho = 0.5, tau1 = 1,tau2 = 1,treat,covValues,data = data,inspection = look) #Simulating the new patients data

  z.v = get.z.v.Current(data,n.looks,look,z.v)
  z = z.v[1:look,1] # Retrieving the Z statistic
  v = z.v[1:look,2] # Retrieving the V statistic
  boundaries = get.boundaries(n.looks = look,v = v, k = c(2,rep(1,look-1)),alpha.star.u = alpha.star.u, alpha.star.l = alpha.star.l) # Calculating stopping boundaries at this look

}
boundaries

