##### This script is to test different parameter settings for the third paper to determine how they affect the outcome of the type I error #####


library(seamlessTrials)
#Running multiple simulations:
save.boundary = save.boundary.values()
n.trt = 3
nsim = 1000
n1 = 20
N1 = 100
N=200
sigma0 = 1
sigma = 1
p1 = .5
p2 = .5
rho = 0.5
tau1 = 1
tau2 = 1
mean.s = NULL
mean.t = NULL
set.seed(10101)
allsims = NULL
###Under null hypothesis that all treatment effects are equal
ptm = proc.time()
reject = simulatetrials(N1 = 250, N = 500, n.trt = n.trt, mean.s = rep(0,n.trt+1), mean.t = rep(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 100, design = "Pocock",tau1= 1,tau2 = 1)
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
###Using boundary.values to calculate the future boundary value:
n.trt = 3


#Pocock's design for Type II error
ptm = proc.time()
pocockreject1 = simulateest(N1 = 100, N = 200, n.trt = n.trt, mean.s = c(rep(0,3),0), mean.t =c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 100, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
pocockreject2 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),.333), mean.t = c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
pocockreject3 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),-.333), mean.t = c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
#### Pocock's design for Power
ptm = proc.time()
pocockpower1 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),0), mean.t = c(rep(0,3),0.333), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
pocockpower2 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),.333), mean.t = c(rep(0,3),0.333), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
pocockpower3 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),-.333), mean.t = c(rep(0,3),0.333)(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "Pocock",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm

### SPBD for type II error
ptm = proc.time()
spbdreject1 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),0), mean.t =c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 100, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
spbdreject2 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),.333), mean.t = c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
spbdreject3 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),-.333), mean.t = c(rep(0,3),0), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
### SPBD for power
ptm = proc.time()
spbdpower1 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),0), mean.t = c(rep(0,3),0.333), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
spbdpower2 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),.333), mean.t = c(rep(0,3),0.333), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
spbdpower3 = simulateest(N1 = 250, N = 500, n.trt = n.trt, mean.s = c(rep(0,3),-.333), mean.t = c(rep(0,3),0.333)(0,n.trt+1), p1 = .5,p2 = .5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000, design = "spbd",tau1= 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm

####No CAR
ptm = proc.time()
nocarreject1 = simulateNoCar(n1 = 20, N1 = 100,N = 200,n.trt = 3,mean.s  = c(rep(0,3),0),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = 0,nsim = 1000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarreject2 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0.333),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarreject3 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),-0.333),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower1 = simulateNoCar(n1 = 20, N1 = 100,N = 200,n.trt = 3,mean.s  = c(rep(0,3),0),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = 0,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower2 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0.333),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower3 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),-0.333),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm

nocarpower1 = simNoCarOrig(n1=20, N1= 100, N=200, n.trt =3, mean.s=rep(0,n.trt+1), mean.t = c(rep(0,3), 0.3333), p1 = .5, p2=.5, sigma0 = 1, sigma = 1, rho = 0.5, nsim = 1000, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
pocockpower1 = simCarOrig(n1= 20, N1= 100, N =200, n.trt = 3, mean.s = rep(0,n.trt+1), mean.t = c(rep(0,3), 0.3333), p1 = .5, p2 = 5, sigma0 = 1, sigma = 1, rho = .5, nsim = 100, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")


allsimerror = function(n1 = 20, N1 = 100, N = 200, n.trt = 3, p1, p2, sigma0 = 1, sigma = 1, rho, nsim = 10000, tau1 = 1, tau2 = 1, save.boundary) {
  pocockerror1 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  pocockerror2 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  pocockerror3 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  spbderror1 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  spbderror2 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  spbderror3 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  nocarerror1 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
  nocarerror2 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
  nocarerror3 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)

  pocockerror = cbind(pocockerror1,pocockerror2,pocockerror3)
  spbderror = cbind(spbderror1,spbderror2,spbderror3)
  nocarerror = cbind(nocarerror1,nocarerror2,nocarerror3)

  allsim = rbind(pocockerror, spbderror, nocarerror)
  return(allsim)


}

allsimpower = function(n1 = 20, N1 = 100, N = 200, n.trt = 3, p1, p2, sigma0 = 1, sigma = 1, rho, nsim = 10000, tau1 = 1, tau2 = 1, save.boundary) {
  pocockpower1 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  pocockpower2 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  pocockpower3 = simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
  spbdpower1 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  spbdpower2 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  spbdpower3 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
  nocarpower1 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
  nocarpower2 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
  nocarpower3 =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)

  pocockpower = cbind(pocockpower1,pocockpower2,pocockpower3)
  spbdpower = cbind(spbdpower1,spbdpower2,spbdpower3)
  nocarpower = cbind(nocarpower1,nocarpower2,nocarpower3)

  allsim = rbind(pocockpower, spbdpower, nocarpower)
  return(allsim)
}
#### Testing 3 different rho values ####
n1 = 20
N1 =100
N = 200
p1 = .5
p2 = .5
nsim = 5000
ptm = proc.time()
error_rho_0 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_5 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_7 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
proc.time() -ptm
allsimerror1 = rbind(error_rho_0,error_rho_5,error_rho_7)
write.csv(allsimerror1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimerror20_100_200_5000.csv")
proc.time() -ptm
power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower20_100_200_5000.csv")

n1 = 20
N1 =100
N = 200
p1 = .5
p2 = .4
nsim = 5000
ptm = proc.time()
error_rho_0 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_5 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_7 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
proc.time() -ptm
allsimerror1 = rbind(error_rho_0,error_rho_5,error_rho_7)
write.csv(allsimerror1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimerror20_100_200_5000_rho54.csv")
proc.time() -ptm
power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower20_100_200_5000_rho54.csv")

p1 = .5
p2 = .5
n1 = 50
N1 =150
N = 300
nsim = 5000
ptm = proc.time()
error_rho_0 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_5 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_7 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
proc.time() -ptm
allsimerror1 = rbind(error_rho_0,error_rho_5,error_rho_7)
write.csv(allsimerror1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimerror50_150_300_5000.csv")
proc.time() -ptm
power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower50_150_300_5000.csv")


#### Trying to change the block size to determine its effect on the type I error
n.trt = 3
n1 = 20
N1 =100
N = 200
nsim = 5000
rho = 0.5
p1 = .5
p2 = .5
spbderror1b12 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 12)
spbderror2b12 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 12)
spbderror3b12 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 12)
spbderror1b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbderror2b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbderror3b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbderror1b20 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 20)
spbderror2b20 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 20)
spbderror3b20 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 20)
spbderror1b40 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 40)
spbderror2b40 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 40)
spbderror3b40 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 40)

spbdpower1b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbdpower2b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbdpower3b16 =   simCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.3333), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd", block.size = 16)
spbdb12 = rbind(spbderror1b12,spbderror2b12,spbderror3b12)
spbdb16 = rbind(spbderror1b16,spbderror2b16,spbderror3b16)
spbdb16power = rbind(spbdpower1b16,spbdpower2b16,spbdpower3b16)
spbdb20 = rbind(spbderror1b20,spbderror2b20,spbderror3b20)
spbdb40 = rbind(spbderror1b40,spbderror2b40,spbderror3b40)

spbderror = cbind(spbdb12$power,spbdb16$power,spbdb20$power, spbdb40$power)
write.csv(spbdm4,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/spbdb12.csv")
write.csv(spbdm6,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/spbdb16.csv")
write.csv(spbdm6power,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/spbdb16power.csv")
colMeans(spbderror)


n1 = 100
N1 =200
N = 400
nsim = 5000
ptm = proc.time()
error_rho_0 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_5 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
error_rho_7 = allsimerror(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
proc.time() -ptm
allsimerror1 = rbind(error_rho_0,error_rho_5,error_rho_7)
write.csv(allsimerror1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimerror100_200_400_5000.csv")
proc.time() -ptm
power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower100_200_400_5000.csv")

##### Testing Bootstrap
n1 = 20
N1 =100
N = 200
p1 = .5
p2 = .5
nsim = 5000
n.trt = 3
rho = 0.5

pocockerror1bs = simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
pocockerror2bs = simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
pocockerror3bs = simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "Pocock")
spbderror1bs =   simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.0000), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
spbderror2bs =   simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
spbderror3bs =   simCarOrigBS(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary, design = "spbd")
nocarerror1bs =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3), 0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
nocarerror2bs =simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)
nocarerror3bs = simNoCarOrig(n1=n1, N1 = N1, N = N, n.trt = 3, mean.s = c(rep(0,3),-0.3333), mean.t = c(rep(0,3), 0.0), p1 = p1, p2 = p2, sigma0 = 1, sigma = 1, rho = rho, nsim = nsim, tau1 = 1, tau2 = 1, save.boundary = save.boundary)

pocockerror = cbind(pocockerror1,pocockerror2,pocockerror3)
spbderror = cbind(spbderror1,spbderror2,spbderror3)
nocarerror = cbind(nocarerror1,nocarerror2,nocarerror3)
bs = rbind(spbderror, nocarerror)
write.csv(bs, "C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/sbpdbs20_100_200.csv")
