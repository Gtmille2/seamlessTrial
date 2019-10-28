##### This script is to test different parameter settings for the third paper to determine how they affect the outcome of the type I error #####


library(seamlessTrials)
#Running multiple simulations:
save.boundary = save.boundary.values()
n.trt = 3
nsim = 100
N1 = 250
N=500
sigma0 = 1
sigma = 1
p1 = .5
p2 = .5
rho = 0.5
tau1 = 1
tau2 = 1
mean.s = NULL
mean.t = NULL

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
nocarreject1 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = 0,nsim = 100,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarreject2 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0.333),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarreject3 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),-0.333),mean.t = rep(0,n.trt+1), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower1 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 1000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower2 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),0.333),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
ptm = proc.time()
nocarpower3 = simulateNoCar(N1 = 250,N = 500,n.trt = 3,mean.s  = c(rep(0,3),-0.333),mean.t = c(rep(0,3),0.333), p1 = .5, p2 =.5,sigma0 = 1,sigma = 1,rho = .5,nsim = 10000,tau1 = 1,tau2 = 1,save.boundary = save.boundary)
proc.time() - ptm
