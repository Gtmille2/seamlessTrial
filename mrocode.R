
n1 = 20
N1 =100
N = 200
p1 = .4
p2 = .4

power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower20_100_200_5000_rho44.csv")


n1 = 20
N1 =100
N = 200
p1 = .5
p2 = .4
nsim = 5000
power_rho_0 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_5 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.5, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
power_rho_7 = allsimpower(n1 = n1, N1 = N1, N = N, n.trt = 3, p1 = p1, p2 = p2,sigma0 = 1, sigma= 1, rho = 0.7, nsim = nsim, tau1  = 1, tau2 = 1, save.boundary = save.boundary)
allsimpower1 = rbind(power_rho_0,power_rho_5,power_rho_7)
write.csv(allsimpower1,"C:/Users/garre/OneDrive/Documents/UTHealth Files/Thesis/allsimpower20_100_200_5000_rho54.csv")
