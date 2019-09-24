#' get.z.v.knownsigma function
#'
#' gets z and v for all treatments at first look and for selected treatment at second look using known values for sigma and rho at first look.
#' Selection is of best treatment from first look
#' @param data data set with primary and secondary endpoint, and treatment
#' @param n1 Number of patients with primary endpoint available at first analysis
#' @param N1 Number of patients with secondary endpoint available at first analysis
#' @param N The total number of patients in the trial
#' @param sigma is the known sigma for the population
#' @param rho is the known correlation between endpoints.


get.z.v.knownsigma <- function(data,n1,N1,N,sigma,rho)
{
  # gets z and v for all treatments at first look and for selected treatment at second look
  # using known values for sigma and rho at first look
  # selection is of best treatment from first look

  n.trt = max(data$treat)

  z = rep(0,2*n.trt)
  dim(z) <- c(n.trt,2)
  v = rep(0,2*n.trt)
  dim(v) <- c(n.trt,2)


  # interim analysis
  # double regression as in Engel and Walstra (Biometrics 1991)

  s.N1 = data$s[data$treat==0][1:N1]
  for (treat in seq(1,n.trt)) s.N1 = c(s.N1,data$s[data$treat==treat][1:N1])
  treat.N1 = rep(seq(0,n.trt),rep(N1,n.trt+1))
  reg1 = lm(s.N1~factor(treat.N1))
  b.hat = reg1$coef[2:(n.trt+1)]

  s.n1 = data$s[data$treat==0][1:n1]
  for (treat in seq(1,n.trt)) s.n1 = c(s.n1,data$s[data$treat==treat][1:n1])
  t.n1 = data$t[data$treat==0][1:n1]
  for (treat in seq(1,n.trt)) t.n1 = c(t.n1,data$t[data$treat==treat][1:n1])
  treat.n1 = rep(seq(0,n.trt),rep(n1,n.trt+1))
  reg2 = lm(t.n1~factor(treat.n1)+s.n1)
  beta.hat = reg2$coef[2:(n.trt+1)]
  gamma.hat = reg2$coef[n.trt+2]

  B.hat = beta.hat + gamma.hat*b.hat

  expected.var.B.hat = 2*sigma^2*( (1-rho^2)/n1 + rho^2/N1)

  z[,1] = B.hat/expected.var.B.hat
  v[,1] = rep(1/expected.var.B.hat,n.trt)



  # final analysis - selecting best treatment

  # make treatment selection
  best = seq(1,n.trt)[z[,1]==max(z[,1])]
  if (length(best) > 1) best = sample(best,1) # breaks ties at random

  t.N = c(data$t[data$treat==0][1:N],data$t[data$treat==best][1:N])
  treat.N = c(rep(0,N),rep(best,N))
  reg = lm(t.N~factor(treat.N))

  B.hat = reg$coef[2]
  var.B.hat = vcov(reg)[2,2]

  z[best,2] = B.hat/var.B.hat
  v[best,2] = 1/var.B.hat



  data.frame(z1=z[,1],v1=v[,1],z2=z[,2],v2=v[,2])

}
