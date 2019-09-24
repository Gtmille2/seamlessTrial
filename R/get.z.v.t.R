#' get.z.v.t function
#'
#' gets z and v for all treatments at first look and for selected treatment at second look selection is prespecified or is of best treatment from first look if selected = 0
#' @param data data set with primary and secondary endpoint, and treatment
#' @param selected The selected treatment from the first look
#' @param n1 Number of patients with primary endpoint available at first analysis
#' @param N1 Number of patients with secondary endpoint available at first analysis
#' @param N The total number of patients in the trial


get.z.v.t <- function(data,selected,n1,N1,N)
{
  # gets z and v for all treatments at first look and for selected treatment at second look
  # selection is prespecified or is of best treatment from first look if selected = 0
  # uses final endpoint data only
  n.trt = max(data$treat)

  z = rep(0,2*n.trt)
  dim(z) <- c(n.trt,2)
  v = rep(0,2*n.trt)
  dim(v) <- c(n.trt,2)


  # interim analysis
  t.n1 = data$t[data$treat==0][1:n1]
  for (treat in seq(1,n.trt)) t.n1 = c(t.n1,data$t[data$treat==treat][1:n1])
  treat.n1 = rep(seq(0,n.trt),rep(n1,n.trt+1))


  reg = lm(t.n1~factor(treat.n1))
  B.hat = reg$coef[2:(n.trt+1)]
  var.B.hat = vcov(reg)[2:(n.trt+1),2:(n.trt+1)]

  z[,1] = B.hat/var.B.hat[1,1]
  v[,1] = rep(1/var.B.hat[1,1],n.trt)

  # final analysis - selecting best treatment

  # make treatment selection
  if (selected >0) best = selected
  else best = seq(1,n.trt)[z[,1]==max(z[,1])]
  if (length(best) > 1) best = sample(best,1) # breaks ties at random

  t.N = c(data$t[data$treat==0][1:N],data$t[data$treat==best][1:N])
  treat.N = c(rep(0,N),rep(best,N))
  reg = lm(t.N~factor(treat.N))
  B.hat = reg$coef[2]
  # var.B.hat = vcov(reg)[2,2]*(2*n-2)/(2*n) # this line replaced by unbiased estimate
  var.B.hat = vcov(reg)[2,2]

  z[best,2] = B.hat/var.B.hat
  v[best,2] = 1/var.B.hat


  data.frame(z1=z[,1],v1=v[,1],z2=z[,2],v2=v[,2])

}
