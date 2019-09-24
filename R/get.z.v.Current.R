#' Get Z and V function for any look
#'
#' This function returns the Z and V for the current look.
#' @param data Data containing the treatment assignments, short-term outcomes, and long-term outcomes
#' @param n1 The number of patients with data available for the long-term endpoint at the current analysis
#' @param N2 The number of patients with data available for the short-term endpoint at the current analysis
#' @param best Best is the selected treatment if beyond the first analysis. If at first analysis, then best equals 0.
#' @export

get.z.v.Current <- function(data,n.looks,look,z.v.prev = NULL)
{
  # gets z and v for all treatments at first look and for selected treatment at second look
  # selection is of best treatment from first look

  n.trt = length(unique(data$treat))-1
  trts = unique(data$treat)
  # keep = unique(c(0,trts[trts!=best]))

  z = rep(0,length(trts)-1)
  # dim(z) <- c(n.trt,2)
  v = rep(0,length(trts)-1)
  # dim(v) <- c(n.trt,2)


  # interim analysis
  # double regression as in Engel and Walstra (Biometrics 1991)
  # N1 is the short term endpoint
  # s.N1 = data$s[data$treat==0][1:N1]
  # for (treat in seq(1,n.trt)) s.N1 = c(s.N1,data$s[data$treat==treat][1:N1])
  s.N1=NULL
  for (treat in trts) s.N1 = c(s.N1,data$s[data$treat==treat])

  treat.N1 = rep(trts,table(data$treat)) # Does this accuratly account for Pocock?
  reg1 = lm(s.N1~factor(treat.N1))
  b.hat = reg1$coef[2:(n.trt+1)]
  var.b.hat = vcov(reg1)[2:(n.trt+1),2:(n.trt+1)]

  sigma.0.hat = summary(reg1)$sigma


  # s.n1 = data$s[data$treat==0][1:n1]
  # for (treat in seq(1,n.trt)) s.n1 = c(s.n1,data$s[data$treat==treat][1:n1])
  #
  s.n1 = na.omit(data)$s
  # t.n1 = data$t[data$treat==0][1:n1]
  t.n1 = na.omit(data)$t
  # for (treat in seq(1,n.trt)) t.n1 = c(t.n1,data$t[data$treat==treat][1:n1])
  treat.n1 = rep(trts,table(na.omit(data)$treat))
  reg2 = lm(t.n1~factor(treat.n1)+s.n1)
  beta.hat = reg2$coef[2:(n.trt+1)]
  var.beta.hat = vcov(reg2)[2:(n.trt+1),2:(n.trt+1)]
  gamma.hat = reg2$coef[n.trt+2]
  var.gamma.hat = vcov(reg2)[n.trt+2,n.trt+2]
  cov.beta.hat.gamma.hat = vcov(reg2)[2:(n.trt+1),n.trt+2]

  sigma.1.hat = summary(reg2)$sigma

  B.hat = beta.hat + gamma.hat*b.hat

  sigma.hat = sqrt( sigma.1.hat^2 + sigma.0.hat^2 * gamma.hat^2)
  rho.hat = gamma.hat*sigma.0.hat/sigma.hat


  expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/table(na.omit(data$treat))[1:n.trt] + rho.hat^2/table(data$treat)[1:n.trt])
  z.v = rep(0,n.looks*n.trt*2)
  dim(z.v) <- c(n.trt,n.looks,2)
  z = B.hat/expected.var.B.hat
  v = (1/expected.var.B.hat)
  if (!is.null(z.v.prev)) z.v[,,1] = z.v.prev[,1]
  if (!is.null(z.v.prev)) z.v[,,2] = z.v.prev[,2]
  z.v[,look,1] = z
  z.v[,look,2] = v

  # final analysis - selecting best treatment

  # make treatment selection

  best = trts[trts!=0][z==max(z)]
  if (length(best) > 1) best = sample(best,1) # breaks ties at random



  data.frame(z = z.v[which(trts[trts!=0]==best),,1],v = z.v[which(trts[trts!=0]==best),,2],treat=best)

}
