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
   if ( look == n.looks) z.v = get.z.v.final(data, n.looks, look, z.v.prev) else z.v = get.z.v.first(data,n.looks,look,z.v.prev = NULL)


}

get.z.v.final = function(data, n.looks,look,z.v.prev) {

  n.trt = length(unique(data$treat))-1
  z = rep(0,2*n.trt)
  dim(z) <- c(n.trt,2)
  v = rep(0,2*n.trt)
  dim(v) <- c(n.trt,2)
  best = z.v.prev$treat[1]
  t.N = data$t
  treat.N = data$treat
  reg = lm(t.N~factor(treat.N))
  B.hat = reg$coef[2]
  var.B.hat = vcov(reg)[2,2]
  z[1,2] = B.hat/var.B.hat
  v[1,2] = 1/var.B.hat
  z.v = z.v.prev
  z.v$z[2] = z[1,2]
  z.v$v[2] = v[1,2]
  # t.N = c(data[data$treat == 0,]$t, data[data$treat == best,]$t)
  z.v
}

get.z.v.first = function(data, n.looks, look, z.v.prev=NULL) {
  n.trt = length(unique(data$treat))-1
  trts = unique(data$treat)
  trts = trts[order(trts)]
  z = rep(0,length(trts)-1)
  # dim(z) <- c(n.trt,2)
  v = rep(0,length(trts)-1)

  s.N1=NULL
  # for (treat in trts) s.N1 = c(s.N1,data$s[data$treat==treat])
  treat.N1 = data$treat
  s.N1 = data$s
  # treat.N1 = rep(trts,table(data$treat)) # Does this accuratly account for Pocock?
  reg1 = lm(s.N1~factor(treat.N1))
  b.hat = reg1$coef[2:(n.trt+1)]
  var.b.hat = vcov(reg1)[2:(n.trt+1),2:(n.trt+1)]

  sigma.0.hat = summary(reg1)$sigma


  # s.n1 = data$s[data$treat==0][1:n1]
  # for (treat in seq(1,n.trt)) s.n1 = c(s.n1,data$s[data$treat==treat][1:n1])
  #
  # s.n1 = na.omit(data)$s
  s.n1=NULL
  s.n1 = na.omit(data)$s
  # for (treat in trts) s.n1 = c(s.n1,na.omit(data)$s[na.omit(data)$treat==treat])

  # t.n1 = data$t[data$treat==0][1:n1]
  # t.n1 = na.omit(data)$t
  t.n1=NULL
  t.n1 = na.omit(data)$t

  # for (treat in trts) t.n1 = c(t.n1,na.omit(data)$t[na.omit(data)$treat==treat])
  # treat.n1 = rep(trts,table(na.omit(data)$treat))
  treat.n1 = na.omit(data)$treat

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


  expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/table(na.omit(data)$treat)[2:(n.trt+1)] + rho.hat^2/table(data$treat)[2:(n.trt+1)])
  z.v = rep(0,n.looks*n.trt*2)
  dim(z.v) <- c(n.trt,n.looks,2)
  z = B.hat/expected.var.B.hat
  v = (1/expected.var.B.hat)
  z.v[,look,1] = z
  z.v[,look,2] = v

  # final analysis - selecting best treatment

  # make treatment selection

  best = trts[trts!=0][z==max(z)]
  if (length(best) > 1) best = sample(best,1) # breaks ties at random



  data.frame(z = z.v[which(trts[trts!=0]==best),,1],v = z.v[which(trts[trts!=0]==best),,2],treat=best)
}


get.z.v.simulate = function(data, n.looks, look, z.v.prev, n1, N1, N) {

  Ns = table(data$treat)
  n1s = round(table(data$treat)*(n1/N))
  N1s = round(table(data$treat)*(N1/N))
  n.trt = length(unique(data$treat))-1
  trts = unique(data$treat)
  trts = trts[order(trts)]
  z = rep(0,length(trts)-1)
  # dim(z) <- c(n.trt,2)
  v = rep(0,length(trts)-1)
  look = 1
  s.N1=NULL
  for (treat in trts) s.N1 = c(s.N1,data$s[data$treat==treat][1:N1s[which(names(N1s)==treat)]])
  treat.N1 = NULL
  for (treat in trts) treat.N1 = c(treat.N1, data$treat[data$treat==treat][1:N1s[which(names(N1s)==treat)]])
  reg1 = lm(s.N1~factor(treat.N1))
  b.hat = reg1$coef[2:(n.trt+1)]
  var.b.hat = vcov(reg1)[2:(n.trt+1),2:(n.trt+1)]
  sigma.0.hat = summary(reg1)$sigma
  s.n1=NULL
  # s.n1 = na.omit(data)$s
  for (treat in trts) s.n1 = c(s.n1,na.omit(data)$s[na.omit(data)$treat==treat][1:n1s[which(names(n1s)==treat)]])

  # t.n1 = data$t[data$treat==0][1:n1]
  # t.n1 = na.omit(data)$t
  t.n1=NULL
  # t.n1 = na.omit(data)$t

  for (treat in trts) t.n1 = c(t.n1,na.omit(data)$t[na.omit(data)$treat==treat][1:n1s[which(names(n1s)==treat)]])
  # treat.n1 = rep(trts,table(na.omit(data)$treat))
  # treat.n1 = na.omit(data)$treat
  treat.n1 = NULL
  for (treat in trts) treat.n1 = c(treat.n1, data$treat[data$treat==treat][1:n1s[which(names(n1s)==treat)]])
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

  expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/n1s[2:(n.trt+1)] + rho.hat^2/N1s[2:(n.trt+1)])
  # expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/n1 + rho.hat^2/N1)

  # expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/table(na.omit(data)$treat)[2:(n.trt+1)] + rho.hat^2/table(data$treat)[2:(n.trt+1)])
  z.v = rep(0,n.looks*n.trt*2)
  dim(z.v) <- c(n.trt,n.looks,2)
  z = B.hat/expected.var.B.hat
  v = (1/expected.var.B.hat)
  z.v[,look,1] = z
  z.v[,look,2] = v

  # final analysis - selecting best treatment

  # make treatment selection

  best = trts[trts!=0][z==max(z)]
  keeps = c(0,best)
  if (length(best) > 1) best = sample(best,1) # breaks ties at random
  # best = z.v.prev$treat[1]
  # t.N = data[data$treat %in% c(0,best),]$t
  # treat.N = data[data$treat%in% c(0,best),]$treat
  # reg = lm(t.N~factor(treat.N))
  t.N = NULL
  for (treat in keeps) t.N = c(t.N, data$t[data$treat==treat][1:Ns[which(names(Ns)==treat)]])
  treat.N = NULL
  for (treat in keeps) treat.N = c(treat.N, data$treat[data$treat==treat][1:Ns[which(names(Ns)==treat)]])
  reg = lm(t.N~factor(treat.N))

  B.hat = reg$coef[2]
  var.B.hat = vcov(reg)[2,2]
  z.v[,,1][best,2] = B.hat/var.B.hat
  z.v[,,2][best,2] = 1/var.B.hat
  # z.v = z.v.prev
  # z.v$z[2] = z[1,2]
  # z.v$v[2] = v[1,2]
  # t.N = c(data[data$treat == 0,]$t, data[data$treat == best,]$t)
  z.v


  data.frame(z = z.v[which(trts[trts!=0]==best),,1],v = z.v[which(trts[trts!=0]==best),,2],treat=best)

}
get.z.v.bootstrap = function(data, n.looks, look, z.v.prev, n1, N1, N) {

  Ns = table(data$treat)
  n1s = round(table(data$treat)*(n1/N))
  N1s = round(table(data$treat)*(N1/N))
  n.trt = length(unique(data$treat))-1
  trts = unique(data$treat)
  trts = trts[order(trts)]
  z = rep(0,length(trts)-1)
  # dim(z) <- c(n.trt,2)
  v = rep(0,length(trts)-1)
  blist = list()
  look = 1
  s.N1=NULL
  for (treat in trts) s.N1 = c(s.N1,data$s[data$treat==treat][1:N1s[which(names(N1s)==treat)]])
  treat.N1 = NULL
  for (treat in trts) treat.N1 = c(treat.N1, data$treat[data$treat==treat][1:N1s[which(names(N1s)==treat)]])
  reg1 = lm(s.N1~factor(treat.N1))
  b.hat = reg1$coef[2:(n.trt+1)]
  var.b.hat = vcov(reg1)[2:(n.trt+1),2:(n.trt+1)]
  sigma.0.hat = summary(reg1)$sigma
  s.n1=NULL
  # s.n1 = na.omit(data)$s
  for (treat in trts) s.n1 = c(s.n1,na.omit(data)$s[na.omit(data)$treat==treat][1:n1s[which(names(n1s)==treat)]])

  # t.n1 = data$t[data$treat==0][1:n1]
  # t.n1 = na.omit(data)$t
  t.n1=NULL
  # t.n1 = na.omit(data)$t

  for (treat in trts) t.n1 = c(t.n1,na.omit(data)$t[na.omit(data)$treat==treat][1:n1s[which(names(n1s)==treat)]])
  # treat.n1 = rep(trts,table(na.omit(data)$treat))
  # treat.n1 = na.omit(data)$treat
  treat.n1 = NULL
  for (treat in trts) treat.n1 = c(treat.n1, data$treat[data$treat==treat][1:n1s[which(names(n1s)==treat)]])
  reg2 = lm(t.n1~factor(treat.n1)+s.n1)
  beta.hat = reg2$coef[2:(n.trt+1)]
  var.beta.hat = vcov(reg2)[2:(n.trt+1),2:(n.trt+1)]
  gamma.hat = reg2$coef[n.trt+2]
  var.gamma.hat = vcov(reg2)[n.trt+2,n.trt+2]
  cov.beta.hat.gamma.hat = vcov(reg2)[2:(n.trt+1),n.trt+2]

  sigma.1.hat = summary(reg2)$sigma

  B.hat = beta.hat + gamma.hat*b.hat
  blist[[1]] = B.hat
  sigma.hat = sqrt( sigma.1.hat^2 + sigma.0.hat^2 * gamma.hat^2)
  rho.hat = gamma.hat*sigma.0.hat/sigma.hat

  expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/n1s[2:(n.trt+1)] + rho.hat^2/N1s[2:(n.trt+1)])
  # expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/n1 + rho.hat^2/N1)

  # expected.var.B.hat = 2*sigma.hat^2*( (1-rho.hat^2)/table(na.omit(data)$treat)[2:(n.trt+1)] + rho.hat^2/table(data$treat)[2:(n.trt+1)])
  z.v = rep(0,n.looks*n.trt*3)
  dim(z.v) <- c(n.trt,n.looks,3)
  z = B.hat/expected.var.B.hat
  v = (1/expected.var.B.hat)
  z.v[,look,1] = z
  z.v[,look,2] = v
  z.v[,look,3] = B.hat
  # final analysis - selecting best treatment

  # make treatment selection

  best = trts[trts!=0][z==max(z)]
  best <<- best
  keeps = c(0,best)
  if (length(best) > 1) best = sample(best,1) # breaks ties at random
  # best = z.v.prev$treat[1]
  # t.N = data[data$treat %in% c(0,best),]$t
  # treat.N = data[data$treat%in% c(0,best),]$treat
  # reg = lm(t.N~factor(treat.N))
  t.N = NULL
  for (treat in keeps) t.N = c(t.N, data$t[data$treat==treat][1:Ns[which(names(Ns)==treat)]])
  treat.N = NULL
  for (treat in keeps) treat.N = c(treat.N, data$treat[data$treat==treat][1:Ns[which(names(Ns)==treat)]])
  reg = lm(t.N~factor(treat.N))

  B.hat = reg$coef[2]
  var.B.hat = vcov(reg)[2,2]
  z.v[,,1][best,2] = B.hat/var.B.hat
  z.v[,,2][best,2] = 1/var.B.hat
  z.v[,,3][best,2] = B.hat
  # z.v = z.v.prev
  # z.v$z[2] = z[1,2]
  # z.v$v[2] = v[1,2]
  # t.N = c(data[data$treat == 0,]$t, data[data$treat == best,]$t)
  z.v


  # data.frame(z = z.v[which(trts[trts!=0]==best),,1],v = z.v[which(trts[trts!=0]==best),,2],treat=best)
  #
}
