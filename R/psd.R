#' Function for returning treatment assignments in Pocock's Design
#'
#' This function takes a vector of factor level assignments, and returns treatment assignments
#' @param x Vector of factor level assignments
#' @param p Default is 3/4.
#' @param n.trt Number of treatments used in the trial
#' @param Best is the selected best treatment if past the first analysis. Default is 0
#' @param tr Tr is the vector of treamtnet assignments. Default is NULL
#' @param n.trt The number of treatments used in the analaysis if at the first stage in the analysis.
#' @export
#' @examples
#' psd(x = seamlessTrials::all)

psd=function(x,p1=3/4,best = 0,tr = NULL,n.trt)
{
  if (best == 0) trts = seq(0,n.trt) else trts = c(0,best)
  if (is.null(tr)) keeps = rep(TRUE,nrow(x)) else keeps = c(tr %in% trts,rep(TRUE,nrow(x)-length(tr)))
  if (is.null(tr)) arrival = as.matrix(ade4::acm.disjonctif(x)) else arrival = as.matrix(ade4::acm.disjonctif(x[keeps,]))
  # arrival = as.matrix(ade4::acm.disjonctif(x))
  weight = rep(1,ncol(arrival))
  s = length(tr[tr %in% trts]) + 1
  n.trt = length(trts)
  N = nrow(arrival)
  base = matrix(rep(0,n.trt*ncol(arrival)),ncol=ncol(arrival))
  for ( i in 1:length(trts)) base[i,] = colSums(arrival[which(tr[tr %in% trts]==trts[i]),])
  if (is.null(tr)) tr=rep(NA,N) else tr = tr[keeps]
  out = list()
  # trin = which(trts==trts)
  delts = rep(0,(n.trt))
  while(s<=N)
  {
    for ( i in 1:length(trts)) {

      inc = base
      inc[i,] = base[i,] + t(arrival[s,])
      delts[i]=t(arrival[s,])%*%(Rfast::colrange(inc)*weight)
    }
    p=g(delts,p1,n.trt)
    tr[s] = sample(trts,1,p,replace=TRUE)
    base[which(trts==tr[s]),] = base[which(trts==tr[s]),] + t(arrival[s,])
    s=s+1
  }
  tr
  }

