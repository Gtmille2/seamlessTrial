#' PBR Function
#'
#' This function is not generally used. This is used in spbd to determine treatment assignments
#' @param n Total n size
#' @param block.size This is the number of subjects in each block




pbr=function(n,block.size=4) #block.size is the number of subjects in each block, assuming fixed
{
  block.num=ceiling(n/block.size)
  cards=NULL
  i=1
  while(i<=block.num)
  {
    cards=c(cards,sample(cbind(rep(1,block.size/2),rep(0,block.size/2)),block.size))
    i=i+1
  }
  cards=cards[1:n]
  return(cards)
}
