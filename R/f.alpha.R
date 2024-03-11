#' Calculation of Incremental Alpha for Group Sequential Design
#' 
#' The function includes the incremental alpha for group sequential design. 
#' The incremental alpha is a vector with sum to the overall alpha. It is may be needed for
#' calling the group sequential design power functions.
#' 
#' @param  overall.alpha  Allocated overall alpha for the GSD
#' @param  side Side of test
#' @param  sf Spending function. "LDOF", "LDPK"
#' @param  timing Timing of IA and FA
#' 
#' @export
#' 
f.alpha = function(overall.alpha=0.025, side = 1, sf="LDOF", timing=c(0.75, 1)){
  K = length(timing)
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  
  ld.obf = function(s){
    if (side == 1){a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))}
    if (side == 2){a = 2*2*(1 - pnorm(qnorm(1-overall.alpha/4)/sqrt(s)))}
    return(a)
  }
  ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
  
  if (sf == "LDOF"){
    gs.alpha = ld.obf(s = timing)
  }
  if (sf == "LDPK") {
    gs.alpha = ld.pk(s = timing)
  }
  if (K == 1){alpha = overall.alpha} else{
    alpha = rep(NA, K)
    
    alpha[1] = gs.alpha[1]
    for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
  }
  
  return(alpha)
}
