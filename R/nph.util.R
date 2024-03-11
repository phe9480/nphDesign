
#' @export
#' 
lr = function(s){1}
fh01 = function(s){(1-s)}
fh11 = function(s){s*(1-s)}
fh55 = function(s){sqrt(s*(1-s))}

#' #stabilized FH tests
sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 
sfh11 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 
sfh55 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(sqrt(s1*(1-s1)))} 

#' #modestly log-rank
mlr = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1/s1)}


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

