## FZ: 
## Replace the accr input list by F.entry and G.ltfu as input parameters
## Implement gamma alpha spending function
## Call the updated function wlr.power.maxcombo.update

finalize.nphDesign.update = function(T, n = 300, r = 1, dist = dist, 
                                     F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                                     G.ltfu = function(t){0}, 
                                     alphabeta = alphabeta, f.ws = list(IA1 = list(lr), FA = list(lr, fh01)),
                                     show.setting="Y", non.centrality = "Heetal2021"){
  
  #Number of analyses
  K = length(f.ws)
  J = lengths(f.ws) 
  
  #Construct distribution functions h0, S0, h1, S1
  if(dist$control$dist=="exponential"){
    lam0 = log(2) / dist$control$median
    if (dist$crossover$status=="N"){
      h0 = function(t){lam0}
      S0 =  function(t){1-pexp(t, rate=lam0)}
    } else {
      h0 = function(t){
        HR = dist$exp$HR; HRx = dist$crossover$HRx
        lt = as.numeric(t < dist$crossover$Tx)
        return(lam0 * lt + HR / HRx * lam0 * (1 - lt))
      } 
      S0 = function(t){
        HR = dist$exp$HR; HRx = dist$crossover$HRx; Tx = dist$crossover$Tx
        lt = as.numeric(t < Tx)
        c0 = exp(-Tx*lam0*(1 - HR / HRx))
        return(exp(-lam0*t)*lt + c0*exp(-HR/HRx * lam0 * t) * (1 - lt))
      }
    }  
  }
  
  h1 = function(t){
    delay = dist$exp$delay; HR = dist$exp$HR; lt = as.numeric(t < delay)
    lam0 * lt + HR * lam0 * (1 - lt)
  }
  S1 = function(t){
    delay = dist$exp$delay; HR = dist$exp$HR; lt = as.numeric(t < delay)
    c1 = exp(-delay*lam0*(1-HR));   
    exp(-lam0 * t) * lt + c1 * exp(-HR * lam0 * t) * (1 - lt)
  }
  f.logHR = function(t){
    delay = dist$exp$delay; HR = dist$exp$HR
    if (dist$crossover$status == "N"){
      return(log(HR)*as.numeric(t >= delay))
    } else{
      HRx = dist$crossover$HRx; Tx = dist$crossover$Tx
      btw = as.numeric(t >= delay & t < Tx)
      ge = as.numeric(t >= Tx)
      return(log(HR) * btw + ge * log(HRx))
    }
  }
  
  #Construct accrual functions
  
  #A = accr$A; xi = accr$xi; 
  #G.ltfu = accr$G.ltfu
  #F.entry = function(t){(t/A)^xi*as.numeric(t <= A) + as.numeric(t > A)}
  
  alpha = alphabeta$alpha; side = alphabeta$side
  overall.alpha = alphabeta$overall.alpha; 
  
  
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(alpha) && !is.null(overall.alpha)){
    ld.obf = function(s){
      if (side == 1){a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))}
      if (side == 2){a = 2*2*(1 - pnorm(qnorm(1-overall.alpha/4)/sqrt(s)))}
      return(a)
    }
    ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
    # add gamma sf by FZ
    ld.gamma = function(s, param = -2){
      overall.alpha * (1 - exp(-s * param))/(1 - exp(-param))
    }
    #allocate alpha based on proportional of events, pseudo information
    nE = rep(NA, K)
    for (j in 1:K){
      nE[j] = f.nEvents(T = T[j], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                        F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events$n.events.total
    }
    timing = nE/nE[K]
    
    if (alphabeta$sf$type == "LDOF"){
      gs.alpha = ld.obf(s = timing)
    }
    if (alphabeta$sf$type == "LDPK") {
      gs.alpha = ld.pk(s = timing)
    }
    ## add gamma sf by FZ
    if (alphabeta$sf$type == "GAMMA") {
      gs.alpha = ld.gamma(s = timing, param=alphabeta$sf$param)
    }
    if (K == 1){alpha = overall.alpha} else{
      alpha[1] = gs.alpha[1]
      for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
    }
  }
  
  #Call Group sequential design power function
  gs.power = wlr.power.maxcombo.update(T = T, events = NULL, 
                                       alpha = alpha, power = NULL, side = side, r = r, n = n, 
                                       h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
                                       f.ws = f.ws, F.entry=F.entry, G.ltfu=G.ltfu,
                                       show.setting=show.setting, non.centrality = non.centrality)
  
  setting = list()
  setting$dist = dist; #setting$accr = accr; 
  setting$alphabeta = alphabeta;
  o = list()
  o$design = gs.power
  #o$setting = setting
  
  return(o)   
}