#' Graphic Exploration of Study Design Setting
#' 
#' This function visualizes the following information: (1) Survival curves;
#' (2) Cumulative event curves; (3) Expected average hazard ratio curve; 
#' (4) Fisher's information curves; (5) Non-centrality parameter curves; 
#' (6) Projected z boundary curves; (7) Power curves; 
#' (8) Projected p value boundary curve for non-combination tests; 
#' 
#' @param maxT  Maximum calendar time in design considerations. Default 50 months.
#' @param dist  An object that includes distributions of two arms in the design assumptions
#' @param accr  A list including accrual distribution specifications. Either power law or user-defined accrual function is supported.
#'              If using power law, the accrual duration A and weight parameter xi need to be specified, which corresponds 
#'              to the accrual function F.entry = function(t){(t/A)^xi*as.numeric(t <= A) + as.numeric(t > A)}
#'              If using self-defined accrual function directly for complex accrual patterns, then F.entry needs to be provided.
#'              G.ltfu is an optional user-customized dropout function. The default is for no dropout with G.ltfu=function(t){0}.
#' @param size  Target sample size
#' @param alphabeta   Type I and type II erros
#' @param AHR.method  Method of calculating average hazard ratio: "geometric.schoenfeld" or "Kalbfleisch and Prentice"
#' @param f.ws  Weight functions for statistical tests.
#' @param show.setting Option whether to show the recommended study design settings
#' @param s.leg Survival plot legend
#' @param s.lab Survival plot labels
#' @param e.leg Events plot legend
#' @param e.lab Events plot labels
#' 
#' @return study design
#'  
#' @examples 
#' 
#' 
#' her2p.dist = list(
#' control = list(dist = "exponential", median = 6.7),
#' exp = list(HR = 6.7/12.5, delay = 6),
#' crossover = list(status = "N", HRx = NULL, Tx = NULL)
#' ) 
#' her2p.accr = list(
#'   A = 23,
#'   xi = 2,
#'   LTFU = 0
#' )
#' her2p.size = list(
#'   n = 300,
#'   r = 1,
#'   maturity = 0.7
#' )
#' 
#' #recommend side = 1 for superiority design
#' #alpha takes priority than sf if alpha is provided.
#' #sf options: "LDOF", "LDPK"
#' her2p.alphabeta = list(
#'   sf = list(type = c("LDOF")),
#'   alpha = c(0.02, 0.03)/2,
#'   overall.alpha = 0.025,
#'   beta = 0.1,
#'   side = 1,
#'   timing = c(0.75, 1.0)
#' )
#' 
#' lr = function(s){1}
#' fh01 = function(s){1-s}
#' fh11 = function(s){s*(1-s)}
#' 
#' her2p.test = list(IA1 = list(lr), FA=list(lr, fh01))
#'  
#' explore.nphDesign(dist=her2p.dist, accr=her2p.accr, size=her2p.size, 
#'      alphabeta=her2p.alphabeta, 
#'      AHR.method="Kalbfleisch and Prentice", f.ws = her2p.test)
#'  
#' @export
explore.nphDesign = function(maxT=50, dist=dist, accr=accr, size=size, 
        alphabeta=alphabeta, AHR.method="geometric.schoenfeld",
        f.ws = list(IA1 = list(lr), FA=list(lr, fh01)), show.setting="Y",
        non.centrality = "Heetal2021",
        s.leg = list(x=30, y=1, txt=c("Control", "Exp. Arm")), 
        s.lab = list(xlab="Survival Time (mo)", ylab="Survival",
                       main="Survival Curve Per Study Design"),
        e.leg = list(txt = c("Control", "Exp. Arm", "Total")), 
        e.lab = list(xlab = "Calendar Time (mo)", ylab="Cumulative Events", main = "")
        ){
  
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
  
  # Construct accrual functions
  # FZ 05/07/2024
  # if the accrual duration A and accrual power xi are provided (xi=1 if not specified), then F.entry will be defined
  # if F.entry is provided instead A and xi, then use F.entry directly
  # if neither A (with xi) nor F.entry is defined, then warning message will be returned.
  if(!is.null(accr$A) & is.null(accr$F.entry)){
    if(!is.null(accr$xi)){
      A = accr$A; xi = accr$xi; 
    }else{
      warning("The accrual power parameter xi is not defined. The default value xi=1 will be used.")
      xi = 1
    }
    F.entry = function(t){(t/A)^xi*as.numeric(t <= A) + as.numeric(t > A)}
  }else if(is.null(accr$A) & !is.null(accr$F.entry)){
    F.entry = accr$F.entry
  }else if(is.null(accr$A) & is.null(accr$F.entry)){
    warning("A accrual list should be defined.")
  }else if(!is.null(accr$A) & !is.null(accr$F.entry)){
    warning("Both power law accrual format (A, xi) and customized accrual function (F.entry) are provided. 
             F.entry will be prioritized for use.")
    F.entry = accr$F.entry
  }
  
  if(is.null(accr$G.ltfu)){
    G.ltfu = function(t){0}
  }else{
    G.ltfu = accr$G.ltfu
  }
  
  t = seq(0, A, by = 1)
  cum.enroll = F.entry(t)
  plot(t, cum.enroll, type="l", lwd=3, xlab = "Calendar Time", ylab="Cumulative Accrual")
  abline(v=seq(0, A, 1), col="gray80", lty=3)
  abline(h=seq(0, 1, 0.1), col="gray80", lty=3)
  
  #(1)Display distributions
  plot_S(S=list(S0, S1), Tmax = maxT,leg=s.leg, param=s.lab)
  
  #(2)Display cumulative events over time
  r = size$r; n = size$n
  plot_events(Tmax = maxT, r=r, h0=h0, S0=S0, h1=h1, S1=S1, 
    F.entry = F.entry, G.ltfu = G.ltfu, n = n, leg = e.leg, param = e.lab)
  
  #(3) Display Expected Average hazard ratio over calendar time
  plot_AHR(r = r, n =n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
      F.entry = F.entry, G.ltfu = G.ltfu, Tmax=maxT, method=AHR.method)
  
  #(4) Fisher's information curve for all tests used in FA
  plot_wlr.info(Tmax = maxT, r = r, h0 = h0, S0= S0, h1 = h1, S1=S1, 
     fws = f.ws[[K]], F.entry = F.entry, G.ltfu = G.ltfu, n=n)
  
  #(5) Non-centrality parameter curves
  plot_wlr.mu (Tmax = maxT, r = r, n = n,h0 = h0, S0= S0, h1 = h1, S1=S1, 
   f.logHR = f.logHR, fws = f.ws[[K]], F.entry = F.entry, G.ltfu = G.ltfu)
 
  #(6) Projected z boundary curves
  ##(6.1) Determine alpha for each analysis
  alpha = alphabeta$alpha; beta = alphabeta$beta; side = alphabeta$side
  timing = alphabeta$timing; overall.alpha = alphabeta$overall.alpha
  m = size$maturity
  
  n0 = size$n; e0 = n0 * m
  
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(alpha) && !is.null(overall.alpha)){
    ld.obf = function(s){
      if (side == 1){a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))}
      if (side == 2){a = 2*2*(1 - pnorm(qnorm(1-overall.alpha/4)/sqrt(s)))}
      return(a)
    }
    ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
    ld.gamma<-function(s,param=-2){overall.alpha * (1-exp(-s*param))/(1-exp(-param))}
    
    if (alphabeta$sf$type == "LDOF"){
      gs.alpha = ld.obf(s = timing)
    }
    if (alphabeta$sf$type == "LDPK") {
      gs.alpha = ld.pk(s = timing)
    }
    if (alphabeta$sf$type == "GAMMA") {
      gs.alpha = ld.gamma(s = timing,param=alphabeta$sf$param)
    }
    if (K == 1){alpha = overall.alpha} else{
      alpha[1] = gs.alpha[1]
      for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
    }
    
  }
  
  ##(6.2) Initial calculation of sample size n based on FA
  init.power = wlr.power.maxcombo(T = NULL, events = e0, 
     alpha = overall.alpha, power = NULL, side = side, r = r, n = n0, 
     h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
     f.ws = f.ws[K], F.entry=F.entry, G.ltfu=G.ltfu, 
     non.centrality = non.centrality)

  pow0 = init.power$overall.power
  while (pow0 < 1 - beta){
    e0 = e0 + 2;  n0 = round(e0 / m);
    while(n0 %% (r+1) > 0){n0 = n0 + 1; e0 = n0 * m}
    pow0 = wlr.power.maxcombo(T = NULL, events = e0, 
               alpha = overall.alpha, power = NULL, side = side, r = r, n = n0, 
               h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
               f.ws = f.ws[K], F.entry=F.entry, G.ltfu=G.ltfu, non.centrality = non.centrality)$overall.power
  }

  #(6.3) Initial group sequential design
   init.gsPow = wlr.power.maxcombo(T = NULL, events = e0 * timing, 
          alpha = alpha, power = NULL, side = side, r = r, n = n0, 
          h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
          f.ws = f.ws, F.entry=F.entry, G.ltfu=G.ltfu, show.setting=show.setting, 
          non.centrality = non.centrality)
   
   #(6.4) Finalizing group sequential design by recursive solving power
   gsPow0 = init.gsPow$overall.power
   gs.power = init.gsPow

   #If there is overshoot over 1% power, then decrease the sample size
   while (gsPow0 > 1 - beta + 0.005){
     e0 = e0 - 2;  n0 = round(e0 / m);
     while(n0 %% (r+1) > 0){n0 = n0 - 1; e0 = n0 * m}
     gs.power = wlr.power.maxcombo(T = NULL, events = e0 * timing, 
                      alpha = alpha, power = NULL, side = side, r = r, n = n0, 
                      h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
                      f.ws = f.ws, F.entry=F.entry, G.ltfu=G.ltfu, 
                      show.setting=show.setting, non.centrality = non.centrality)
     gsPow0 = gs.power$overall.power
   }
   
   #If the power is not enough, increase N
   while (gsPow0 < 1 - beta){
     e0 = e0 + 1;  n0 = round(e0 / m);
     while(n0 %% (r+1) > 0){n0 = n0 + 1; e0 = n0 * m}
     gs.power = wlr.power.maxcombo(T = NULL, events = e0 * timing, 
               alpha = alpha, power = NULL, side = side, r = r, n = n0, 
               h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
               f.ws = f.ws, F.entry=F.entry, G.ltfu=G.ltfu, 
               show.setting=show.setting, non.centrality = non.centrality)
     gsPow0 = gs.power$overall.power
   }
   
   #(6.4) Display the z-boundary
   plot(gs.power$events$events, gs.power$bounds[,1], type="b", lwd=2, 
        xlab = "Analysis Timing Based on Events", ylab = "z Rejection boundary", 
        main = "Projected Rejection Boundary (Z)")
   abline(v = seq(0, gs.power$events$events[K]*2, 10), col="gray80", lty=3)
   abline(h = seq(0, 10, 0.05), col="gray80", lty=3)
   
   #(7) Power curve vs events 
   e.seq = seq(max(10, e0/50), min(e0*1.2, n0*0.85), length.out = 30)
   overall.pow.seq = rep(NA, length(e.seq))
   pow.seq = matrix(NA, nrow=length(e.seq), ncol=K)
   
   for (i in 1:length(e.seq)){
     wlr.seq = wlr.power.maxcombo(T = NULL, events = e.seq[i] * timing, 
            alpha = alpha, power = NULL, side = side, r = r, n = n0, 
            h0 = h0, S0=S0, h1 = h1, S1= S1, f.logHR = f.logHR, 
            f.ws = f.ws, F.entry=F.entry, G.ltfu=G.ltfu, non.centrality = non.centrality)
     overall.pow.seq[i] = wlr.seq$overall.power
     pow.seq[i,] = wlr.seq$power
   }
   #Overall power
   plot(e.seq, overall.pow.seq, type="l", lwd=3, xlab="Number of Events at Final Analysis", 
        ylab="Power", main="Overall Power")
   
   abline(h = seq(0, 1, 0.05), col="gray80", lty=3)
   abline(v = seq(0, max(e.seq), by=10), col="gray80", lty=3)
   
   abline(v=e0*timing, col="green", lwd=3)
   
   #Power of each analysis
   for(j in 1:K){
     plot(e.seq*timing[j], pow.seq[,j], type="l", lwd=3, xlab = paste("Number of Events at Analysis", j), 
          ylab="Power", main=paste("Power of analysis", j))
     
     abline(h = seq(0, 1, 0.05), col="gray80", lty=3)
     abline(v = seq(0, max(e.seq*timing[j]), by=10), col="gray80", lty=3)
     
     abline(v=e0*timing[j], col="green", lwd=3)
   }
   
   #(8) P-value boundary for non-combination test at each analysis
   if (sum(J) == K){
     if(side == 1){
       p.bd = 1 - pnorm(gs.power$bounds[,1])
     }
     if (side == 2){
       p.bd = 2*(1 - pnorm(gs.power$bounds[,1]))
     }
     plot(gs.power$events$events, p.bd, type="b", lwd=2, 
          xlab = "Analysis Timing Based on Events", ylab="P value", 
          main = "Projected Rejection Boundary (p-value)")
     abline(v = seq(0, gs.power$events$events[K]*2, 10), col="gray80", lty=3)
     abline(h = seq(0, 1, 0.001), col="gray80", lty=3)
     
   }
  
  setting = list()
  setting$dist = dist; setting$accr = accr; setting$size = size; 
  setting$alphabeta = alphabeta;
  o = list()
  o$design = gs.power
  #o$setting = setting
  
  return(o)   
}

