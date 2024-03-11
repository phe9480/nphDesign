#' Non-centrality Parameter Under H1 For A Weighted Log-rank Test
#' 
#'  This function calculates the non-centrality parameter under H1 using (He et al 2021)
#'  method and evaluated at a calendar time with observable data. Both Schoenfeld (1981) 
#'  and this method use the same asymptotic variance estimation (Fisher's information), as shown 
#'  consistent with Tsiatis (1982) when considered under H0. One can choose the asymptotic variance
#'  estimation under H0 or under H1. The default is under H1 which is more important for
#'  weighted log-rank test because the weight function is usually based on the pooled data under H1.
#'  
#'  Different from Schoenfeld (1981) method, (He et al 2021) method considers 
#'  the effect as weighted average of each arm under H1. The difference to 
#'  Schoenfeld method is negligible under most common situations with HR large with log-rank test.
#'  However, when the assumed HR under H1 is small for weighted log-rank test, there is 
#'  notable small difference, in which scenario (He et al 2021) method can potentially 
#'  reduce sample size. The function calculates the asymptotic non-centrality 
#'  parameter for observable survival data at time of analysis, which has a maximum 
#'  of trial observational period of T, 
#'  calculated from first subject entry to the trial data cutoff date. This 
#'  function allows flexible alternative hypothesis in terms of HR(t), the 
#'  hazard ratio function over time. For delayed effect scenario under H1,
#'  one can define HR(t) as a piecewise constant function of survival time t.
#'  In addition, the function can handle user-defined flexible non-uniform enrollment 
#'  distribution function and independent time to lost-to-followup process which 
#'  is user-defined function of any lost-to-followup pattern such as constant
#'  lost-to-followup rate or Weibull distribution. For most common setting
#'  in practice, assuming the same lost-to-followup pattern in both arms.
#'  If the total number of subjects n is not provided, the function returns
#'  the non-centrality parameter of n^(-1/2)*Z where Z is the normalized 
#'  weighted log-rank test statistic Z = U/sqrt(var(U)) with U as the weighted
#'  log-rank score statistic. 
#'  
#' @param T  Calendar times, calculated from first subject randomization date, 
#'           when the weighted log-rank test will be evaluated at the same calendar time. 
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param n Total sample size for two arms. Default is NULL. 
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param  rho Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#' @param  gamma Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  s.tau  Survival rate cut S(tau) at t = tau1; default 0.5, ie. cut at median.
#'       s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  f.ws  Self-defined weight function of survival rate. 
#'         For example, f.ws = function(s){1/max(s, 0.25)}
#'         When f.ws is specified, the weight function takes them as priority.
#' @param F.entry Distribution function of enrollment. For uniform enrollment, 
#' F.entry(t) = (t/A) where A is the enrollment period, i.e., F.entry(t) = t/A for 0<=t<=A, and 
#' F.entry(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#' survival time is min(survival time, lost-to-follow-up time). Default G.ltfu = 0 (no lost-to-followup)
#' @param hypoth Hypothesis for the asymptotic variance estimation. The default is under H1 
#' which is more important for weighted log-rank test because the weight function 
#' is usually based on the pooled data under H1.
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  mu: Non-centrality parameter at calendar time, if n is not available.
#'  \item  R:  mu = n^(-1/2)*R
#'  \item  wt: Weight function used
#'  \item  hypoth: Hypothesis used
#'  }
#'  
#' @examples 
#' #Example (1) Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' #6 months after last patient randomized.
#' 
#' HR = 0.65; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' lambda1 = lambda0 * HR
#' f.logHR.PH = function(t){log(HR)}
#' h1.PH = function(t){lambda1}; S1.PH= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' rho = 0; gamma = 0; tau = NULL; s.tau = 0; f.ws = NULL; G.ltfu = function(t){0}
#' T = 24; r = 1; n = 450
#' 
#' wlr.mu.schoenfeld(T = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)
#'       
#' wlr.mu(T = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)
#'      
#' #Example (2) Same trial set up as example (1) but assuming delayed effect for 
#' experimental arm. The delayed period is assumed 6 months, and after delay the
#' hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' h1.D6 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1.D6 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HR)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G.ltfu = function(t){0}
#' 
#' wlr.mu.schoenfeld(T = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)
#'      
#' wlr.mu(T = 24, r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)
#' #Example (3). Draw a plot of non-centrality parameter over time    
#' 
#' plot.mu = function(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1, S1=S1, s.tau=0.5){
#'   maxT = 48; mu.lr = mu.fh = mu.fh11 = mu.sfh = rep(NA, maxT)
#'   for(i in 1:maxT){
#'     mu.lr[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     mu.fh[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'     mu.sfh[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = s.tau, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu    
#'     mu.fh11[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 1, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu          
#'   }
#'
#'   maxY = max(c(mu.lr, mu.fh, mu.sfh))
#'   plot(1:maxT, mu.lr, type="n", ylim=c(0, 4.5), 
#'   xlab="Months", ylab = "Non-centrality parameter")   
#'   lines(1:maxT, mu.lr, lty = 1, col=1)
#'   lines(1:maxT, mu.fh, lty = 2, col=2)
#'   lines(1:maxT, mu.sfh, lty = 3, col=3)
#'   lines(1:maxT, mu.fh11, lty = 4, col=4)
#'   abline(h = 1.96, lty = 2, col="gray")
#'      
#'   legend(0, 4.5, c("Log-rank", "FH(0,1)", paste("sFH(0,1,",s.tau,")", sep=""),"FH(1,1)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.8)
#' }
#' 
#' #(a)Proportional Hazards
#' plot.mu(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1.PH, S1=S1.PH)
#' 
#' #(b)Delayed effect of 6 months
#' plot.mu(f.logHR = f.logHR.D6, h0=h0, S0=S0, h1=h1.D6, S1=S1.D6, s.tau=0.5)
#' 
#' #(c)Crossover to effective subsequent IO therapy and delay 6 months
#' crossT = 30; HRx = 0.9; #after crossover, assuming the tail piecewise HR = 0.9.
#' h0x = function(t){lambda0*as.numeric(t < crossT) + HR/HRx*lambda0*as.numeric(t >= crossT)}; 
#' c0 = exp(-crossT*lambda0*(1-HR/HRx)); 
#' S0x = function(t){exp(-lambda0*t)*as.numeric(t<crossT) + c0*exp(-HR/HRx*lambda0*t)*as.numeric(t>=crossT)}
#' h1.D6x = function(t){lambda0*as.numeric(t < delay) + HR*lambda0*as.numeric(t >= delay)}
#' c1 = exp(-delay*lambda0*(1-HR)); 
#' S1.D6x = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c1*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' f.logHR.D6x = function(t){log(as.numeric(t < delay) 
#' + as.numeric(t >= delay & t < crossT)*HR + as.numeric(t >= crossT)*HRx)}
#' 
#' plot.mu(f.logHR = f.logHR.D6x, h0=h0x, S0=S0x, h1=h1.D6x, S1=S1.D6x, s.tau=0.5)
#' 
#' #Comparison of Schoenfeld(1981) and the new method.
#' plot.mu.comp = function(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1, S1=S1, s.tau=0.5){
#'   maxT = 48; mu.lr = mu.fh = sch.mu.lr = sch.mu.fh = rep(NA, maxT)
#'   for(i in 1:maxT){
#'     mu.lr[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     mu.fh[i] = wlr.mu(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'     sch.mu.lr[i] = wlr.mu.schoenfeld(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     sch.mu.fh[i] = wlr.mu.schoenfeld(T = i, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1, S1=S1, f.logHR = f.logHR,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'   }
#'
#'   plot(1:maxT, mu.lr, type="n", ylim=c(0, 4.5), 
#'   xlab="Months", ylab = "Non-centrality parameter")   
#'   lines(1:maxT, mu.lr, lty = 1, col=1)
#'   lines(1:maxT, mu.fh, lty = 2, col=1)
#'   lines(1:maxT, sch.mu.lr, lty = 1, col=2)
#'   lines(1:maxT, sch.mu.fh, lty = 2, col=2)
#'   abline(h = 1.96, lty = 2, col="gray")
#'   
#'   legend(0, 4.5, c("New: Log-rank", "New: FH(0,1)","Schoenfeld: Log-rank", "Schoenfeld: FH(0,1)"), 
#'   col=c(1,1,2,2), lty=c(1,2,1,2), bty="n", cex=0.8)
#' }
#' #(a)Proportional Hazards
#' plot.mu.comp(f.logHR = f.logHR.PH, h0=h0, S0=S0, h1=h1.PH, S1=S1.PH)
#' 
#' #(b)Delayed effect of 6 months
#' plot.mu.comp(f.logHR = f.logHR.D6, h0=h0, S0=S0, h1=h1.D6, S1=S1.D6, s.tau=0.5)
#' 
#' #(c)Crossover to effective subsequent IO therapy and delay 6 months 
#' plot.mu.comp(f.logHR = f.logHR.D6x, h0=h0x, S0=S0x, h1=h1.D6x, S1=S1.D6x, s.tau=0.5)
#' 
#' #Comparison of Schoenfeld(1981) and the new method at fixed time
#' #(d)Proportional Hazards
#'   HRseq = seq(0.4, 0.9, by=0.01)
#'   mu.lr = mu.fh = sch.mu.lr = sch.mu.fh = rep(NA, length(HRseq))
#'   for(i in 1:length(HRseq)){
#'     f.logHR.PH = function(t){log(HRseq[i])}
#'     h1.PH = function(t){lambda0 * HRseq[i]}; 
#'     S1.PH= function(t){exp(-lambda0 * HRseq[i] * t)}
#'     
#'     mu.lr[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     sch.mu.lr[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'        
#'     mu.fh[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'     sch.mu.fh[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.PH, S1=S1.PH, f.logHR = f.logHR.PH,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'   }
#'
#'   plot(HRseq, mu.lr, type="n", ylim=c(0, 10), 
#'   xlab="Hazard Ratio", ylab = "Non-centrality parameter at 48 months")   
#'   lines(HRseq, mu.lr, lty = 1, col=1)
#'   lines(HRseq, mu.fh, lty = 2, col=2)
#'   lines(HRseq, sch.mu.lr, lty = 3, col=3)
#'   lines(HRseq, sch.mu.fh, lty = 4, col=4)
#'   abline(h = 1.96, lty = 2, col="gray")
#'   
#'   legend(0.65, 10, c("New: Log-rank", "New: FH(0,1)","Schoenfeld: Log-rank", "Schoenfeld: FH(0,1)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.8)
#' 
#'   #(e)Delayed effect of 6 mo
#'   HRseq = seq(0.4, 0.9, by=0.01)
#'   mu.lr = mu.fh = sch.mu.lr = sch.mu.fh = rep(NA, length(HRseq))
#'   for(i in 1:length(HRseq)){
#'     h1.D6 = function(t){lambda0*as.numeric(t < delay)+HRseq[i]*lambda0*as.numeric(t >= delay)}
#'     c = exp(-delay*lambda0*(1-HRseq[i])); 
#'     S1.D6 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HRseq[i]*lambda0*t)*as.numeric(t>=delay)}
#'     f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HRseq[i])}
#'     mu.lr[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     sch.mu.lr[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'        
#'     mu.fh[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'     sch.mu.fh[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0, S0=S0,
#'        h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'   }
#'
#'   plot(HRseq, mu.lr, type="n", ylim=c(0, 10), 
#'   xlab="Hazard Ratio", ylab = "Non-centrality parameter at 48 months")   
#'   lines(HRseq, mu.lr, lty = 1, col=1)
#'   lines(HRseq, mu.fh, lty = 2, col=2)
#'   lines(HRseq, sch.mu.lr, lty = 3, col=3)
#'   lines(HRseq, sch.mu.fh, lty = 4, col=4)
#'   abline(h = 1.96, lty = 2, col="gray")
#'   
#'   legend(0.65, 10, c("New: Log-rank", "New: FH(0,1)","Schoenfeld: Log-rank", "Schoenfeld: FH(0,1)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.8)
#' 
#'   #(f)Crossover to effective subsequent IO therapy and delay 6 months 
#'   HRseq = seq(0.4, 0.9, by=0.01)
#'   mu.lr = mu.fh = sch.mu.lr = sch.mu.fh = rep(NA, length(HRseq))
#'   for(i in 1:length(HRseq)){
#'     crossT = 30; HRx = 0.9; #after crossover, assuming the tail piecewise HR = 0.9.
#'     h0x = function(t){lambda0*as.numeric(t < crossT) + HRseq[i]/HRx*lambda0*as.numeric(t >= crossT)}; 
#'     c0 = exp(-crossT*lambda0*(1-HRseq[i]/HRx)); 
#'     S0x = function(t){exp(-lambda0*t)*as.numeric(t<crossT) + c0*exp(-HRseq[i]/HRx*lambda0*t)*as.numeric(t>=crossT)}
#'     h1.D6x = function(t){lambda0*as.numeric(t < delay) + HRseq[i]*lambda0*as.numeric(t >= delay)}
#'     c1 = exp(-delay*lambda0*(1-HRseq[i])); 
#'     S1.D6x = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c1*exp(-HRseq[i]*lambda0*t)*as.numeric(t>=delay)}
#'     f.logHR.D6x = function(t){log(as.numeric(t < delay) 
#'        + as.numeric(t >= delay & t < crossT)*HRseq[i] + as.numeric(t >= crossT)*HRx)}
#'   
#'     mu.lr[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'     sch.mu.lr[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'        
#'     mu.fh[i] = wlr.mu(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'     sch.mu.fh[i] = wlr.mu.schoenfeld(T = 48, r = 1, n = 450, h0 = h0x, S0=S0x,
#'        h1 = h1.D6x, S1=S1.D6x, f.logHR = f.logHR.D6x,
#'        rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'        F.entry = F.entry, G.ltfu = G.ltfu)$mu      
#'   }
#'
#'   plot(HRseq, mu.lr, type="n", ylim=c(0, 10), 
#'   xlab="Hazard Ratio", ylab = "Non-centrality parameter at 48 months")   
#'   lines(HRseq, mu.lr, lty = 1, col=1)
#'   lines(HRseq, mu.fh, lty = 2, col=2)
#'   lines(HRseq, sch.mu.lr, lty = 3, col=3)
#'   lines(HRseq, sch.mu.fh, lty = 4, col=4)
#'   abline(h = 1.96, lty = 2, col="gray")
#'   
#'   legend(0.65, 10, c("New: Log-rank", "New: FH(0,1)","Schoenfeld: Log-rank", "Schoenfeld: FH(0,1)"), 
#'   col=1:4, lty=1:4, bty="n", cex=0.8)
#' 
#' @export
#' 
wlr.mu = function(T = 24, r = 1, n = 450, 
       h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
       h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
       f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
       rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
       F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
       G.ltfu = function(t){0}){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 

  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}  

  #UNDER H1
  f.bar = function(t){r0 * f0(t) + r1 * f1(t)}
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)}
  pi1 = function(t){r1 * S1(t) / S.bar(t)}
  
  #Weight function
  f.w = function(t, f.S = S0, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma){
    s = f.S(t)
    #First priority: f.ws
    if(!is.null(f.ws)){
      w = f.ws(s)
    }else {
      #Second priority: s.tau
      if (!is.null(s.tau)){
        s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max);
      } else {
        s.til = apply(cbind(s, f.S(tau)), MARGIN=1,FUN=max);        
      }
      w = s.til^rho*(1-s.til)^gamma
    }
    return(w)
  }
  
  #eta
  I.eta0 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (exp(f.logHR(t))-1)*F.entry(T-t) * (1 - G.ltfu(t)) * f0(t))
  }
  I.eta1 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (1-exp(-f.logHR(t)))*F.entry(T-t) * (1 - G.ltfu(t)) * f1(t))
  }
  #Use -log(HR(t)) as treatment effect to make mu in positive direction for 1-sided test
  D0= integrate(I.eta0, lower=0, upper=T, abs.tol=1e-8)$value
  D1= integrate(I.eta1, lower=0, upper=T, abs.tol=1e-8)$value
  eta = -sqrt(r1*r0)*(r0*D0 + r1*D1)
  
  #var
  I.var = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w^2 * F.entry(T-t) * (1 - G.ltfu(t)) * f.bar(t))
  }

  v = integrate(I.var, lower=0, upper=T, abs.tol=1e-8)$value

  #################
  #Exact calculation
  #################
  #eta
  I.eta0.2 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (exp(f.logHR(t))-1)*pi1(t)*(1-pi1(t))*F.entry(T-t) * (1 - G.ltfu(t)) * f0(t))
  }
  I.eta1.2 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w * (1-exp(-f.logHR(t)))*pi1(t)*(1-pi1(t))*F.entry(T-t) * (1 - G.ltfu(t)) * f1(t))
  }
  #Use -log(HR(t)) as treatment effect to make mu in positive direction for 1-sided test
  D0.2= integrate(I.eta0.2, lower=0, upper=T, abs.tol=1e-8)$value
  D1.2= integrate(I.eta1.2, lower=0, upper=T, abs.tol=1e-8)$value
  eta2 = -(r0*D0.2 + r1*D1.2)
  
  #var
  I.var2 = function(t){
    w = f.w(t, f.S = S.bar, f.ws=f.ws, tau=tau, s.tau=s.tau, rho=rho, gamma=gamma)
    return(w^2*pi1(t)*(1-pi1(t)) * F.entry(T-t) * (1 - G.ltfu(t)) * f.bar(t))
  }
  
  v2 = integrate(I.var2, lower=0, upper=T, abs.tol=1e-8)$value
  
  #R = n^(-1/2)*mu
  R = eta/sqrt(v)
  R2 = eta2/sqrt(v2)
  if (!is.null(n)){mu = sqrt(n * F.entry(T)) * R} else {mu = NULL}
  if (!is.null(n)){mu2 = sqrt(n * F.entry(T)) * R2} else {mu2 = NULL}

  o = list()
  o$mu = mu
  o$mu2 = mu2
  o$R = R
  o$R2 = R2
  
  if(!is.null(f.ws)){wt = f.ws} else{
    wt = data.frame(cbind(rho, gamma, tau, s.tau))
  }
  o$wt = wt
  return(o)
}
