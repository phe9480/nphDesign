#'  Estimated Rejection Boundary Per Design of Weighted Log-rank Test at a Calendar Time in Group Sequential Design
#' 
#'  This function calculates the rejection boundary at a calendar time based on the asymptotic
#'  distribution of the weighted log-rank test statistic under H0, with provided
#'  alpha spending approach to achieve the desired power. For group sequential design, 
#'  the power will also be calculated for each analysis for the provided calendar time of analysis.
#'   
#'  
#'  The non-centrality parameter can be chosen from Schoenfeld (1981) method, or (He et al 2021) method (default)
#'  The difference is usually negligible under most common situations with HR large.
#'  However, when the assumed HR under H1 is small for weighted log-rank test, there is 
#'  notable small difference, in which scenario (He et al 2021) method can potentially 
#'  reduce sample size. Regarding the timing of analyses, either the calendar time
#'  or required number of total events can be specified. By default, the required
#'  number of total events will be used as determination of analysis cut off. 
#'  This function allows flexible alternative hypothesis in terms of HR(t), the 
#'  hazard ratio function over time. For delayed effect scenario under H1,
#'  one can define HR(t) as a piecewise constant function of survival time t.
#'  In addition, the function can handle user-defined flexible non-uniform enrollment 
#'  distribution function and independent time to lost-to-followup process which 
#'  is user-defined function of any lost-to-followup pattern such as constant
#'  lost-to-followup rate or Weibull distribution. For most common setting
#'  in practice, assuming the same lost-to-followup pattern in both arms.
#'  
#' @param T  Calendar times, calculated from first subject randomization date, 
#'           for the interim analyses and final analysis. 
#' @param events Required numbers of events for all analyses         
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
#'         
#'         Either f.ws or (rho, gamma, tau, s.tau) must be specified. For K analyses,
#'         must provide f.ws as a list of K-component weight functions or (rho, gamma, tau, s.tau) with K components for each parameter.
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
#'  \item  power:  Power at each analysis
#'  \item  overall.power Overall power for the study
#'  \item  T:      Calendar times of analyses
#'  \item  events: Target events of analyses
#'  \item  corr:   Correlation matrix of weighted log-rank tests
#'  \item  wt: Weight function used
#'  \item  hypoth: Hypothesis used
#'  \item  non-centrality: Method of non-centrality parameter
#'  }
#'  
#' @examples 
#' #Example (1) Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' #6 months after last patient randomized.
#' #IA and FA are planned at 225 and 330 events respectively.
#' 
#' r = 1; 
#' HR = 0.65; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' lambda1 = lambda0 * HR
#' h1.PH = function(t){lambda1}; 
#' S1.PH= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G.ltfu = function(t){0}
#' 
#' #(a) log-rank test is used for both IA and FA 
#' rho = gamma = s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' wlr.bounds.design(T = c(30, 42), r = r, alpha=c(0.02, 0.03),
#' h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu)
#' 
#' #(b) log-rank is used for IA and FH01 is used for FA
#' rho = c(0,0); gamma = c(0, 1); s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' 
#' wlr.bounds.design(T = c(30, 42), r = r, alpha=c(0.02, 0.03),
#' h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu)
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
#' 
#' #(a) log-rank test is used for both IA and FA 
#' rho = gamma = s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' tmp = wlr.bounds.design(T = c(30, 42), r = r, alpha=c(0.02, 0.03),
#' h0 = h0, S0= S0, h1= h1.D6, S1=S1.D6, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu)
#' b = tmp$bounds.pooled.H0
#' 
#' wlr.power(T = NULL, events = c(397,496), alpha=c(0.024, 0.025), 
#' b = NULL, r = r, n = n,
#' h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, f.logHR = f.logHR.PH, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, side=2,
#' G.ltfu = G.ltfu, H0.type="pooled.H0", non.centrality = "Heetal2021")
#' 
#' #(b) log-rank is used for IA and FH01 is used for FA
#' rho = c(0,0); gamma = c(0, 1); s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' 
#' wlr.bounds.design(T = c(30, 42), r = r, alpha=c(0.02, 0.03),
#' h0 = h0, S0= S0, h1= h1.D6, S1=S1.D6, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu)
#' 
#' #Example (3). Plot of boundary over time for FA given IA fixed with time and alpha level
#' rho=c(0,0); gamma = c(0,1);
#' 
#' #(a)Proportional Hazards
#' h1=h1.PH; S1 = S1.PH; 
#' 
#' #(b)Delayed effect of 6 months
#' h1=h1.D6; S1 = S1.D6; 
#' 
#' maxT = 60; T.IA = 30; lrfh.bd.PH = lrfh.bd.D6 =  lrlr.bd.PH = lrlr.bd.D6 = rep(NA, maxT - T.IA)
#' fh.PH.mu = fh.D6.mu = lr.PH.mu = lr.D6.mu = rep(NA, maxT - T.IA)
#' f.logHR.PH = function(t){log(HR)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HR)}
#' 
#'   for(j in (T.IA+1):maxT){
#'     fh.PH.mu[j-T.IA] = wlr.mu(T = j, r = r, n = 450, h0 = h0, S0=S0, h1 = h1.PH, S1=S1.PH, 
#'     f.logHR = f.logHR.PH, rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'      
#'     fh.D6.mu[j-T.IA] = wlr.mu(T = j, r = r, n = 450, h0 = h0, S0=S0, h1 = h1.D6, S1=S1.D6, 
#'     f.logHR = f.logHR.D6, rho = 0, gamma = 1, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'      
#'     lr.PH.mu[j-T.IA] = wlr.mu(T = j, r = r, n = 450, h0 = h0, S0=S0, h1 = h1.PH, S1=S1.PH, 
#'     f.logHR = f.logHR.PH, rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'      
#'     lr.D6.mu[j-T.IA] = wlr.mu(T = j, r = r, n = 450, h0 = h0, S0=S0, h1 = h1.D6, S1=S1.D6, 
#'     f.logHR = f.logHR.D6, rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu)$mu
#'      
#'     lrfh.bd.PH[j-T.IA] = wlr.bounds.design(T = c(T.IA, j), r = r, alpha=c(0.02, 0.03),
#'                   h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, rho = c(0,0), gamma = c(0,1), 
#'                   tau =NULL, s.tau =s.tau, f.ws = NULL,
#'                   F.entry = F.entry, G.ltfu = G.ltfu)$bounds.pooled.H0[2]
#'      
#'     lrfh.bd.D6[j-T.IA] = wlr.bounds.design(T = c(T.IA, j), r = r, alpha=c(0.02, 0.03),
#'                   h0 = h0, S0= S0, h1= h1.D6, S1=S1.D6, rho = c(0,0), gamma = c(0,1), 
#'                   tau =NULL, s.tau =s.tau, f.ws = NULL,
#'                   F.entry = F.entry, G.ltfu = G.ltfu)$bounds.pooled.H0[2]
#'     lrlr.bd.PH[j-T.IA] = wlr.bounds.design(T = c(T.IA, j), r = r, alpha=c(0.02, 0.03),
#'                   h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, rho = c(0,0), gamma = c(0,0), 
#'                   tau =NULL, s.tau =s.tau, f.ws = NULL,
#'                   F.entry = F.entry, G.ltfu = G.ltfu)$bounds.pooled.H0[2]
#'     lrlr.bd.D6[j-T.IA] = wlr.bounds.design(T = c(T.IA, j), r = r, alpha=c(0.02, 0.03),
#'                   h0 = h0, S0= S0, h1= h1.D6, S1=S1.D6, rho = c(0,0), gamma = c(0,0), 
#'                   tau =NULL, s.tau =s.tau, f.ws = NULL,
#'                   F.entry = F.entry, G.ltfu = G.ltfu)$bounds.pooled.H0[2]
#'   }
#'  
#'   plot((T.IA+1):maxT, lrlr.bd.PH, type="n", ylim=c(0, 5), 
#'   xlab="Months", ylab = "FA Rejection Boundary")   
#'   lines((T.IA+1):maxT, lrlr.bd.PH, lty = 1, col=1)
#'   lines((T.IA+1):maxT, lrlr.bd.D6, lty = 2, col=1)
#'   lines((T.IA+1):maxT, lrfh.bd.PH, lty = 1, col=2)
#'   lines((T.IA+1):maxT, lrfh.bd.D6, lty = 2, col=2)
#'   
#'   lines((T.IA+1):maxT, lr.PH.mu, lty = 1, col=1)
#'   lines((T.IA+1):maxT, lr.D6.mu, lty = 2, col=1)
#'   lines((T.IA+1):maxT, fh.PH.mu, lty = 1, col=2)
#'   lines((T.IA+1):maxT, fh.D6.mu, lty = 2, col=2)
#'   
#'   legend(45, 2.1, c("lr-lr: PH", "lr-lr: Delay 6", 
#'   "lr-fh: PH", "lr-fh: Delay 6"), 
#'   col=c(1,1,1,2,2,2), lty=c(1,2,3,1,2,3), bty="n", cex=0.8)
#' 
#' @export
#' 
#' 
wlr.bounds.design = function(T = c(24, 36), r = 1, alpha = c(0.01, 0.04), side = 2,
    h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
    h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
    rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0, 0), f.ws = NULL,
    F.entry = function(t){(t/21)*as.numeric(t <= 21) + as.numeric(t > 21)}, 
    G.ltfu = function(t){0}, hypoth="H.Pooled"){
  
  #K analyses
  K=length(T)

  #Correlation matrix 
  corr.Hp = corr.strict.H0 =matrix(1, nrow=K, ncol=K)
  if(K > 1){
    #correlation matrix of (z1, ..., zK)
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        info = wlr.info(T = c(T[i], T[j]), r = r, n = NULL, h0 = h0, S0= S0, h1 = h1, S1=S1, 
                rho=c(rho[i], rho[j]), gamma=c(gamma[i], gamma[j]), tau=c(tau[i],tau[j]), s.tau=c(s.tau[i],s.tau[j]),
                f.ws=c(f.ws[[i]],f.ws[[j]]), F.entry = F.entry, G.ltfu = G.ltfu)
        corr.Hp[i, j] = info$corr.Hp
        corr.strict.H0[i, j] = info$corr.H0
        corr.Hp[j, i] = corr.Hp[i, j]
        corr.strict.H0[j, i] = corr.strict.H0[i, j]
      }
    }
  }  

  #Recursively solve for each boundary
  f.solver.b = function(corr = corr){
   b = rep(NA, K)
   if (side == 1) {b[1] = qnorm(1-alpha[1])} else {b[1] = qnorm(1-alpha[1]/2)}
  
   if(K > 1){
    if(side == 1){
      for (i in 2:K){
        f.b = function(x){
          mvtnorm::pmvnorm(lower = c(rep(-Inf, i-1), x), 
                           upper = c(b[1:(i-1)], Inf), mean=rep(0, i), 
                           corr = corr[1:i, 1:i], abseps = 1e-8, maxpts=100000)[1] - alpha[i]
        }
        b[i] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
      }
    } else {
      for (i in 2:K){
        f.b = function(x){
          I1 = mvtnorm::pmvnorm(lower = c(-b[1:(i-1)], x), upper = c(b[1:(i-1)], Inf), mean=rep(0, i), corr = corr[1:i, 1:i], abseps = 1e-8, maxpts=100000)[1]
          I2 = mvtnorm::pmvnorm(lower = c(-b[1:(i-1)], -Inf), upper = c(b[1:(i-1)], -x), mean=rep(0, i), corr = corr[1:i, 1:i], abseps = 1e-8, maxpts=100000)[1]
          return(I1+I2-alpha[i])
        }
        b[i] = uniroot(f=f.b, lower=1, upper=20, tol = 1e-8)$root
      }
    }  
   }  
   return(b)
  }
  #H0 (pooled samples), use weighted survival distribution (f.bar = p0*f0+p1*f1) 
  #to estimate the pooled samples
  b.pooled.H0 = f.solver.b (corr = corr.Hp)
  #Strict H0: assuming f1 = f0.
  b.strict.H0 = f.solver.b (corr = corr.strict.H0)
  o = list()
  o$bounds.pooled.H0 = b.pooled.H0
  o$bounds.strict.H0 = b.strict.H0
  o$T = T
  o$corr.pooled.H0 = corr.Hp
  o$corr.strict.H0 = corr.strict.H0
  if(!is.null(f.ws)){wt = f.ws} else{
    wt = data.frame(cbind(rho, gamma, tau, s.tau))
  }
  o$wt = wt
  
  return(o)
}

