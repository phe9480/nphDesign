#'  Power Calculation at a Calendar Time Using A Max-combo Test
#' 
#'  This function calculates the power at a calendar time based on the asymptotic
#'  distribution of the weighted log-rank test statistic under H1, with provided
#'  rejection boundaries.
#'  For group sequential design, the power will be calculated for each analysis
#'  and overall study. 
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
#' @param T  Calendar time for analysis 
#' @param events Target number of events
#' @param alpha Allocated alpha level. Default 0.025 (1-sided).
#' @param power Power, default 0.9
#' @param side  Side of test 1 or 2. Default side = 1. 
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param n Total sample size for two arms. Default n = NULL. If n is not provided,
#'          n will be calculated based on single-time point analysis using 
#'          the first weighted log-rank test as a starting point with alpha level. 
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param  f.logHR Log(hazard ratio) function, = h1(t) / h0(t).
#' @param  rho Vector of rho parameters for multiple Fleming-Harrington (rho, gamma) 
#'             weighted log-rank tests.
#' @param  gamma Vector of parameters for Fleming-Harrington (rho, gamma) 
#'               weighted log-rank tests. For log-rank test, set rho = gamma = 0.
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
#'         must provide f.ws as a list of K-component weight functions or 
#'         (rho, gamma, tau, s.tau) with K components for each parameter.
#'         
#'         Note: The actual math formula for weighting function is based on 
#'         the left of each event time t, i.e., w(t) = f.ws(s(t-)). 
#'         For FH(rho, gamma) test, when gamma = 0, then the first event time has 
#'         weight 1;  when gamma > 0, the first event weight is 0. This can 
#'         ensure consistency with FH(0,0) = logrank; and FH(1,0) = generalized Wilcoxon.
#'         
#' @param F.entry Distribution function of enrollment. For example, 
#'        F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). 
#'        Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#'        survival time is min(survival time, lost-to-follow-up time). 
#'        Default G.ltfu = 0 (no lost-to-followup)
#' @param non.centrality  Option of non-centrality parameter either "Heetal2021"
#'        or "Schoenfeld". Default "Heetal2021". The default method is more efficient
#'        for HR < 0.65 in non-proportional hazards scenarios. No difference when HR > 0.65.
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  mu: Non-centrality parameter for each weighted log-rank test included in the max-combo.
#'  \item  events: Expected number of events at the calendar time T
#'  \item  T:      Expected calendar time for analysis (Data Cutoff Date)
#'  \item  power:  Power of the max-combo test
#'  \item  power.piece Power of each weighted log-rank test included in the max-combo test.
#'  \item  n:       Total number of subjects for two arms
#'  \item  medians: Median of each treatment group
#'  \item  corr.Hp: Correlation matrix of the weighted log-rank tests included 
#'                  in the max-combo test for the pooled distribution (H0).
#'  \item  corr.strict.H0 Correlation matrix of the weighted log-rank tests included 
#'                  in the max-combo test assuming the data strictly follow
#'                  the control arm distribution (Strict H0).
#'  \item  wt: Weight functions used for the weighted log-rank tests included
#'             in the max-combo test.
#'  \item  bounds: Rejection boundary for the max-combo test after multiplicity
#'                 adjustment based on asymptotic multivariate normal distribution.
#'                 b.Hp: bound based on pooled distribution (preferred); 
#'                 b.H0: bound based on strict H0
#'  \item  setting Setting of the distribution for each arm (hazard function, 
#'                 survival function, and log(HR) function), enrollment distribution
#'                 function, and lost-to-followup distribution.
#'  \item  non-centrality: Method of non-centrality parameter
#'  }
#'  
#' @examples 
#' #Distributions for both arms
#' m0 = 12; #median RFS for control arm
#' lambda0 = log(2) / m0
#' h0 = function(t){lambda0}; 
#' S0 = function(t){exp(-lambda0 * t)}
#' HRd = 0.60 #hazard ratio after delay
#' 
#' h.D3=function(t){lambda0*as.numeric(t<3)+HRd*lambda0*as.numeric(t>=3)}
#' c3 = exp(-3*lambda0*(1-HRd)); 
#' S.D3 = function(t){S0(t)*as.numeric(t<3)+c3*exp(-HRd*lambda0*t)*as.numeric(t>=3)}
#' f.logHR.D3 = function(t){log(as.numeric(t<3) + as.numeric(t>= 3)*HRd)}
#' 
#' h.D6=function(t){lambda0*as.numeric(t<6)+HRd*lambda0*as.numeric(t>=6)}
#' c6 = exp(-6*lambda0*(1-HRd)); 
#' S.D6 = function(t){exp(-lambda0*t)*as.numeric(t<6)+c6*exp(-HRd*lambda0*t)*as.numeric(t>=6)}
#' f.logHR.D6 = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*HRd)}
#' 
#' #Define weight funtions for weighted log-rank tests
#' lr = function(s){1}
#' fh01 = function(s){(1-s)}
#' fh11 = function(s){s*(1-s)}
#' #stabilized FH(0, 1; 0.5)
#' sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 
#' #modestly log-rank
#' mfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1/s1)}
#' 
#' #Define max-combo test
#' maxcombo1 = list(lr)
#' maxcombo2 = list(lr, fh01)
#' maxcombo3 = list(lr, fh01, fh11)
#' 
#' maxcombo4 = list(lr, sfh01)
#' maxcombo5 = list(lr, mfh01)
#' maxcombo6 = list(lr, sfh01, mfh01, fh01)
#' 
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.05), 
#'   power = 0.9, side = 2, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo1)
#' 
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo2)
#'   
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo3)
#'   
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo4)
#'   
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo5)
#'   
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, f.ws = maxcombo6)
#'   
#' @export
#' 
#' 
wlr.power.maxcombo0 = function(T = 36, events = NULL, 
                     alpha=0.025, power = 0.9, side = 1, r = 1, n = NULL, 
            h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
            h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
            f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
            f.ws = list(lr, fh01, fh11),
            F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
            G.ltfu = function(t){0}, non.centrality = "Heetal2021"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Determine sample size n using log-rank test if not provided but requires T as initial value
  if(is.null(n)){
    if(non.centrality != "Schoenfeld"){
      R = wlr.mu(T = T, r = r, n = NULL, h0 = h0, S0=S0,
               h1 = h1, S1=S1, f.logHR = f.logHR,
               rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
               F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }else{
      R = wlr.mu.schoenfeld(T = T, r = r, n = NULL, h0 = h0, S0=S0,
              h1 = h1, S1=S1, f.logHR = f.logHR,
              rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
              F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }
    if(side == 1){za = qnorm(1-sum(alpha))} else {za = qnorm(1-sum(alpha)/2)}
    n = ceiling(((za + qnorm(power))/R)^2 )
  }
  n0=n
  #Determine the calendar times for the required number of events
  if(!is.null(events)){
    T = rep(NA, length(events))
    for (i in 1:length(events)){
      f.root = function(t){
        events[i] - f.nEvents(T = t, r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events$n.events.total
      }
      T[i] = uniroot(f.root, interval= c(1, 1000), tol = 1e-8)$root
    }
  }
  
  #K analyses
  K=length(T)

  #calculate the events  
  events0 = events1 = rep(NA, K)
  if (is.null(events)){events = rep(NA, K)}
  
  for (i in 1:length(T)){
    n.tmp = f.nEvents(T = T[i], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events
    
    if (is.na(events[i])){events[i] = ceiling(n.tmp$n.events.total)}
    events0[i] = ceiling(n.tmp$n.events0)
    events1[i] = ceiling(n.tmp$n.events1)
  }
    
  #Determine number of WLRTs at each analysis
  #J[1]: number of WLRTs at 1st analysis
  J = length(f.ws)
  
  #Non-centrality parameters for all WLRTs in all analyses
  mu = matrix(NA, nrow=K, ncol=max(J))
  
  for(i in 1:K){
    for(j in 1:J[i]){
      if (non.centrality != "Schoenfeld") {
        non.centrality = "Heetal2021"
        mu[i, j] = wlr.mu(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
            rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[j]],
            F.entry = F.entry, G.ltfu = G.ltfu)$mu2
      }else{
        mu[i, j] = wlr.mu.schoenfeld(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
            rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[j]],
            F.entry = F.entry, G.ltfu = G.ltfu)$mu2
      }
    }
  }

  #Rejection boundary
  #Only consider 1 analysis now
    if(K==1){
      corr.Hp = corr.strict.H0 =matrix(1, nrow=J[1], ncol=J[1])
      if (J[1] >= 2){
        #correlation matrix of (z11, ..., z_1J[1])
        for (i in 1:(J[1]-1)){
          for (j in (i+1):J[1]){
            info = wlr.info(T = T[1], r = r, n = NULL, h0 = h0, S0= S0, h1 = h1, S1=S1, 
                     rho=NULL, gamma=NULL, tau=NULL, s.tau=NULL,
                     f.ws=c(f.ws[[i]],f.ws[[j]]), F.entry = F.entry, G.ltfu = G.ltfu)
            corr.Hp[i, j] = info$corr.Hp
            corr.strict.H0[i, j] = info$corr.H0
            corr.Hp[j, i] = corr.Hp[i, j]
            corr.strict.H0[j, i] = corr.strict.H0[i, j]
          }
        }
      }
    }
  
    #maxcombo has at least 2 components.
    if (J[1] >= 2){
      if(side == 1){
        #1-sided test
        f.b.Hp = function(x){
          1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.Hp, abseps = 1e-8, maxpts=100000)[1] - alpha
        }
        f.b.H0 = function(x){
          1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.strict.H0, abseps = 1e-8, maxpts=100000)[1] - alpha
        }      
        b.Hp = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0 = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      } else {
        #2-sided test
        f.b.Hp = function(x){
          1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
                             mean=0, corr = corr.Hp, abseps = 1e-8, maxpts=100000)[1] - alpha
        }
        f.b.H0 = function(x){
          1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
                             mean=0, corr = corr.strict.H0, abseps = 1e-8, maxpts=100000)[1] - alpha
        }      
        b.Hp = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0 = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      }
    } else{
      if (side == 1) {b.Hp = b.H0 = qnorm(1-alpha)} else {b.Hp = b.H0 = qnorm(1-alpha/2)}
    }
  
  
  #Power calculation: Only consider superiority of 1-sided fashion
  power = rep(NA, K); power.piece = matrix(NA, nrow=K, ncol=max(J))
  for (i in 1:K){
    if(J[1] >= 2){
      power[i] = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(b.Hp, J[1]), 
                                mean=mu[i,], corr = corr.Hp, abseps = 1e-8, maxpts=100000)[1]
    } else {power[i] = 1- pnorm(b.Hp, mean=mu[i, 1])}  
    for (j in 1:J[i]){
      if(side == 1){power.piece[i, j] = 1- pnorm(qnorm(1-alpha), mean=mu[i, j])}else{
        power.piece[i, j] = 1- pnorm(qnorm(1-alpha/2), mean=mu[i, j])
      }
    }
  }
  
  #Calculate the medians
  f.m0 = function(t){S0(t) - 0.5}
  f.m1 = function(t){S1(t) - 0.5}
  
  m0 = uniroot(f.m0, interval= c(1, 100), tol = 1e-8)$root
  m1 = uniroot(f.m1, interval= c(1, 100), tol = 1e-8)$root
  
  all.events = data.frame(cbind(events0, events1, events))
  
  maturity0 = events0/(r0*n);
  maturity1 = events1/(r1*n);
  overall.maturity = events/n
  
  maturity = data.frame(cbind(maturity0, maturity1, overall.maturity))
  
  o = list()
  o$mu = mu
  o$events = all.events
  o$maturity = maturity
  
  o$T = T
  o$power = power
  o$power.piece = power.piece
  o$n = n
  
  o$medians = data.frame(cbind(m0, m1))
  
  o$corr.Hp = corr.Hp
  o$corr.strict.H0 = corr.strict.H0
  o$wt = f.ws
  o$bounds = data.frame(cbind(b.Hp, b.H0))
  #H1 setting of distributions;
  setting = list(r=r, alpha=alpha, side=side, h0=h0, S0=S0, h1=h1, S1=S1, log.HR = f.logHR, F.entry=F.entry, G.ltfu=G.ltfu)
  o$setting = setting
  o$non.centrality = non.centrality
  return(o)
}

