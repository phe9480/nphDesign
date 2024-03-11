#'  Power Calculation at a Calendar Time Using A Weighted Log-rank Test
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
#' @param T  Calendar times for analysis, calculated from first subject randomization date, 
#'           for the interim analyses and final analysis. If T is not provided, 
#'           events must be provided, so the calendar times T can be calculated.
#' @param events Target numbers of events for all analyses. If events are provided,
#'               the analysis time T will be calculated. 
#' @param alpha Allocated alpha levels. Required parameter if b is not provided. 
#'              The alphas add up to the total type I error for the group sequential 
#'              design of the endpoint. 
#'           If alpha spending function a(t) is used for information time c(t1, ..., tK),
#'           then alpha1 = a(t1), alpha2 = a(t2)-a(t1), ..., alphaK = a(tK)-a(t_{K-1}),
#'           and the total alpha for all analyses is a(tK).
#' @param b  Rejection boundary in normalized Z. If b is NULL, then the boundaries 
#'           will be calculated based on alpha, side of test at each analysis time T. Default b = NULL.
#'           If b is provided, then alpha is ignored. Default NULL.
#' @param r  Randomization ratio of experimental arm : control arm as r:1. 
#'           When r = 1, it is equal allocation. Default r = 1.
#' @param n Total sample size for two arms. Default n = NULL. If n is not provided,
#'          n will be calculated based on single-time point analysis using 
#'          the first weighted log-rank test as a starting point
#'          with total alpha level, sum(alpha). 
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
#'         s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
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
#'  \item  incremental.power Incremental power for each analysis. 
#'                  The sum of all incremental powers is the overall power.
#'  \item  T:      Calendar times of analyses
#'  \item  n:      Total sample size
#'  \item  bounds: Rejection boundary of Z statistic for power calculation
#'  \item  events: Target events of analyses
#'  \item  Critical.Values.Medians CV of median
#'  \item  Critical.Values.HR CV of HR
#'  \item  Critical.Values.Weighted CV of weighted HR from weighted Cox Regression
#'  \item  setting: Distribution setting under alternative hypothesis
#'  \item  corr:   Correlation matrix of weighted log-rank tests
#'  \item  wt: Weight function used
#'  \item  H0.type: Type of H0 used. Default is "pooled.H0". 
#'                  The other option is "strict.H0", which assumes the experimental
#'                  arm exactly follows the control arm.
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
#' r = 1; n = 450; events = c(264, 330)
#' HR = 6.7/12.5; lambda0 = log(2) / 6.7; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' lambda1 = lambda0 * HR
#' f.logHR.PH = function(t){log(HR)}
#' h1.PH = function(t){lambda1}; 
#' S1.PH= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' F.entry.uni = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G.ltfu = function(t){0}
#' 
#' #(a) log-rank test is used for both IA and FA 
#' rho = gamma = s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' wlr.power(T = c(24, 30), events = NULL, b = NULL, 
#' alpha=c(0.02, 0.03), power = 0.9, side = 2, r = 1, n = NULL,
#' h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, f.logHR = f.logHR.PH, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry.uni, G.ltfu = G.ltfu, non.centrality = "Heetal2021")
#' 
#' #(b) log-rank is used for IA and FH01 is used for FA
#' rho = c(0,0); gamma = c(0, 1); s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' 
#' wlr.power(T = NULL, events = events, b = c(qnorm(1-0.024/2), qnorm(1-0.042/2)), 
#' r = r, n = n, alpha=c(0.02, 0.03),
#' h0 = h0, S0= S0, h1= h1.PH, S1=S1.PH, f.logHR = f.logHR.PH, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu, non.centrality = "Heetal2021")
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
#' 
#' #(a) log-rank test is used for both IA and FA 
#' rho = gamma = s.tau = rep(0, 2); tau = NULL; f.ws = NULL; events=c(397, 496)
#' wlr.power(T = NULL, events = events, b = NULL, r = 1, n = 672,
#' h0 = h0, S0= S0, h1= h.D6, S1=S.D6, f.logHR = f.logHR.D6, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu, non.centrality = "Heetal2021")
#' 
#' #(b) log-rank is used for IA and FH01 is used for FA
#' rho = c(0,0); gamma = c(0, 1); s.tau = rep(0, 2); tau = NULL; f.ws = NULL; 
#' 
#' wlr.power(T = NULL, events = events, b = c(qnorm(1-0.024/2), qnorm(1-0.042/2)), r = r, n = n,
#' h0 = h0, S0= S0, h1= h1.D6, S1=S1.D6, f.logHR = f.logHR.D6, 
#' rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
#' F.entry = F.entry, G.ltfu = G.ltfu, non.centrality = "Heetal2021")
#' 
#' @export
#' 
#' 
wlr.power = function(T = c(24, 36), events = NULL, b = NULL, 
                     alpha=c(0.02, 0.03), power = 0.9, side = 2, r = 1, n = NULL, 
            h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
            h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
            f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
            rho = c(0,0), gamma = c(0,0), tau = c(NULL,NULL), s.tau = c(0, 0), f.ws = NULL,
            F.entry = function(t){(t/21)*as.numeric(t <= 21) + as.numeric(t > 21)}, 
            G.ltfu = function(t){0}, non.centrality = "Heetal2021", H0.type="pooled.H0"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1  
  
  #Determine sample size n if not provided but requires T as initial value
  if(is.null(n)){
    if(non.centrality != "Schoenfeld"){
      R = wlr.mu(T = T[length(T)], r = r, n = NULL, h0 = h0, S0=S0,
               h1 = h1, S1=S1, f.logHR = f.logHR,
               rho = rho[1], gamma = gamma[1], tau = tau[1], s.tau = s.tau[1], f.ws = f.ws[[1]],
               F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }else{
      R = wlr.mu.schoenfeld(T = T[length(T)], r = r, n = NULL, h0 = h0, S0=S0,
                 h1 = h1, S1=S1, f.logHR = f.logHR,
                 rho = rho[1], gamma = gamma[1], tau = tau[1], s.tau = s.tau[1], f.ws = f.ws[[1]],
                 F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }
    if(side == 1){za = qnorm(1-sum(alpha))} else {za = qnorm(1-sum(alpha)/2)}
    n = ceiling(((za + qnorm(power))/R)^2 )
  }
  
  #Determine the calendar times T for the required number of events
  if(!is.null(events)){
    T = rep(NA, length(events))
    for (i in 1:length(events)){
      f.root = function(t){
        events[i] - f.nEvents(T = t, r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                            F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events$n.events.total
      }
      T[i] = uniroot(f.root, interval= c(1, 200), tol = 1e-8)$root
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
    
  #Non-centrality parameters
  mu = rep(NA, K)
  
  for(i in 1:K){
    if (non.centrality != "Schoenfeld") {
      non.centrality = "Heetal2021"
      mu[i] = wlr.mu(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
              rho = rho[i], gamma = gamma[i], tau = tau[i], s.tau = s.tau[i], f.ws = f.ws[[i]],
              F.entry = F.entry, G.ltfu = G.ltfu)$mu2
    }else{
      mu[i] = wlr.mu.schoenfeld(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
              rho = rho[i], gamma = gamma[i], tau = tau[i], s.tau = s.tau[i], f.ws = f.ws[[i]],
              F.entry = F.entry, G.ltfu = G.ltfu)$mu2
    }
  }

  #Rejection boundary
  if(is.null(b)){
    bd.tmp = wlr.bounds.design(T = T, r = r, alpha=alpha, side=side,
            h0 = h0, S0= S0, h1= h1, S1=S1, 
            rho = rho, gamma = gamma, tau =tau, s.tau =s.tau, f.ws = f.ws,
            F.entry = F.entry, G.ltfu = G.ltfu)
    if(H0.type != "strict.H0"){
      b = bd.tmp$bounds.pooled.H0
    } else {b = bd.tmp$bounds.strict.H0}
  }
  
  #Power calculation
  power = rep(NA, K)
  for (i in 1:K){
    power[i] = 1-pnorm(b[i], mean=mu[i])
  }
  
  #Correlation matrix under H1
  corr = matrix(1, nrow=K, ncol=K)
  if(K > 1){
    #correlation matrix of (z1, ..., zK)
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        corr[i, j] = wlr.info(T = c(T[i], T[j]), r = r, n = n, h0 = h0, S0= S0, h1 = h1, S1=S1, 
                rho=c(rho[i], rho[j]), gamma=c(gamma[i],gamma[j]), 
                tau = c(tau[i],tau[j]), s.tau=c(s.tau[i],s.tau[j]),
                f.ws=list(f.ws[[i]], f.ws[[j]]), F.entry = F.entry, G.ltfu = G.ltfu)$corr.Hp
        corr[j, i] = corr[i, j]
      }
    }
  }  
  #Overall power and incremental power
  overall.power = power[1]; incr.power = rep(0, K)
  if(K > 1) {
    for(i in 2:K){
      incr.power[i] = mvtnorm::pmvnorm(lower = c(rep(-Inf, i-1), b[i]), 
                                       upper = c(b[1:(i-1)], Inf), mean=mu[1:i], 
                                       corr = corr[1:i, 1:i], abseps = 1e-8, maxpts=100000)[1]
      overall.power = overall.power + incr.power[i]
    }  
  }
  
  #Calculate the medians
  f.m0 = function(t){S0(t) - 0.5}
  f.m1 = function(t){S1(t) - 0.5}
  
  m0 = uniroot(f.m0, interval= c(1, 100), tol = 1e-8)$root
  m1 = uniroot(f.m1, interval= c(1, 100), tol = 1e-8)$root
  
  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}  
  
  #UNDER H1
  f.bar = function(t){r0 * f0(t) + r1 * f1(t)}
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)}
  
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
  
  CV.HR.H0 = CV.HR.H1 = rep(NA, K)
  for (i in 1:K){  
    #Asymptotic variance of AHR
    I.I0 = function(t){
      w = f.w(t, f.S = S.bar, f.ws=f.ws[[i]], tau=tau[i], s.tau=s.tau[i], rho=rho[i], gamma=gamma[i])
      return(w^2 * F.entry(T[i]-t) * (1 - G.ltfu(t)) * f.bar(t))
    }
    I0 = integrate(I.I0, lower=0, upper=T[i], abs.tol=1e-8)$value
  
    #d: Prob. of event
    I.d = function(t){
      w = f.w(t, f.S = S.bar, f.ws=f.ws[[i]], tau=tau[i], s.tau=s.tau[i], rho=rho[i], gamma=gamma[i])
      return(w * F.entry(T[i]-t) * (1 - G.ltfu(t)) * f.bar(t))
    }
    d = integrate(I.d, lower=0, upper=T[i], abs.tol=1e-8)$value
  
    V.H0 = I0 / (n*F.entry(T[i])*r0*r1*d^2)
    
    r0.H1 = events0[i] / (events0[i]+events1[i])
    r1.H1 = 1 - r0.H1
    V.H1 = I0 / (n*F.entry(T[i])*r0.H1*r1.H1*d^2)
    
    CV.HR.H0[i] = exp(-b[i]*sqrt(V.H0))
    CV.HR.H1[i] = exp(-b[i]*sqrt(V.H1))
    
  }
  
  #only applicable to proportional hazards and using log-rank test
  ph = (f.logHR(1) == f.logHR(10) && f.logHR(10) == f.logHR(100))
  
  f.ws.lr = rho.gamma.lr = 1
  for (i in 1:K) {
    f.ws.lr.i = (!is.null(f.ws) && f.ws[[i]](1) == 1 && f.ws[[i]](10) == 1 && f.ws[[i]](100) == 1)
    f.ws.lr = f.ws.lr * f.ws.lr.i
    
    rho.gamma.lr.i = (!is.null(rho) && rho[i]==0 && gamma[i] == 0 )
    rho.gamma.lr = rho.gamma.lr * rho.gamma.lr.i
  }
  
  if (ph && (f.ws.lr || rho.gamma.lr)){
    CV.median.H0 = m0/CV.HR.H0
    CV.median.H1 = m0/CV.HR.H1
  } else{CV.median.H0 = CV.median.H1 = NULL}
  

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
  o$incremental.power = incr.power
  o$overall.power = overall.power
  
  o$n = n
  
  o$medians = data.frame(cbind(m0, m1))
  
  #only applicable to proportional hazards and using log-rank test
  
  o$Critical.Values.Medians.H0 = CV.median.H0
  o$Critical.Values.Medians.H1 = CV.median.H1

  #o$Critical.Values.HR.H0 = CV.HR.H0
  o$Critical.Values.Weighted.H0 = CV.HR.H0
  
  o$Critical.Values.Weighted.HR.H1 = CV.HR.H1
  
  o$corr = corr
  if(!is.null(f.ws)){wt = f.ws} else{
    wt = data.frame(cbind(rho, gamma, tau, s.tau))
  }
  
  o$wt = wt
  o$bounds = b
  
  #H1 setting of distributions;
  setting = list(alpha = alpha, side = side, h0=h0, S0=S0, h1=h1, S1=S1, log.HR = f.logHR, F.entry=F.entry, G.ltfu=G.ltfu)
  o$setting = setting
  o$H0.type = H0.type
  o$non.centrality = non.centrality
  return(o)
}

