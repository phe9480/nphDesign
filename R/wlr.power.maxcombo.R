#'  Power Calculation for Group Sequential Design Using A Max-combo Test
#' 
#'  This function calculates the power for group sequential design 
#'  based on the asymptotic
#'  distribution of the weighted log-rank test statistic under H1, with provided
#'  rejection boundaries.
#'  For group sequential design, the power will be calculated for each analysis
#'  and overall study. This function allows different max-combo tests for 
#'  for different analyses. 
#'  
#'  The non-centrality parameter can be chosen from Schoenfeld (1981) method, 
#'  or (He et al 2021) method (default).
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
#' @param  f.logHR Log(hazard ratio) function, = log(h1(t) / h0(t)).
#' @param  f.ws  Weight functions of survival rate used for maxcombo test at each analysis.
#'               For example (1), if there are 3 analyses planned including IA1, IA2, and FA.
#'               IA1 uses log-rank test, IA2 uses a maxcombo test of (logrank, FH01),
#'               and FA uses a maxcombo test of (logrank, FH01, FH11). Then specify as
#'               f.ws = list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)),
#'               where define lr = function(s){1}; fh01=function(s){1-s}; fh11 = function(s){s*(1-s)};
#'               For example (2), if only fh01 is used for all three analyses IA1, IA2 and FA.
#'               f.ws = list(IA1=list(fh01), IA2=list(fh01), FA1=list(fh01)).
#'               For example (3), if only logrank is used for a single time analysis, then
#'               f.ws = list(IA1=list(lr)).
#' @param F.entry Distribution function of enrollment. For example, 
#'        F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). 
#'        Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#'        survival time is min(survival time, lost-to-follow-up time). 
#'        Default G.ltfu = 0 (no lost-to-followup)
#' @param non.centrality  Option of non-centrality parameter either "Heetal2021"
#'        or "Schoenfeld". Default "Heetal2021". The default method is more efficient
#'        for HR < 0.65 in non-proportional hazards scenarios. No difference when HR > 0.65.
#' @param show.setting Option whether to show the recommended study design settings
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  mu:     Non-centrality parameter for each weighted log-rank test 
#'                 included in the max-combo.
#'  \item  events: Expected number of events at the calendar time T
#'  \item  maturity: Maturity of each arm for each analysis
#'  \item  T:      Expected calendar time for analysis (Data Cutoff Date)
#'  \item  power:  Power of the max-combo test
#'  \item  power.piece Power of each weighted log-rank test included in the max-combo test.
#'  \item  overall.power Overall power of the study
#'  \item  incremental.power Incremental power for each analysis. 
#'                  The sum of all incremental powers is the overall power.
#'  \item  n:       Total number of subjects for two arms
#'  \item  medians: Median of each treatment group
#'  \item  corr.Hp: Correlation matrix of the weighted log-rank tests included 
#'                  in the max-combo test for the pooled distribution (H0).
#'  \item  corr.strict.H0 Correlation matrix of the weighted log-rank tests included 
#'                  in the max-combo test assuming the data strictly follow
#'                  the control arm distribution (Strict H0).
#'  \item  wt: Weight functions used for each analysis
#'  \item  bounds: Rejection boundary for the max-combo test after multiplicity
#'                 adjustment based on asymptotic multivariate normal distribution.
#'                 b.Hp: bound based on pooled distribution (preferred); 
#'                 b.H0: bound based on strict H0
#'  \item  setting Setting of the distribution for each arm (hazard function, 
#'                 survival function, and log(HR) function), enrollment distribution
#'                 function, and lost-to-followup distribution.
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
#' #Define weight functions for weighted log-rank tests
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
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' G.ltfu = function(t){0}
#' 
#' wlr.power.maxcombo(T = 36, events = NULL, alpha=c(0.05), 
#'   power = 0.9, side = 2, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(FA=list(lr)), F.entry=F.entry, G.ltfu=G.ltfu)
#' 
#' #Equivalent to
#' 
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.05), 
#'   power = 0.9, side = 2, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(FA=lr), F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' #Also equivalent to
#' wlr.power(T = 36, events = NULL, alpha=c(0.05), 
#'   power = 0.9, side = 2, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(lr), F.entry=F.entry, G.ltfu=G.ltfu)
#' 
#' wlr.power.maxcombo(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(FA=list(lr, fh01)), F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' wlr.power.maxcombo(T = NULL, events = 230, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = 400, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(FA=list(lr, fh01)), F.entry=F.entry, G.ltfu=G.ltfu)   
#'   
#'   #Equivalent to
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(lr, fh01), F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' wlr.power.maxcombo(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(FA=list(lr, fh01, fh11)), F.entry=F.entry, G.ltfu=G.ltfu)
#'
#'   #Equivalent to   
#' wlr.power.maxcombo0(T = 36, events = NULL, alpha=c(0.025), 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(lr, fh01, fh11), F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#'   #Group sequential design
#'   #IA1: logrank; IA2: max(logrank, fh01); FA: max(lr, fh01, fh11)
#' wlr.power.maxcombo(T = c(24, 36, 48), events = NULL, 
#'   alpha=c(0.01, 0.02, 0.02)/2, 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(IA1 = list(lr), IA2 = list(lr, fh01), FA=list(lr, fh01, fh11)), 
#'   F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#'   #log-rank test for IA1, IA2, and FA
#' wlr.power.maxcombo(T = c(24, 36, 48), events = NULL, 
#'   alpha=c(0.01, 0.02, 0.02)/2, 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(IA1 = list(lr), IA2 = list(lr), FA=list(lr)), 
#'   F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#'   #Equivalent to 
#' wlr.power(T = c(24, 36, 48), events = NULL, 
#'   alpha=c(0.01, 0.02, 0.02)/2, 
#'   power = 0.9, side = 1, r = 1, n = NULL, rho=NULL, gamma=NULL, s.tau=NULL,
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(lr, lr, lr), 
#'   F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#'   #log-rank test for IA1, IA2, and FH01 for FA
#' wlr.power.maxcombo(T = c(24, 36, 48), events = NULL, 
#'   alpha=c(0.01, 0.02, 0.02)/2, 
#'   power = 0.9, side = 1, r = 1, n = NULL, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(IA1 = list(lr), IA2 = list(lr), FA=list(fh01)), 
#'   F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#' wlr.power.maxcombo(T = NULL, events = c(40, 60, 80), 
#'   alpha=c(0.01, 0.02, 0.02)/2, 
#'   power = NULL, side = 1, r = 1, n = 100, 
#'   h0 = h0, S0=S0,h1 = h.D3, S1= S.D3, f.logHR = f.logHR.D3, 
#'   f.ws = list(IA1 = list(lr), IA2 = list(lr), FA=list(fh01)), 
#'   F.entry=F.entry, G.ltfu=G.ltfu)
#'   
#'   
#' @export
#' 
wlr.power.maxcombo = function(T = c(24, 36, 48), events = NULL, 
            alpha=c(0.01, 0.02, 0.02)/2, power = 0.9, side = 1, r = 1, n = NULL, 
            h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
            h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
            f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
            f.ws = list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)),
            F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
            G.ltfu = function(t){0}, non.centrality = "Heetal2021", show.setting="N"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 

  #K analyses
  K=length(f.ws)

  #Number of test components in each analysis
  J = lengths(f.ws)
  
  #If n not provided but T provided. length(T) must be K.
  #Determine an initial n using log-rank test at FA
  if(is.null(n)){
    if(non.centrality != "Schoenfeld"){
      R = wlr.mu(T = T[K], r = r, n = NULL, h0 = h0, S0=S0,
             h1 = h1, S1=S1, f.logHR = f.logHR,
             rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
             F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }else{
      R = wlr.mu.schoenfeld(T = T[K], r = r, n = NULL, h0 = h0, S0=S0,
              h1 = h1, S1=S1, f.logHR = f.logHR,
              rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
              F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }
    if(side == 1){za = qnorm(1-sum(alpha))} else {za = qnorm(1-sum(alpha)/2)}
    n = ceiling(((za + qnorm(power))/R)^2 )
  }
  
  #If T is not provided but the events are provided. length(events) must be K.
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
  
  #If events are not provided, but T is provided. length(T) must be K.
  #Calculate the events according to T.
  events0 = events1 = rep(NA, K)
  if (is.null(events)){events = rep(NA, K)}
  
  for (i in 1:K){
    n.tmp = f.nEvents(T = T[i], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events
    
    if (is.na(events[i])){events[i] = n.tmp$n.events.total}
    events0[i] = n.tmp$n.events0
    events1[i] = n.tmp$n.events1
  }
    
  #Non-centrality parameters for all WLRTs in all analyses
  mu = matrix(NA, nrow=K, ncol=max(J))
  
  for(i in 1:K){
    for(j in 1:J[[i]]){
      if (non.centrality != "Schoenfeld") {
        non.centrality = "Heetal2021"
        mu[i, j] <- wlr.mu(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
            rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[i]][[j]],
            F.entry = F.entry, G.ltfu = G.ltfu)$mu2
      }else{
        tt = wlr.mu.schoenfeld(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
                                rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[i]][[j]],
                                F.entry = F.entry, G.ltfu = G.ltfu)
        mu[i, j] = c(tt$mu2)
      }
    }
  }

  #Complete correlation matrix structure: J[1]+J[2]+...+J[K] dimensions
  corr.Hp = corr.strict.H0 =matrix(1, nrow=sum(J), ncol=sum(J))
  #calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[[i]]){
      for (ip in i:K){
        for (jp in 1:J[ip]){
          row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the corr matrix
            
          #incremental location pointer for column compared to row of the corr matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            #information matrix for Zij and Zi'j'
            info = wlr.info(T = c(T[i], T[ip]), r = r, n = NULL, h0 = h0, S0= S0, h1 = h1, S1=S1, 
                     rho=NULL, gamma=NULL, tau=NULL, s.tau=NULL,
                     f.ws=c(f.ws[[i]][[j]],f.ws[[ip]][[jp]]), F.entry = F.entry, G.ltfu = G.ltfu)
            corr.Hp[row, col] = info$corr.Hp
            corr.strict.H0[row, col] = info$corr.H0
            corr.Hp[col, row]  = corr.Hp[row, col]
            corr.strict.H0[col, row] = corr.strict.H0[row, col]
          }
        }
      }
    }
  }
  
  #Rejection boundary: recursively solve for the rejection boundary for each analysis
  b.Hp = b.H0 = rep(NA, K)
  
  #First Analysis
  #maxcombo has at least 2 components.
  if (J[[1]] >= 2){
    if(side == 1){
      #1-sided test
      f.b.Hp = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
            mean=0, corr = corr.Hp[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      f.b.H0 = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
            mean=0, corr = corr.strict.H0[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }      
      b.Hp[1] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
      b.H0[1] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
    } else {
      #2-sided test
      f.b.Hp = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
            mean=0, corr = corr.Hp[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      f.b.H0 = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
            mean=0, corr = corr.strict.H0[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }      
      b.Hp[1] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
      b.H0[1] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
    }
  } else{
    if (side == 1) {b.Hp[1] = b.H0[1] = qnorm(1-alpha[1])} else {b.Hp[1] = b.H0[1] = qnorm(1-alpha[1]/2)}
  }
  
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
      if(side == 1){
        #1-sided test
        f.b.Hp = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b.Hp[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}

          idx1 = 1:sum(J[1:(i-1)])
          if (length(UU1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                    mean=rep(0, length(idx1)),corr = corr.Hp[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                  mean=rep(0, length(idx)), corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        f.b.H0 = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b.Hp[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}
          
          idx1 = 1:sum(J[1:(i-1)])
          if (length(LL1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                    mean=rep(0, length(idx1)), 
                    corr = corr.strict.H0[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                                mean=rep(0, length(idx)), 
                                corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        b.Hp[i] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0[i] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      } else {
        #2-sided test
        f.b.Hp = function(x){
          LL = rep(-b.Hp[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.Hp[j], J[j]))}}
          LL = c(LL, rep(x, J[i]))

          UU = rep(b.Hp[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.Hp[j], J[j]))}}
          UU = c(UU, rep(Inf, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          p1 = mvtnorm::pmvnorm(lower=LL, upper=UU, mean=rep(0, length(idx)), 
                                corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          LL2 = LL; UU2 = UU
          LL2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-Inf, J[i])
          UU2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-x, J[i])
          p2 = mvtnorm::pmvnorm(lower=LL2, upper=U2, mean=rep(0, length(idx)), 
                                corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(p1+p2-alpha[i])
        }
        f.b.H0 = function(x){
          LL = rep(-b.H0[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.H0[j], J[j]))}}
          LL = c(LL, rep(x, J[i]))
          
          UU = rep(b.H0[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.H0[j], J[j]))}}
          UU = c(UU, rep(Inf, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          p1 = mvtnorm::pmvnorm(lower=LL, upper=UU, mean=rep(0, length(idx)), corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          LL2 = LL; UU2 = UU
          LL2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-Inf, J[i])
          UU2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-x, J[i])
          p2 = mvtnorm::pmvnorm(lower=LL2, upper=U2, mean=rep(0, length(idx)), corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(p1+p2-alpha[i])
        }
        b.Hp[i] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0[i] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      }
    }
  }
  
  #Power calculation: Only consider superiority of 1-sided fashion
  power = rep(NA, K); #Power for each analysis

  #piecewise power for each weighted logrank test by each analysis using maxcombo rejection boundary
  power.piece = matrix(NA, nrow=K, ncol=max(J)); 
  
  for (i in 1:K){
    if(J[i] >= 2){
      mui = c(mu[i,]); mui = mui[!is.na(mui)]
      ix = (as.numeric(i>1)*sum(J[1:(i-1)])+1) : sum(J[1:i])
      power[i] = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,J[i]),upper=rep(b.Hp[i], J[i]), 
                                mean=mui, corr = corr.Hp[ix, ix], abseps = 1e-8, maxpts=100000)[1]
    } else {power[i] = 1 - pnorm(b.Hp[i], mean=mu[i, 1])}  
    
    for (j in 1:J[i]){
      power.piece[i, j] = 1-pnorm(b.Hp[i], mean=mu[i, j])
    }
  }
  
  #Overall power and incremental power
  overall.power = power[1]; incr.power = rep(0, K)
  
  if (K > 1){
    if (side == 1){
      for (i in 2:K){
        LL1 = rep(-Inf, sum(J[1:(i-1)]))
        UU1 = rep(b.Hp[1], J[1])
        if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}
        mu.i1 = mu[1,] 
        if (i > 2) {for (j in 2:(i-1)){mu.i1 = c(mu.i1, mu[j,])}}
        mu.i1 = mu.i1[!is.na(mu.i1)] #mean vector up to (i-1)th analysis
        
        idx1 = 1:sum(J[1:(i-1)])
        if (length(LL1) == 1) {
          P1 = pnorm(UU1, mean=mu.i1)
        } else {
          P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, mean=mu.i1, 
                                corr = corr.Hp[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
        }

        LL2 = rep(-Inf, sum(J[1:i]))
        UU2 = rep(b.Hp[1], J[1])
        for (j in 2:i){UU2 = c(UU2, rep(b.Hp[j], J[j]))}
        
        mu.i = mu[1,] 
        for (j in 2:i){mu.i = c(mu.i, mu[j,])}
        mu.i = mu.i[!is.na(mu.i)] #mean vector up to ith analysis
        
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, mean=mu.i, 
                              corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
        
        incr.power[i] = P1 - P2
        overall.power = overall.power + incr.power[i]
      }
    }else {
      for (i in 2:K){
        LL = rep(-b.Hp[1], J[1])
        if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.Hp[j], J[j]))}}
        LL = c(LL, rep(b.Hp[i], J[i]))
        
        UU = rep(b.Hp[1], J[1])
        if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.Hp[j], J[j]))}}
        UU = c(UU, rep(Inf, J[i]))
        
        mu.uptoi = mu[1,] 
        for (j in 2:i){mu.uptoi = c(mu.uptoi, mu[j,])}
        mu.uptoi = mu.uptoi[!is.na(mu.uptoi)] #mean vector up to ith analysis
        
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        incr.power[i] = mvtnorm::pmvnorm(lower = LL, upper = UU, mean=mu.uptoi, 
                                         corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
        overall.power = overall.power + incr.power[i]
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
  
  #######Critical Values#########
  CV.HR.H0 = CV.HR.H1 = rep(NULL, K)
  if (sum(J) == K){
    lr0 = 0
    for (i in 1:K){for (j in 1:J[i]){lr0 = f.ws[[i]][[j]](rnorm(1)) + lr0}}
    if (lr0 == K){
      #logrank test in all analyses

      for (i in 1:K){
        CV.HR.H0[i] = exp(-b.Hp[i]/sqrt(F.entry(T[i])*r0*r1*events[i]))
        r0.H1 = events0[i] / events[i]
        r1.H1 = 1 - r0.H1
        CV.HR.H1[i] = exp(-b.Hp[i]/sqrt(F.entry(T[i])*r0.H1*r1.H1*events[i]))
      }
    }
  }
  
  #only applicable to proportional hazards and using log-rank test
  ph = (f.logHR(1) == f.logHR(10) && f.logHR(10) == f.logHR(100))
  
  if (ph){
    CV.median.H0 = m0/CV.HR.H0
    CV.median.H1 = m0/CV.HR.H1
  } else{CV.median.H0 = CV.median.H1 = NULL}
  
  o$power.piece = power.piece
  o$overall.power = overall.power
  o$incremental.power = incr.power
  o$n = n
  
  o$medians = data.frame(round(cbind(m0, m1),2))
  o$Critical.Values.Medians.H0 = CV.median.H0
  o$Critical.Values.Medians.H1 = CV.median.H1
  
  o$Critical.Values.HR.H0 = CV.HR.H0
  o$Critical.Values.HR.H1 = CV.HR.H1
  
  o$corr.Hp = corr.Hp
  o$corr.strict.H0 = corr.strict.H0
  o$wt = f.ws
  o$bounds = data.frame(cbind(b.Hp, b.H0))
  #H1 setting of distributions;
  setting = list(r=r, alpha=alpha, side=side, h0=h0, S0=S0, h1=h1, S1=S1, log.HR = f.logHR, F.entry=F.entry, G.ltfu=G.ltfu, non.centrality=non.centrality)
  if (show.setting!="N"){o$setting = setting}
  
  #AHR using un-weighted Cox regression
  ahr = rep(NULL, K)
  for (i in 1:K){
  ahr[i] = wlr.AHR(T=T[i], r=r, n = n, h0=h0, S0=S0,h1=h1, S1 = S1,f.logHR = f.logHR,
                   rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws=NULL,
          F.entry = F.entry, G.ltfu = G.ltfu)$AHR
  }
  
  o$AHR = ahr
  return(o)
}

