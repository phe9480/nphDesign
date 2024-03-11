#'  Display Expected Average Weighted Hazard Ratio Over Calendar Time
#' 
#'  This function plots the expected average weighted HR under H1 over 
#'  calendar time. Two methods are implemented: (He et al 2021) method and 
#'  Kalbfleisch and Prentice (1981) method.
#'  Both methods are very close and the difference is usually negligible.
#'  
#'  This function allows flexible weight function as used in the weighted log-rank test. This 
#'  function also allows flexible alternative hypothesis in terms of HR(t), the 
#'  instaneous hazard ratio function over time. For delayed effect scenario under H1,
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
#' @param Tmax  Maximum range of calendar times. Default Tmax = 50.
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
#' @param  method Methods to calculate the AHR: Options include "geometric.schoenfeld",
#'                and "Kalbfleisch and Prentice". Default is "geometric.schoenfeld".
#'  
#' @examples 
#' #############
#' #Example (1) 
#' #############
#' #Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' #enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' #Control arm ~ exponential distribution with median 12 months, and 
#' #Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' #Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' #6 months after last patient randomized.
#' 
#' #Assuming delayed effect 6 months, and after delay the hazard ratio is assumed 0.65.
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
#' plot_AHR(r = 1, n = 450, h0 = h0, S0=S0,
#'      h1 = h1.D6, S1=S1.D6, f.logHR = f.logHR.D6,
#'      rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
#'      F.entry = F.entry, G.ltfu = G.ltfu, ylim=c(0.6, 1))
#'      
#' 
#' @export
#' 
plot_AHR = function(Tmax = 50, r = 1, n = 450, 
                   h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
                   h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                   f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
                   rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
                   F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                   G.ltfu = function(t){0}, method="geometric.schoenfeld", ...){
  
  t = seq(0, Tmax, by = 1); ahr = rep(NA, length(t))
  for (i in 1:length(t)){
    ahr.i=wlr.AHR(T=t[i], r=r, n = n, h0=h0, S0=S0,h1=h1, S1 = S1,f.logHR = f.logHR,
                  rho=rho, gamma=gamma, tau=tau, s.tau=s.tau, f.ws=f.ws,
                  F.entry = F.entry, G.ltfu = G.ltfu)
    if (method == "geometric.schoenfeld") {ahr[i] = ahr.i$AHR}
    if (method == "geometric.heetal2021") {ahr[i] = ahr.i$AHR2}
    if (method == "Kalbfleisch and Prentice") {ahr[i] = ahr.i$AHR.KP}
  }
  ahr = round(ahr, 3)
  plot(t, ahr, type="n", xlab="Calendar Time (mo)", ylab="Expected Average Hazard Ratio", ...)
  lines(t, ahr, lwd=3, col=1, lty=1)
  abline(h = seq(0, 1, 0.05), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  
}
