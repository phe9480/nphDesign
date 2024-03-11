#'  Display Non-centrality Parameter Under H1 For A Weighted Log-rank Test
#' 
#'  This function plots the non-centrality parameter under H1 evaluated at a calendar time with observable data. Both Schoenfeld (1981) 
#'  and this method use the same asymptotic variance estimation (Fisher's information), as shown 
#'  consistent with Tsiatis (1982) when considered under H0. 
#'   
#'  
#' @param Tmax  Range of calendar time
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
#'  
#' @examples 
#' 
#' plot_wlr.mu (Tmax = 50, r = 1, n = 450, 
#' h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
#' h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
#' f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
#' fws = list(function(s){1-s}), F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
#' G.ltfu = function(t){0})
#' 
#' @export
#' 
plot_wlr.mu = function(Tmax = 50, r = 1, n = 450, 
       h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
       h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
       f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
       fws = list(function(s){1-s}), F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
       G.ltfu = function(t){0}){
  
  t = seq(0.5, Tmax, by = 1)
  mu = matrix(NA, nrow=length(t), ncol=length(fws))
  for (i in 1:length(t)){
    for (j in 1:length(fws)){
      mu[i, j] = wlr.mu(T = t[i], r = r, h0 = h0, S0= S0, h1 = h1, S1=S1, 
            rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = fws[[j]],
            F.entry = F.entry, G.ltfu = G.ltfu, n = n)$mu
    }
  }
  plot(t, mu[,1], type="n", ylim=c(0, max(mu)), xlab="Calendar Time (mo)", 
       ylab = "Non-centrality Parameter")
  for (j in 1:length(fws)){
    lines(t, mu[,j], lwd=3, col=j)
  }
  abline(h = seq(0, max(mu), 0.1), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  
}
