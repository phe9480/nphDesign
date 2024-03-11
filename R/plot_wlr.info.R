#' Display Fisher's Information Weighted Log-rank Tests Over Calendar Time
#' 
#' This function plots the Fisher's information for 
#' weighted log-rank tests alternative hypotheses based on the provided 
#' enrollment distribution function and random lost-to-followup distribution if applicable.
#' 
#' 
#' @param T  A vector of two calendar times, calculated from first subject randomization date. 
#'           If the two calendar times are the same or just one calendar time is provided,
#'            the two weighted log-rank tests will be evaluated at the same calendar time. 
#' @param r  Randomization ratio of experimental arm : control arm as r:1. When r = 1, it is equal allocation. Default r = 1.
#' @param h0 Hazard function of control arm. h0(t) = log(2)/m0 means T~exponential distribution with median m0.
#' @param S0 Survival function of control arm. In general, S0(t) = exp(- integral of h0(u) for u from 0 to t).
#'           but providing S0(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f0(t) = h0(t) * S0(t).
#' @param h1 Hazard function of experimental arm. h1(t) = log(2)/m1 means T~exponential distribution with median m0.
#' @param S1 Survival function of experimental arm. In general, S1(t) = exp(- integral of h1(u) for u from 0 to t).
#'           but providing S1(t) can improves computational efficiency and 
#'           usually the survival function is known in study design. The density function f1(t) = h1(t) * S1(t).
#' @param F.entry Distribution function of enrollment. For uniform enrollment, 
#' F.entry(t) = (t/A) where A is the enrollment period, i.e., F.entry(t) = t/A for 0<=t<=A, and 
#' F.entry(t) = 1 when t > A. For more general non-uniform enrollment with weight psi, 
#' F.entry(t) = (t/A)^psi*I(0<=t<=A) + I(t>A). Default F.entry is uniform distribution function.
#' @param G.ltfu Distribution function of lost-to-follow-up censoring process. The observed
#' survival time is min(survival time, lost-to-follow-up time). Default G.ltfu = 0 (no lost-to-followup)
#' @param  f.ws  Self-defined weight function of survival rate. 
#'         For example, f.ws1 = function(s){1/max(s, 0.25)}
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  cov.H0: Covariance matrix of n**(- 1 / 2) * (U1, U2) under H1, 
#'  where U1 and U2 are weighted log-rank score statistics.
#'  \item  info.H0: Fisher's information under H0. H0: survival distribution follows control arm.
#'  \item  corr.H0: Correlation between U1 and U2 under H0
#'  \item  cov.H1: Covariance matrix of n**(- 1 / 2) * (U1, U2) under H1, 
#'  where U1 and U2 are weighted log-rank score statistics.
#'  \item  info.H1: Fisher's information under H1. H1: survival distribution 
#'  follows weighted average of control arm and experimental arm, i.e., 
#'  S = p1 * S1 + p0 * S0.
#'  \item  corr.H1: Correlation between U1 and U2 under H1
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
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' #(a) 1 weighted log-rank test at 1 analysis time
#' plot_wlr.info(Tmax = 50, r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' fws = list(function(s){1-s}, function(s){1}),
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' 
#' @export

plot_wlr.info = function(Tmax = 50, r = 1, n = 450, h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)}, 
          h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
          fws = list(function(s){1-s}), F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, G.ltfu = function(t){0}){
  t = seq(0.5, Tmax, by = 1)
  info = matrix(NA, nrow=length(t), ncol=length(fws))
  for (i in 1:length(t)){
    for (j in 1:length(fws)){
      info[i, j] = wlr.info(T = t[i], r = r, h0 = h0, S0= S0, h1 = h1, S1=S1, 
              rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = list(fws[[j]], NULL),
              F.entry = F.entry, G.ltfu = G.ltfu, n = n)$info.Hp
    }
  }
  plot(t, info[,1], type="n", ylim=c(0, max(info)), xlab="Calendar Time (mo)", ylab = "Fisher's Information")
  for (j in 1:length(fws)){
    lines(t, info[,j], lwd=3, col=j)
  }
  
  abline(h = seq(0, max(info), 5), col="gray80", lty=3)
  abline(v = seq(0, Tmax, 2), col="gray80", lty=3)
  
}
