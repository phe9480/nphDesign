#' Fisher's Information Matrix of Two Weighted Log-rank Tests At Different Calendar Time
#' 
#' This function calculates the Fisher's information matrix for two 
#' weighted log-rank tests under null and alternative hypotheses
#' at two calendar times, which are counted from first subject randomized. The function
#' returns the information matrix, based on the provided enrollment distribution function 
#' and random lost-to-followup distribution if applicable.
#' 
#' \itemize{
#' \item If the only one weighted log-rank test is provided, the Fisher information
#' will be calculated only for the provided weighted log-rank test. A weighted
#' log-rank test specification requires the explicit weight function.
#' \item If the total sample size is not provided, then the covariance matrix of n^-1/2*(U1, U2)
#' will be provided, where U1 and U2 are the corresponding weighted log-rank scores.
#' U1 and U2 can be based on two different calendar times. A combination test can be 
#' constructed accordingly based on the asymptotic multivariate
#' normal distribution of n^-1/2*(U1, U2) when U1 and U2 are from 
#' the same calendar time, i.e. the same analysis.
#' }
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
#' @param  rho1 Parameter for Fleming-Harrington (rho1, gamma1) weighted log-rank test.
#' @param  gamma1 Parameter for Fleming-Harrington (rho1, gamma1) weighted log-rank test.
#'         For log-rank test, set rho1 = gamma1 = 0.
#' @param  tau1  Cut point for stabilized FH test, sFH(rho1, gamma1, tau1); with weight
#'       function defined as w1(t) = s_tilda1^rho1*(1-s_tilda1)^gamma1, where
#'       s_tilda1 = max(s(t), s.tau1) or max(s(t), s(tau1)) if s.tau1 = NULL
#'       tau1 = Inf reduces to regular Fleming-Harrington test(rho1, gamma1)
#' @param  s.tau1  Survival rate cut S(tau1) at t = tau1; default 0.5, ie. cut at median.
#'       s.tau1 = 0 reduces to regular Fleming-Harrington test(rho1, gamma1)
#' @param  f.ws1  Self-defined weight function of survival rate. 
#'         For example, f.ws1 = function(s){1/max(s, 0.25)}
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#' @param  rho2 Parameter for Fleming-Harrington (rho2, gamma2) weighted log-rank test.
#' @param  gamma2 Parameter for Fleming-Harrington (rho2, gamma2) weighted log-rank test.
#'         For log-rank test, set rho2 = gamma2 = 0.
#' @param  tau2  Cut point for stabilized FH test, sFH(rho2, gamma2, tau2); with weight
#'       function defined as w2(t) = s_tilda2^rho2*(1-s_tilda2)^gamma2, where
#'       s_tilda2 = max(s(t), s.tau2) or max(s(t), s(tau2)) if s.tau2 = NULL
#'       tau2 = Inf reduces to regular Fleming-Harrington test(rho2, gamma2)
#' @param  s.tau2  Survival rate cut S(tau2) at t = tau2; default 0.5, ie. cut at median.
#'       s.tau2 = 0 reduces to regular Fleming-Harrington test(rho2, gamma2)
#' @param  f.ws2  Self-defined weight function of survival rate. 
#'         For example, f.ws2 = function(s){1/max(s, 0.25)}. 
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#' @param n Total sample size for two arms. Default is NULL. 
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
#' wlr.info(T = c(24), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0), gamma = c(0), tau = NULL, s.tau = c(0), f.ws = NULL,
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' #(b) 2 weighted log-rank tests at 1 analysis time
#' wlr.info(T = c(24), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' #(c) 1 weighted log-rank test at 2 analysis times. If run 
#' #the same weighted log-rank test at each analysis time separately, 
#' #then no correation will be produced.
#' wlr.info(T = c(24, 36), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0), gamma = c(0), tau = NULL, s.tau = c(0), f.ws = NULL,
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' #(d) 2 weighted log-rank tests at 2 analysis times. 
#' #Log-rank and FH01 at 24 and 36 months respectively.
#' wlr.info(T = c(24,36), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' #Same results using f.ws2() function to specify the F-H weight
#' wlr.info(T = c(24,36), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0,0), f.ws = list(function(s){1},function(s){(1-s)}),
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' #(e) Draw plot of information for log-rank test and FH01 for T in (0, 48)
#' maxT = 48
#' info.lr.H0 = info.lr.H1 =info.fh.H0 =info.fh.H1 =info.sfh.H0 =info.sfh.H1 = rep(NA, maxT)
#' for (i in 1:maxT) {
#'   lr = wlr.info(T = i, r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(0), tau = NULL, s.tau = c(0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   info.lr.H0[i] = lr$info.H0; info.lr.H1[i] = lr$info.Hp;
#'   
#'   fh = wlr.info(T = i, r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(1), tau = NULL, s.tau = c(0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   info.fh.H0[i] = fh$info.H0; info.fh.H1[i] = fh$info.Hp;
#'   
#'   sfh = wlr.info(T = i, r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(1), tau = NULL, s.tau = c(0.5), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   info.sfh.H0[i] = sfh$info.H0; info.sfh.H1[i] = sfh$info.Hp;
#' }
#' 
#' plot(1:maxT, info.lr.H0/info.lr.H0[maxT], type="n", ylim=c(0, 1), 
#' xlab="Months", ylab = "Information Fraction")   
#' lines(1:maxT, info.lr.H0/info.lr.H0[maxT], lty = 1, col=1)
#' lines(1:maxT, info.lr.H1/info.lr.H1[maxT], lty = 2, col=1)
#' lines(1:maxT, info.fh.H0/info.fh.H0[maxT], lty = 1, col=2)
#' lines(1:maxT, info.fh.H1/info.fh.H1[maxT], lty = 2, col=2)
#' lines(1:maxT, info.sfh.H0/info.sfh.H0[maxT], lty = 1, col=3)
#' lines(1:maxT, info.sfh.H1/info.sfh.H1[maxT], lty = 2, col=3) 
#' legend(0, 1, c("Log-rank(H0)", "Log-rank(H1)", "FH01(H0)", "FH01(H1)", "sFH01(H0)", "sFH01(H1)"), col=c(1,1,2,2,3,3), lty=c(1,2,1,2,1,2), bty="n", cex=0.8)
#' 
#' #(f) Draw plot of correlation between IA and FA over the timing of IA
#' maxT = 48
#' corr.lr.H0 = corr.lr.H1 =corr.fh.H0 =corr.fh.H1 =corr.sfh.H0 =corr.sfh.H1 = rep(NA, maxT)
#' for (i in 1:maxT) {
#'   lr = wlr.info(T = c(i, maxT), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(0), tau = NULL, s.tau = c(0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.lr.H0[i] = lr$corr.H0; corr.lr.H1[i] = lr$corr.Hp;
#'   
#'   fh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(1), tau = NULL, s.tau = c(0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.fh.H0[i] = fh$corr.H0; corr.fh.H1[i] = fh$corr.Hp;
#'   
#'   sfh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0), gamma = c(1), tau = NULL, s.tau = c(0.5), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.sfh.H0[i] = sfh$corr.H0; corr.sfh.H1[i] = sfh$corr.Hp;
#' }
#' 
#' plot(1:maxT, corr.lr.H0, type="n", ylim=c(0, 1), 
#' xlab="Months", ylab = "Correlation btw IA and FA")   
#' lines(1:maxT, corr.lr.H0, lty = 1, col=1)
#' lines(1:maxT, corr.lr.H1, lty = 2, col=1)
#' lines(1:maxT, corr.fh.H0, lty = 1, col=2)
#' lines(1:maxT, corr.fh.H1, lty = 2, col=2)
#' lines(1:maxT, corr.sfh.H0, lty = 1, col=3)
#' lines(1:maxT, corr.sfh.H1, lty = 2, col=3)
#' legend(0, 1, c("Log-rank(H0)", "Log-rank(H1)", "FH01(H0)", "FH01(H1)", "sFH01(H0)", "sFH01(H1)"), col=c(1,1,2,2,3,3), lty=c(1,2,1,2,1,2), bty="n", cex=0.8)
#' 
#' #(g) Draw a plot to show the theoretical calculation between correlation and sqrt(information fraction)
#' #Before enrollment complete, the correlation > sqrt(Information Fraction).
#' #After enrollment complete, they are equivalent.
#' plot(1:maxT, corr.lr.H0, type="n", ylim=c(0, 1), 
#' xlab="Months", ylab = "Correlation btw IA and FA")   
#' lines(1:maxT, corr.lr.H0, lty = 1, col=1)
#' lines(1:maxT, sqrt(info.lr.H0/info.lr.H0[maxT]), lty = 2, col=1)
#' lines(1:maxT, corr.fh.H0, lty = 1, col=2)
#' lines(1:maxT, sqrt(info.fh.H0/info.fh.H0[maxT]), lty = 2, col=2)
#' lines(1:maxT, corr.sfh.H0, lty = 1, col=3)
#' lines(1:maxT, sqrt(info.sfh.H0/info.sfh.H0[maxT]), lty = 2, col=3)
#' legend(0, 1, c("Log-rank(H0): corr", "Log-rank(H0): sqrt(IF)", "FH01(H0): corr", "FH01(H0): sqrt(IF)", "sFH01(H0): corr", "sFH01(H0): sqrt(IF)"), col=c(1,1,2,2,3,3), lty=c(1,2,1,2,1,2), bty="n", cex=0.8)
#' 
#' #(h) Draw plot of correlation between IA and FA for different choices of tests at IA and FA
#' maxT = 48
#' corr.lr.lr.H0 =corr.lr.fh.H0 = corr.lr.sfh.H0 =  rep(NA, maxT)
#' corr.fh.lr.H0 =corr.fh.fh.H0 =corr.fh.sfh.H0 = rep(NA, maxT)
#' corr.sfh.lr.H0 =corr.sfh.fh.H0 =corr.sfh.sfh.H0 = rep(NA, maxT)
#' for (i in 1:maxT) {
#'   lr.lr = wlr.info(T = c(i, maxT), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(0,0), tau = NULL, s.tau = c(0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.lr.lr.H0[i] = lr.lr$corr.H0;
#'   
#'   lr.fh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.lr.fh.H0[i] = lr.fh$corr.H0; 
#'   
#'   lr.sfh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0, 0.5), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.lr.sfh.H0[i] = lr.sfh$corr.H0; 
#'      
#'   fh.lr = wlr.info(T = c(i, maxT), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(1,0), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.fh.lr.H0[i] = fh.lr$corr.H0;
#'   
#'   fh.fh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(1,1), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.fh.fh.H0[i] = fh.fh$corr.H0; 
#'   
#'   fh.sfh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(1,1), tau = NULL, s.tau = c(0,0.5), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.fh.sfh.H0[i] = fh.sfh$corr.H0; 
#'      
#'   sfh.lr = wlr.info(T = c(i, maxT), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,0), gamma = c(1,0), tau = NULL, s.tau = c(0.5, 0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.sfh.lr.H0[i] = sfh.lr$corr.H0;
#'   
#'   sfh.fh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,1), gamma = c(0,1), tau = NULL, s.tau = c(0.5, 0), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.sfh.fh.H0[i] = sfh.fh$corr.H0; 
#'   
#'   sfh.sfh = wlr.info(T = c(i, 48), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#'   rho = c(0,1), gamma = c(0,1), tau = NULL, s.tau = c(0.5, 0.5), f.ws = NULL,
#'   F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#'   corr.sfh.sfh.H0[i] = sfh.sfh$corr.H0;     
#' }
#' 
#' plot(1:maxT, corr.lr.H0, type="n", ylim=c(0, 1), 
#' xlab="Months", ylab = "Correlation btw IA and FA")   
#' lines(1:maxT, corr.lr.lr.H0, lty = 1, col=1)
#' lines(1:maxT, corr.lr.fh.H0, lty = 2, col=1)
#' lines(1:maxT, corr.lr.sfh.H0, lty = 3, col=1)
#' 
#' lines(1:maxT, corr.fh.lr.H0, lty = 1, col=2)
#' lines(1:maxT, corr.fh.fh.H0, lty = 2, col=2)
#' lines(1:maxT, corr.fh.sfh.H0, lty = 3, col=2)
#' 
#' lines(1:maxT, corr.sfh.lr.H0, lty = 1, col=3)
#' lines(1:maxT, corr.sfh.fh.H0, lty = 2, col=3)
#' lines(1:maxT, corr.sfh.sfh.H0, lty = 3, col=3)
#' legend(0, 1, c("LR-LR(H0)", "LR-FH(H0)", "LR-sFH(H0)", 
#' "FH-LR(H0)", "FH-FH(H0)", "FH-sFH(H0)",
#' "sFH-LR(H0)", "sFH-FH(H0)", "sFH-sFH(H0)"), col=c(1,1,1,2,2,2,3,3,3), 
#' lty=c(1:3,1:3,1:3), bty="n", cex=0.8)
#' 
#' 
#' #Example (2) Same trial set up as example (1) but assuming delayed effect for 
#' experimental arm. The delayed period is assumed 6 months, and after delay the
#' hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' #Hazard function and survival function for experimental arm
#' h1 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' #2 weighted log-rank tests at 2 analysis times. 
#' #Log-rank and FH01 at 24 and 36 months respectively.
#' wlr.info(T = c(24,36), r = 1, h0 = h0, S0= S0, h1 = h1, S1=S1, 
#' rho = c(0,0), gamma = c(0,1), tau = NULL, s.tau = c(0,0), f.ws = NULL,
#' F.entry = F.entry, G.ltfu = function(t){0}, n = 450)
#' 
#' @export

wlr.info = function(T = c(24,36), r = 1, n = 450, h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)}, 
          h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
          rho = c(0,0), gamma = c(0,0), tau = NULL, s.tau = c(0,0), f.ws = NULL,
          F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, G.ltfu = function(t){0}){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}

  #Weight function
  f.w = function(t, f.S = S0, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1]){
    s = f.S(t)
    #First priority: f.ws
    if(!is.null(f.ws0)){
      w = f.ws0(s)
    }else {
      #Second priority: s.tau
      if (!is.null(s.tau0)){
        s.til = apply(cbind(s, s.tau0), MARGIN=1,FUN=max);
      } else {
        s.til = apply(cbind(s, f.S(tau0)), MARGIN=1,FUN=max);        
      }
      w = s.til^rho0*(1-s.til)^gamma0
    }
    return(w)
  }

  #sigma matrix (covariance)
  f.cov = function(f = f0, f.S = S0) {
   #(1) Only 1 analysis time
   if(length(T)==1 || (length(T)==2 && T[2] == T[1])){
  
    #(1.1) Only 1 weighted log-rank test at 1 time
    if(is.null(f.ws[[2]]) && (length(rho)==1 || is.null(rho))) {
      I11 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        return(w1^2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      sigma11 = r1*r0*integrate(I11, lower=0, upper=T[1], abs.tol=1e-8)$value
      cov = sigma11
      out = list()
      out$cov = cov; 
      info=NULL
      if (!is.null(n)){info = n * F.entry(T[1]) * cov}
      out$info = info
      return(out)
    } else {
      
      #(1.2) Two weighted log-rank tests at 1 time
      I11 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        return(w1^2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      I12 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        w2 = f.w(t, f.S = f.S, f.ws0=f.ws[[2]], tau0=tau[2], s.tau0=s.tau[2], rho0=rho[2], gamma0=gamma[2])
        return(w1 * w2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      I22 = function(t){
        w2 = f.w(t, f.S = f.S, f.ws0=f.ws[[2]], tau0=tau[2], s.tau0=s.tau[2], rho0=rho[2], gamma0=gamma[2])
        return(w2^2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      sigma11 = r1*r0*integrate(I11, lower=0, upper=T[1], abs.tol=1e-8)$value
      sigma12 = r1*r0*integrate(I12, lower=0, upper=T[1], abs.tol=1e-8)$value
      sigma22 = r1*r0*integrate(I22, lower=0, upper=T[1], abs.tol=1e-8)$value
      
      cov = matrix(NA, nrow=2, ncol=2)
      cov[1,1]=sigma11; cov[1,2] = cov[2,1] = sigma12; cov[2,2]=sigma22
      out = list()
      out$cov = cov
      info=NULL
      if (!is.null(n)){info = n * F.entry(T[1]) * cov}
      out$info = info
      return(out)
    }
  } 
  #(2) Two analysis times scenarios
  if (length(T)==2 && T[2] > T[1]) {
    #(2.1) One weighted log-rank test at 2 times
    if(is.null(f.ws[[2]]) && (length(rho)==1 || is.null(rho))){
      I11 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        return(w1^2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      I22 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        return(w1^2 * F.entry(T[2]-t) * (1 - G.ltfu(t)) * f(t))
      }
      sigma11 = r1*r0*integrate(I11, lower=0, upper=T[1], abs.tol=1e-8)$value
      sigma12 = sigma11
      sigma22 = r1*r0*integrate(I22, lower=0, upper=T[2], abs.tol=1e-8)$value
    } else {
    #(2.2) Two weighted log-rank tests at 2 times
      I11 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        return(w1^2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      I12 = function(t){
        w1 = f.w(t, f.S = f.S, f.ws0=f.ws[[1]], tau0=tau[1], s.tau0=s.tau[1], rho0=rho[1], gamma0=gamma[1])
        w2 = f.w(t, f.S = f.S, f.ws0=f.ws[[2]], tau0=tau[2], s.tau0=s.tau[2], rho0=rho[2], gamma0=gamma[2])
        return(w1 * w2 * F.entry(T[1]-t) * (1 - G.ltfu(t)) * f(t))
      }
      I22 = function(t){
        w2 = f.w(t, f.S = f.S, f.ws0=f.ws[[2]], tau0=tau[2], s.tau0=s.tau[2], rho0=rho[2], gamma0=gamma[2])
        return(w2^2 * F.entry(T[2]-t) * (1 - G.ltfu(t)) * f(t))
      }
      sigma11 = r1*r0*integrate(I11, lower=0, upper=T[1], abs.tol=1e-8)$value
      sigma12 = r1*r0*integrate(I12, lower=0, upper=T[1], abs.tol=1e-8)$value
      sigma22 = r1*r0*integrate(I22, lower=0, upper=T[2], abs.tol=1e-8)$value
    }
    cov = matrix(NA, nrow=2, ncol=2)
    cov[1,1]=sigma11; cov[1,2] = cov[2,1] = sigma12; cov[2,2]=sigma22
    out = list()
    out$cov = cov
    info=NULL
    if (!is.null(n)){
      info = n * F.entry(T[1]) * cov
      info[2, 2] = n * F.entry(T[2]) * cov[2,2]
    }
    out$info = info
    return(out)
   }
  }
  #UNDER H0
  ##Cov of n^(-1/2)*(U1, U2)
  cov.info.H0 = f.cov(f = f0, f.S = S0)
  cov.H0 = cov.info.H0$cov; info.H0 = cov.info.H0$info
  corr.H0 = 1
  if(is.matrix(cov.H0)) {corr.H0 = cov.H0[1,2]/sqrt(cov.H0[1,1]*cov.H0[2,2])}

  #UNDER Hp: Pooled Samples
  f.bar = function(t){r0 * f0(t) + r1 * f1(t)}
  S.bar = function(t){r0 * S0(t) + r1 * S1(t)}
  
  cov.info.Hp = f.cov(f = f.bar, f.S = S.bar)
  cov.Hp = cov.info.Hp$cov; info.Hp = cov.info.Hp$info
  corr.Hp = 1
  if(is.matrix(cov.Hp)) {corr.Hp = cov.Hp[1,2]/sqrt(cov.Hp[1,1]*cov.Hp[2,2])}
  
  o = list()
  o$cov.H0 = cov.H0
  o$info.H0 = info.H0
  o$corr.H0 = corr.H0
  o$cov.Hp = cov.Hp  
  o$info.Hp = info.Hp  
  o$corr.Hp = corr.Hp
  
  if(!is.null(f.ws[[1]])){wt1 = f.ws[[1]]} else{
    wt1 = data.frame(cbind(rho[1], gamma[1], tau[1], s.tau[1]))
  }
  if(!is.null(f.ws[[2]])){wt2 = f.ws[[2]]} else{
    wt2 = data.frame(cbind(rho[2], gamma[2], tau[2], s.tau[2]))
  }
  o$wt1 = wt1
  o$wt2 = wt2   
  return(o)
}
