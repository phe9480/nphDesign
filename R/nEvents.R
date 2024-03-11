#' Expected Number of Events At Calendar Time Since First Subject Randomized
#' 
#' This function calculates the expected number of events under alternative hypothesis
#' at calendar time t, which is calculated from first subject randomized. The function
#' returns the expected number of events for each arm at time t, based on the provided
#' enrollment distribution function and random lost-to-followup distribution if applicable.
#' If the total sample size is not provided, then only the corresponding probability of event
#' for each arm is provided.
#' 
#' @param T  Calendar time calculated from first subject randomization date.
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
#' @param n Total sample size for two arms. Default is NULL. 
#'
#' @return An object with dataframes below.
#'  \describe{
#'  \item{p.event}{Probability of event for each arm (p.event0: control group; p.event1: experimental group)}
#'  \item{n.events}{Expected number of events}
#'       \itemize{
#'       \item n.events0: number of events for control group
#'       \item n.events1: number of events for experimental group
#'       \item n.events.total: total number of events for two groups
#'       }
#'  \item{param}{Parameters specified: T, r, and n}
#'  \item{param.fun}{f0, f1, F.entry, G.ltfu}
#'  }
#'  
#' @examples 
#' #Example (1) Trial scenario: 1:1 randomization, n = 450, enrollment follows non-uniform 
#' enrollment distribution with weight 1.5 and enrollment period is 18 months. 
#' Control arm ~ exponential distribution with median 12 months, and 
#' Experimental arm ~ exponential distribution (Proportional Hazards) with median 12 / 0.7 months.
#' Assuming no lost-to-followup. Find the expected number of events at calendar time 24 months, i.e.
#' 6 months after last patient randomized.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' maxT = 30
#' nE = matrix(NA, nrow = maxT, ncol=3)
#' for (T in 1:maxT) {
#'   o = nEvents(T = T, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, n = 450)
#'   nE[T, 1] = o$n.events$n.events0; nE[T, 2] = o$n.events$n.events1; nE[T, 3] = o$n.events$n.events.total
#' }
#' plot(1:maxT, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:maxT, nE[, 1], lty = 1, col=1)
#' lines(1:maxT, nE[, 2], lty = 2, col=2)
#' lines(1:maxT, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' #Example (2) Same trial set up as example (1) but assuming delayed effect for 
#' experimental arm. The delayed period is assumed 6 months, and after delay the
#' hazard ratio is assumed 0.65.
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' h1 = function(t){lambda0*as.numeric(t < delay)+HR*lambda0*as.numeric(t >= delay)}
#' c = exp(-delay*lambda0*(1-HR)); 
#' S1 = function(t){exp(-lambda0*t)*as.numeric(t<delay) + c*exp(-HR*lambda0*t)*as.numeric(t>=delay)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' maxT = 30
#' nE = matrix(NA, nrow = maxT, ncol=3)
#' for (T in 1:maxT) {
#'   o = nEvents(T = T, r = 1, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, n = 450);
#'   nE[T, 1] = o$n.events$n.events0; nE[T, 2] = o$n.events$n.events1; nE[T, 3] = o$n.events$n.events.total
#' }
#' plot(1:maxT, nE[,3], type="n", xlab="Months", ylab = "Num of events")   
#' lines(1:maxT, nE[, 1], lty = 1, col=1)
#' lines(1:maxT, nE[, 2], lty = 2, col=2)
#' lines(1:maxT, nE[, 3], lty = 3, col=3)
#' legend(0, max(nE), c("Control", "Experimental", "Total"), col=1:3, lty=1:3, bty="n", cex=0.8)
#' 
#' @export

nEvents = function(T = 24, r = 1, h0 = function(t){log(2)/12}, S0=function(t){exp(-log(2)/12 * t)}, 
                   h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)},
          F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, G.ltfu = function(t){0}, n = 450){

  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #Density function
  f0 = function(t) {return(h0(t) * S0(t))}
  f1 = function(t) {return(h1(t) * S1(t))}
  
  #Integrand of control arm and experimental arm for calculation of prob. of event
  I0 = function(t){F.entry(T-t) * (1 - G.ltfu(t)) * f0(t)}
  I1 = function(t){F.entry(T-t) * (1 - G.ltfu(t)) * f1(t)}
  
  #prob. of event for control and experimental arm
  p.event0 = integrate(I0, lower=0, upper=T)$value
  p.event1 = integrate(I1, lower=0, upper=T)$value

  p.event = data.frame(cbind(p.event0, p.event1))
  
  #expected number of events
  if(!is.null(n)){
    n.events0 = n * F.entry(T) * r0 * p.event0
    n.events1 = n * F.entry(T) * r1 * p.event1
    n.events.total = n.events0 + n.events1
  }
  
  n.events = data.frame(cbind(n.events0, n.events1, n.events.total))
  
  param = data.frame(cbind(T, r, n))
  param.fun = list()
  param.fun$f0 = f0
  param.fun$f1 = f1
  param.fun$F.entry = F.entry
  param.fun$G.ltfu = G.ltfu
  
  o = list()
  o$p.event = p.event
  o$n.events = n.events
  o$param = param
  o$param.fun = param.fun
  return(o)
}
