#' Display of Cumulative Number of Events Over Time per Study Design
#' 
#' This function plots the cumulative number of events for each arm and total per study design
#' 
#' @param t  A sequence of calendar time for the plot (x-axis)
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
#' @param ... Other graphic parameters passed to the plot
#'
#' @return Display of the graph
#'  
#' @examples 
#' 
#' HR = 0.65; delay = 6; lambda0 = log(2) / 12; 
#' h0 = function(t){lambda0}; S0 = function(t){exp(-lambda0 * t)}
#' Hazard function and survival function for experimental arm
#' lambda1 = lambda0 * HR
#' h1 = function(t){lambda1}; S1= function(t){exp(-lambda1 * t)}
#' F.entry = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}
#' 
#' leg = list(x=0, y=80, txt=c("Control", "Exp. Arm", "Total"))
#' par = list(main="OS: Number of Events Over time",
#'            xlab="Calendar Time (mo)", 
#'            ylab="Events")
#'            
#' plot_events(Tmax = 50, r=1, h0=h0, S0=S0, h1=h1, S1=S1, 
#' F.entry = F.entry, G.ltfu = G.ltfu, n = 100, leg = leg, par = par)
#' 
#' @export
#' 
plot_events = function(Tmax = 50, r=1, h0=h0, S0=S0, h1=h1, S1=S1, 
                F.entry = F.entry, G.ltfu = G.ltfu, n = 100, 
                leg=list(x=0, y=n/2,txt=c("Control", "Exp. Arm", "Total")), 
                param=list(xlab = "Calendar Time (mo)", 
                         ylab = "Cumulative Events",
                         main = "Cumulative Events")){
  t = seq(1, Tmax, by = 1)
  nE = matrix(NA, nrow = length(t), ncol=3)
  for (i in 1:length(t)) {
    o = f.nEvents(T = t[i], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, 
                  F.entry = F.entry, G.ltfu = G.ltfu, n = n);
    nE[i, 1] = o$n.events$n.events0; 
    nE[i, 2] = o$n.events$n.events1; 
    nE[i, 3] = o$n.events$n.events.total
  }

  plot(t, nE[, 3], type = "n", xlab=param$xlab, ylab=param$ylab, main=param$main) 
  lines(t, nE[, 1], lty = 1, col=1, lwd=3)
  lines(t, nE[, 2], lty = 2, col=2, lwd=3)
  lines(t, nE[, 3], lty = 3, col=3, lwd=3)
  abline(h = seq(0, max(nE), 25), col="gray80", lty=3)
  abline(v=seq(0, Tmax, by=2), col="gray80", lty=3)
  
  legend(0, max(nE[, 3]), leg$txt, col=1:3, lty=1:3, bty="n", cex=0.8)
  
}

