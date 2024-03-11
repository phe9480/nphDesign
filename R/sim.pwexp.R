#' Simulate data following piece-wise exponential distribution with non-uniform accrual pattern and lost to follow-up
#'
#' Simulate Randomized two-arm trial data with the following characteristics:  
#' (1) randomization time (entry time) is generated according to the specified non-uniform accrual pattern, 
#' i.e. the cumulative recruitment at calendar time t is (t/A)^w with weight w and enrollment complete in A months.
#' w = 1 means uniform enrollment, which is usually not realistic due to graduate sites activation process.
#' (2) Survival time follows piece-wise exponential distribution for each arm.
#' (3) N total patients with r:1 randomization ratio
#' (4) Random drop off can be incorporated into the censoring process.
#' (5) Data cutoff dates are determined by specified vector of target events for all analyses.
#' (6) A dataset is generated for each analysis according to the specified number of target events. 
#'     Multiple analyses can be specified according to the vector of targetEvents, eg, targetEvents = c(100, 200, 300) 
#'     defines 3 analyses at 100, 200, and 300 events separately.
#'
#' @param  nSim Number of trials
#' @param  N Total number patients in two arms.
#' @param  A Total accrual period in months
#' @param  w Weight parameter in cumulative enrollment pattern. 
#' The cumulative enrollment at month t is (t / A)^w, eg, at month 6, 
#'   the enrollment is N*(6/24)^2 = N/16 for 24 months planned accrual period.
#' @param  r Randomization ratio r:1, where r refers to the experimental arm, eg, r=2 in 2:1 ratio
#' @param  lambda0 Hazard rates for control arm of intervals defined by cuts; for exponential(lambda0) distribution,
#'         lambda0 = log(2) / median;
#' @param  lambda1 Hazard rates for experimental arm for intervals; for exponential(lambda1) distribution,
#'         lambda1 = log(2) / median; For delayed effect under H1, lambda1 is a vector (below).
#' @param  cuts Timepoints to form intervals for piecewise exponential distribution. For example,
#   \itemize{
#   \item Proportional hazards with hr = 0.65. Then lambda0 = log(2)/m0, lambda1 = log(2)/m0*hr, cuts = NULL. 
#   \item Delayed effect at month 6, and control arm has constant hazard (median m0) and 
#'       experimental arm has hr = 0.6 after delay, then cuts = 6, and 
#'       lamda0 = log(2) / m0 or lambda0 = rep(log(2) / m0, 2), 
#'       lamda1 = c(log(2)/m0, log(2)/m0*hr). 
#   \item Delayed effect at month 6, and control arm has crossover to subsequent IO 
#'       treatment after 24 mo, so its hazard decreases 20%. Then, 
#'       lambda0 = c(log(2)/m0, log(2)/m0, log(2)/m0*0.8), 
#'       lambda1 = c(log(2)/m0, log(2)/m0*hr, log(2)/m0*hr), and
#'       cuts = c(6, 24), which forms 3 intervals (0, 6), (6, 24), (24, infinity)
#       }
#' @param dropOff0 Drop Off rate per month, eg, 1%, for control arm
#' @param dropOff1 Drop Off rate per month, eg, 1%, for experimental arm
#' @param targetEvents A vector of target events is used to determine DCOs. For example, 
#'              397 target events are used to determine IA DCO; and 496 events are used 
#'              to determine the FA cutoff.
#' @return An object with a dataframe for each analysis including the following variables:
#' \describe{
#' \item{sim}{sequence number of simulated dataset;}
#' \item{treatment}{treatment group with values of "control" and "experimental"}
#' \item{enterTime}{Time of randomization in calendar time}
#' \item{calendarTime}{the time when event/censoring occurred in calendar time}
#' \item{survTime}{Survival time for analysis, = calendarTime - enterTime}
#' \item{cnsr}{censor status (0=event; 1=censor) before administrative censoring due to data cut}
#' \item{calendarCutOff}{Data CutOff Time (DCO);}
#' \item{survTimeCut}{Survival time after cut}
#' \item{cnsrCut}{Censor status after cut}
#' }
#' @examples
#' #Example (1): Simulate 10 samples from proportional hazards scenario. 
#' #Total 672 pts, 1:1 randomization, control median OS 11.7 mo; 
#' #HR = 0.65, enrollment 21 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 397 and 496 events respectively.
#' 
#' N=672; m0 = 11.7; A=21; r=1; hr = 0.65; w = 1.5; dropOff0 = dropOff1 = 0; 
#' targetEvents = c(397, 496); lambda0 = log(2) / m0; lambda1 = log(2)/m0*hr; cuts = NULL;
#' sim.ph = sim.pwexp(nSim=10, N = N, A = A, w=w, r=r, lambda0=lambda0, lambda1= lambda1, cuts=cuts, dropOff0= dropOff0, dropOff1= dropOff1, targetEvents = targetEvents)
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.ph[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Example (2): Simulate 10 samples with delayed effect at month 6;
#' #control arm has constant hazard (median 11.7 mo) and experimental arm has hr = 0.65 after delay.
#' #Total 672 pts, 1:1 randomization, control median OS 11.7 mo. 
#' #HR = 0.65, enrollment 21 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 397 and 496 events respectively.
#' 
#' N=672; m0 = 11.7; A=21; r=1; hr = 0.65; w = 1.5; dropOff0 = dropOff1 = 0; 
#' targetEvents = c(397, 496); cuts = 6
#' lambda0 = log(2) / m0; lambda1 = c(log(2)/m0, log(2)/m0*hr)
#' 
#' sim.delay6 = sim.pwexp(nSim=1, N = N, A = A, w=w, r=r, lambda0=lambda0, lambda1=lambda1, cuts=cuts, dropOff0= dropOff0, dropOff1= dropOff1, targetEvents = targetEvents)
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Example (3): Simulate 10 samples with delayed effect at month 6 
#' #Control arm has crossover to subsequent IO after 24 mo, so its hazard decreases 20%.
#' #control arm has constant hazard (median 11.7 mo) and experimental arm has 
#' #hr = 1 and 0.65 at intervals (0, 6) and (6, 24) respectively.
#' #Total 672 pts, 1:1 randomization, control median OS 11.7 mo; 
#' #HR = 0.65, enrollment 21 months, weight 1.5, no drop offs; 
#' #IA and FA are performed at 397 and 496 events respectively.
#' 
#' N=672; m0 = 11.7; A=21; r=1; hr = 0.65; w = 1.5; dropOff0 = dropOff1 = 0; 
#' targetEvents = c(397, 496); cuts = c(6, 24); crossTime = 24; crossEffect = 0.8
#' lambda0 = log(2)/m0*c(1, 1, crossEffect); lambda1 = log(2)/m0*c(1, hr, hr)
#' 
#' sim.delay6crs=sim.pwexp(nSim=10,N=N,A=A,w=w,r=r,lambda0=lambda0, lambda1=lambda1,cuts=cuts,dropOff0=dropOff0,dropOff1=dropOff1, targetEvents=targetEvents)
#' km.IA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=sim.delay6crs[[1]][sim==1,])
#' plot(km.IA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data= sim.delay6crs[[2]][sim==1,])
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' @export 
sim.pwexp = function(nSim=100, N = 672, A = 21, w=1.5, r=1, lambda0=log(2)/11.7, lambda1=log(2)/11.7*0.745, cuts=NULL, dropOff0=0, dropOff1=0, targetEvents = c(397, 496)) {
  
  nEachMonth = f.nEachMonth(N=N, A=A, w=w, r=r)
  
  gamma = nEachMonth$n0 + nEachMonth$n1
  eta0 = -log(1-dropOff0)
  eta1 = -log(1-dropOff1)
  
  o = nphsim::nphsim(nsim=nSim,lambdaC=lambda0,lambdaE=lambda1, ssC=N/(r+1), ssE=N*r/(r+1),
             intervals=cuts, gamma=gamma, R=rep(1, A),eta=eta0, etaE=eta1,fixEnrollTime = FALSE)
  dat = o$simd
  data.out = NULL
  L = length(targetEvents) #number of analyses
  
  for (j in 1:nSim) {
    
    dataj = dat[dat$sim == j,]
    
    #rename variables
    dataj <- dplyr::rename(dataj, enterTime = enterT, calendarTime = ct, survTime= survival)

    #cut data according to the specified target events    
    dataj.cut <- lapply(1:L, function(k) {
      f.dataCut(data=dataj, targetEvents[k])
    })
    
    data.out=lapply(1:L, function(k) {
      rbind(data.out[[k]], dataj.cut[[k]])
    })
    
  }
  return(data.out)
}