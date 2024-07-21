#' Stratified Weighted Log-rank Test
#' 
#' The weight function includes Stabilized Fleming-Harrington class (rho, gamma, tau, s.tau),
#' and any user-defined weight function. For the stabilized Fleming-Harrington class,
#' to produce stabilized weighting function from pooled survival curve, specify either 
#' tau or s.tau, which are thresholds in survival time and survival rate, respectively. 
#' The weight function is based on the pooled survival curve within each strata by default.
#' 
#' @param  time  Survival time
#' @param  event Event indicator; 1 = event, 0 = censor
#' @param  group Treatment group; 1 = experimental group, 0 = control
#' @param  rho Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#' @param  gamma Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t), s.tau) or max(s(t), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  s.tau  Survival rate cut S(tau) at t = tau; default 0.
#'       s.tau = 0 reduces to regular Fleming-Harrington test(rho, gamma)
#' @param  f.ws  Self-defined weight function of survival rate. 
#'         For example, f.ws = function(s){1/max(s, 0.25)}
#'         When f.ws is specified, the weight function takes it as priority.
#'         
#'         Note: The actual math formula for weighting function is based on 
#'         the left of each event time t, i.e., w(t) = f.ws(s(t-)). 
#'         For FH(rho, gamma) test, when gamma = 0, then the first event time has 
#'         weight 1;  when gamma > 0, the first event weight is 0. This can 
#'         ensure consistency with FH(0,0) = logrank; and FH(1,0) = generalized Wilcoxon.
#'         
#' @param  strata1 Stratification variable 1
#' @param  strata2 Stratification variable 2
#' @param  strata3 Stratification variable 3
#' @param  side  Type of test. one.sided or two.sided. default = two.sided
#' 
#' @return An object with dataframes below.
#' \describe{
#' \item{data}{dataframe including the following variables:}
#'          \itemize{
#'          \item time: Survival time 
#'          \item event: Event indicator; 0 = censor; 1 = event
#'          \item group: Group indicator; 0 = control; 1 = experimental treatment
#'          }
#' \item{uni.event.time}{dataframe of unique event times and intermediate statistics
#'          including the following variables:}
#'          \itemize{
#'          \item u.eTime: Unique event times;
#'          \item Y0:      Risk set of control arm at each of unique event times;
#'          \item Y1:      Risk set of experimental arm at each of unique event times;
#'          \item Y:       Risk set of pooled data at each of unique event times;
#'          \item dN0:     Event set of control arm at each of unique event times;
#'          \item dN1:     Event set of experimental arm at each of unique event times;
#'          \item dN:      Event set of pooled data at each of unique event times;
#'          \item s:       Survival time of pooled data by KM method;
#'          \item s.til:   s.tilda defined as, s.til = max(s, s.tau);
#'          \item w:       Weight function w(t) at each of unique event times;
#'          \item U:       Score statistic at each of unique event times;
#'          \item V:       Variance statistic at each of unique event times;
#'          \item z:       Normalized z-statistic at each of unique event times; z = sum(U)/sqrt(sum(V));
#'          }
#' \item{test.results.strata}{dataframe of intermediate calculation per stratum}
#'          \itemize{
#'          \item rho:     Stabilized Fleming-Harrington test parameter
#'          \item gamma:   Stabilized Fleming-Harrington test parameter
#'          \item tau:     Stabilized Fleming-Harrington test parameter if provided
#'          \item s.tau:   Stabilized Fleming-Harrington test parameter if provided
#'          \item test.side: two.sided or one.sided
#'          \item chisq:   Chi-square statistic, chisq = z^2
#'          \item z:       Z statistic
#'          \item p:       P value per stratum
#'          \item strata1 Strata 1 value
#'          \item strata1 Strata 2 value
#'          \item strata1 Strata 3 value
#'          }
#' \item{test.results}{dataframe including the following variables:}
#'          \itemize{
#'          \item rho:     Stabilized Fleming-Harrington test parameter
#'          \item gamma:   Stabilized Fleming-Harrington test parameter
#'          \item tau:     Stabilized Fleming-Harrington test parameter if provided
#'          \item s.tau:   Stabilized Fleming-Harrington test parameter if provided
#'          \item test.side: two.sided or one.sided
#'          \item chisq:   Chi-square statistic, chisq = z^2
#'          \item z:       Z statistic
#'          \item p:       P value
#'          }
#' }
#' @examples
#' wlr(time=rexp(100), event=sample(c(0,1), 100, replace = TRUE), group=c(rep(0, 50), rep(1, 50)), rho=0, gamma=1, tau = NULL, s.tau=0, strata1=sample(c(1,2), 100, replace = TRUE),strata2=NULL, strata3=NULL)
#' wlr(time=rexp(100), event=sample(c(0,1), 100, replace = TRUE), group=c(rep(0, 50), rep(1, 50)), rho=0, gamma=1, tau = NULL, s.tau=0, strata1=sample(c(1,2), 100, replace = TRUE), strata2=sample(c(1,2), 100, replace = TRUE))
#' wlr(time=rexp(100), event=sample(c(0,1), 100, replace = TRUE), group=c(rep(0, 50), rep(1, 50)), rho=0, gamma=1, tau = NULL, s.tau=0, strata1=sample(c(1,2), 100, replace = TRUE), strata2=sample(c(1,2), 100, replace = TRUE), strata3=sample(c(1,2), 100, replace = TRUE))
#' wlr(time=rexp(100), event=sample(c(0,1), 100, replace = TRUE), group=c(rep(0, 50), rep(1, 50)), rho=0, gamma=1, tau = NULL, s.tau=0, strata1=sample(c(1,2), 100, replace = TRUE), strata2=sample(c(1,2), 100, replace = TRUE), strata3=sample(c(3,4), 100, replace = TRUE), f.ws=function(s){1/max(s^2, 0.25)})
#' wlr(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = 3, s.tau=NULL)
#' wlr(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = 20, s.tau=NULL)
#' wlr(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = Inf, s.tau=NULL)
#' wlr(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=1, tau = 10, s.tau=NULL)
#' wlr(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=1, tau = 10, s.tau=0.5, side="one.sided")
#' wlr(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=0, tau = 10, s.tau=0.5, side="one.sided")
#' wlr(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=0, tau = 10, s.tau=0, side="one.sided")
#' t = rexp(100); e = sample(c(0,1), 100, replace = TRUE); g = c(rep(0, 50), rep(1, 50)); str1 = sample(c(1,2), 100, replace = TRUE)
#' #FH(0, 1); s.tau = 0 means no threshold, so reduces to FH(0, 1)
#' wlr(time=t, event=e, group=g, rho=0, gamma=1, tau = NULL, s.tau=0, strata1=str1)
#' 
#' #FH(0, 1); when tau = Inf, so reduces to FH(0, 1)
#' wlr(time=t, event=e, group=g, rho=0, gamma=1, tau = Inf, s.tau=NULL, strata1=str1)
#' 
#' #FH(0, 1, tau=2)
#' wlr(time=t, event=e, group=g, rho=0, gamma=1, tau = 2, s.tau=NULL, strata1=str1)
#' 
#' #FH(0, 1, s.tau=0.5)
#' wlr(time=t, event=e, group=g, rho=0, gamma=1, tau = NULL, s.tau=0.5, strata1=str1)
#' 
#' #FH(0, 1, s.tau=0.5); tau is ignored if s.tau is available.
#' wlr(time=t, event=e, group=g, rho=0, gamma=1, tau = 2, s.tau=0.5, strata1=str1)
#' 
#' #logrank test; s.tau set as 0 but it is actually not used in calculation.
#' wlr(time=t, event=e, group=g, rho=0, gamma=0, tau = NULL, s.tau=0, strata1=str1)
#' 
#' #logrank test; s.tau value doesn't make any difference, because it is actually not used in calculation.
#' wlr(time=t, event=e, group=g, rho=0, gamma=0, tau = NULL, s.tau=0.5, strata1=str1)
#' @export
#' 
wlr = function(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1),
               group=c(0,1,0,1,0,1,0,1), strata1=NULL, strata2=NULL, strata3=NULL, 
               rho=0, gamma=1, tau = NULL, s.tau=0,
               f.ws=NULL, side = c("two.sided", "one.sided")) {
  #Initialize strata1-3
  n = length(time)
  if(is.null(strata1)){strata1 = rep(1,n)}
  if(is.null(strata2)){strata2 = rep(1,n)}
  if(is.null(strata3)){strata3 = rep(1,n)}
  u.strata1 = unique(strata1)
  u.strata2 = unique(strata2)
  u.strata3 = unique(strata3)
  
  n.strata1 = length(u.strata1)
  n.strata2 = length(u.strata2)
  n.strata3 = length(u.strata3)
  n.strata = n.strata1*n.strata2*n.strata3
  
  #unstratified analysis
  if (n.strata == 1) {
    o = wlr0(time=time, event=event, group=group, rho=rho, gamma=gamma, tau = tau, 
             s.tau=s.tau, f.ws=f.ws, side = side)
  }else{
    uni.event.time=data=test.results.strata=NULL
    for (i in 1:n.strata1) {
      for (j in 1:n.strata2){
        for (k in 1:n.strata3){
          ix = (strata1==u.strata1[i] & strata2 == u.strata2[j] & strata3 == u.strata3[k])
          #ignore empty strata
          if(sum(ix) > 0){
            timei = time[ix]
            eventi = event[ix]
            groupi = group[ix]
           
            wlr.i = wlr0(time=timei, event=eventi, group=groupi, rho=rho, gamma=gamma, tau = tau, 
                     s.tau=s.tau, f.ws=f.ws, side = side)
          
            wlr.i$uni.event.time$strata1 = u.strata1[i]
            wlr.i$uni.event.time$strata2 = u.strata2[j]
            wlr.i$uni.event.time$strata3 = u.strata3[k]
            uni.event.time = rbind(uni.event.time, wlr.i$uni.event.time)
          
            wlr.i$data$strata1 = u.strata1[i];
            wlr.i$data$strata2 = u.strata2[j]
            wlr.i$data$strata3 = u.strata3[k]
            data = rbind(data, wlr.i$data);
          
            wlr.i$test.results$strata1 = u.strata1[i]
            wlr.i$test.results$strata2 = u.strata2[j]
            wlr.i$test.results$strata3 = u.strata3[k]
            test.results.strata = rbind(test.results.strata, wlr.i$test.results)
          }  
        } 
      }       
    }  
    
    U=uni.event.time$U
    V=uni.event.time$V
    z=sum(U[!is.nan(V)])/sqrt(sum(V[!is.nan(V)]))
    
    if(side[1] == "one.sided") {
      test.side = 1; p = 1-pnorm(z)
      } else {
        test.side = 2
        p = 2*(1-pnorm(abs(z)))
      }
    
    #create a dataframe to output the parameters
    #chisq = z*z; 
    test.results = data.frame(cbind(z, p, test.side))
    
    o=list()
    
    o$uni.event.time = uni.event.time
    #o$data = data
    
    o$test.results.strata = test.results.strata
    o$test.results = test.results

    if(!is.null(f.ws)){wt = f.ws} else{
      wt = data.frame(cbind(rho, gamma, tau, s.tau))
    }
    o$wt = wt
  }
  return(o)
}
