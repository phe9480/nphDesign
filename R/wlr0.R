#' Unstratified Weighted Log-rank Test
#' 
#' The weight function includes Stabilized Fleming-Harrington class (rho, gamma, tau, s.tau), 
#' and any user-defined weight function. For the stabilized Fleming-Harrington class,
#' to produce stabilized weighting function from pooled survival curve, specify either 
#' tau or s.tau, which are thresholds in survival time and survival rate, respectively. 
#'
#' 
#' @param  time  Survival time
#' @param  event Event indicator; 1 = event, 0 = censor
#' @param  group Treatment group; 1 = experimental group, 0 = control
#' @param  rho Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#' @param  gamma Parameter for Fleming-Harrington (rho, gamma) weighted log-rank test.
#'         For log-rank test, set rho = gamma = 0.
#' @param  tau  Cut point for stabilized FH test, sFH(rho, gamma, tau); with weight
#'       function defined as w(t) = s_tilda^rho*(1-s_tilda)^gamma, where
#'       s_tilda = max(s(t-), s.tau) or max(s(t-), s(tau)) if s.tau = NULL
#'       tau = Inf reduces to regular Fleming-Harrington test(rho, gamma). 
#' @param  s.tau  Survival rate cut S(tau) at t = tau; default 0.5, ie. cut at median.
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
#' wlr0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = NULL, s.tau=0.5, f.ws=function(s){1/max(s, 0.25)})
#' wlr0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = 10, s.tau=NULL)
#' wlr0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = 3, s.tau=NULL)
#' wlr0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = 20, s.tau=NULL)
#' wlr0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = Inf, s.tau=NULL)
#' wlr0(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=1, tau = 10, s.tau=NULL)
#' wlr0(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=1, tau = 10, s.tau=0.5, side="one.sided")
#' wlr0(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=0, tau = 10, s.tau=0.5, side="one.sided")
#' wlr0(time=c(12,7,10,5,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(1,1,0,0,0,1,0,1), rho=0, gamma=0, tau = 10, s.tau=0, side="one.sided")
#' ######################################
#' #Validation Note:
#' #Fully validated for unstratified version using nph::logrank.test() function.
#' ex1 = nph::logrank.test(time=time,   event=event,  group=group,  rho = 0,  gamma = 1.5)
#' time=rexp(100); event=sample(c(0,1), 100, replace = TRUE); group=c(rep(0, 50), rep(1, 50)); 
#' rho=0; gamma=1; tau = NULL; s.tau=0; strata1=sample(c(1,2), 100, replace = TRUE)
#' wlr0(time=time,   event=event,  group=group,  rho = 0,  gamma = 1.5, tau = NULL, s.tau=0)
#' wlr0(time=time,   event=event,  group=group,  rho = 0,  gamma = 1.5, tau = NULL, s.tau=NULL, f.ws=function(s){1/max(s, 0.25)})
#' 
#' @export
#' 
#' 
wlr0 = function(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1),
                group=c(0,1,0,1,0,1,0,1), rho=0, gamma=1, tau = NULL, s.tau=0.5, 
                f.ws=NULL, side = c("two.sided", "one.sided")) {
  
  u.eTime = unique(time[event==1]) #vector of unique event times
  u.Ne = length(u.eTime) #number of unique event times
  if (u.Ne == 0) {return(NULL); stop("No events")}
  
  u.eTime = sort(u.eTime) #sort by increasing order of unique event times
  
  Y0 = Y1 = Y = dN = dN0 = dN1 = rep(NA, u.Ne) #risk-set, death-set
  
  for (i in 1:u.Ne) {
    Y0[i] = sum(time>=u.eTime[i] & group == 0)
    Y1[i] = sum(time>=u.eTime[i] & group == 1)
    dN0[i] = sum(time == u.eTime[i] & group == 0 & event == 1)
    dN1[i] = sum(time == u.eTime[i] & group == 1 & event == 1)
  }
  
  Y=Y0+Y1
  dN = dN0 + dN1
  
  s = rep(NA, u.Ne)
  
  s[1] = 1 - dN[1]/Y[1]
  if (u.Ne > 1) {
    for (i in 2:u.Ne) {
      s[i] = s[i-1] * (1 - dN[i]/Y[i])
    }
  }
  w = rep(NA, u.Ne)
  #If f.ws() function is provided, then use directly. 
  if(!is.null(f.ws)){
    w[1] = f.ws(1)
    for(i in 2:u.Ne){
      w[i] = f.ws(s[i-1])
      s.til=NULL
    }
  } else {
    #Find S(tau): = S(max(ti)|ti <= tau), where ti is unique event time.
    if (is.null(s.tau)) {
      if (is.null(tau)){stop("tau or s.tau or f.ws is required for weight function.")}else{
        s.tau = 1
        for (i in 1:u.Ne) {
          if (u.eTime[i] <= tau) {s.tau = s[i]} else {break}
        }
      }
    } #if s.tau is missing, find s.tau by tau.
    s.til = apply(cbind(s, s.tau),MARGIN=1,FUN=max)

    #w = s.til^rho*(1-s.til)^gamma
    #Special handling to ensure FH00 = logrank; FH10 = wilcoxon
    if(gamma == 0){w[1] = 1} else{w[1] = 0}
    for (i in 2:u.Ne){w[i] = s.til[i-1]^rho*(1-s.til[i-1])^gamma}
  }
  
  U=w*(dN0-Y0/Y*dN)
  V=w^2*Y0*Y1/(Y^2)*(Y-dN)/(Y-1)*dN
  
  z=sum(U[!is.nan(V)])/sqrt(sum(V[!is.nan(V)]))
  
  if(side[1] == "one.sided") {
    test.side = 1; p = 1-pnorm(z)
  } else {
      test.side = 2; p = 2*(1-pnorm(abs(z)))
  }
  
  #create a dataframe to output the parameters
  chisq = z*z
  test.results = data.frame(cbind(z, chisq, p, test.side))
  
  #create a dataframe to output the list of unique event times and statistics
  uni.event.time = data.frame(cbind(u.eTime,Y0, Y1, Y, dN0, dN1, dN, s, s.til, w, U, V, z))
  
  #create a dataframe to output the original data
  data = data.frame(cbind(time, event, group))
  
  o=list()
  o$uni.event.time = uni.event.time
  #o$data = data
  o$test.results = test.results
  if(!is.null(f.ws)){wt = f.ws} else{
    wt = data.frame(cbind(rho, gamma, tau, s.tau))
  }
  o$wt = wt
  return(o)
}
