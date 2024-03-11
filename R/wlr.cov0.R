#' Covariance Matrix of Two Weighted Log-rank Tests (Unstratified) At The Same Analysis
#'
#' This function calculates the covariance matrix of two weighted log-rank tests at the same analysis time.
#' The two weight functions are specified by stabilized Fleming-Harrington class
#' with parameters (rho, gamma, tau, s.tau), where tau and s.tau are thresholds for
#' survival time and survival rates, respectively. Either tau or s.tau can be specified.
#' tau = Inf or s.tau = 0 reduces to the Fleming-Harrington test (rho, gamma).
#' User-defined weight functions f.ws1 and f.ws2 can be used as well. For example,
#' f.ws1 = function(s){s^rho*(1-s)^gamma} is equivalent to Fleming-Harrington (rho, gamma)
#' test with parameters specified for rho and gamma.
#' 
#' @param  time  Survival time
#' @param  event Event indicator; 1 = event, 0 = censor
#' @param  group Treatment group; 1 = experimental group, 0 = control
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
#'         For example, f.ws2 = function(s){1/max(s, 0.25)}
#'         When f.ws1 or f.ws2 is specified, the weight function takes them as priority.
#'         
#' @return An object with dataframes below.
#' \describe{
#' \item{data}{dataframe with variables: time, event, and group.}
#' \item{uni.event.time}{dataframe with variables}
#'          \itemize{
#'          \item u.eTime: Unique event times;
#'          \item Y0:      Risk set of control arm at each of unique event times;
#'          \item Y1:      Risk set of experimental arm at each of unique event times;
#'          \item Y:       Risk set of pooled data at each of unique event times;
#'          \item dN0:     Event set of control arm at each of unique event times;
#'          \item dN1:     Event set of experimental arm at each of unique event times;
#'          \item dN:      Event set of pooled data at each of unique event times;
#'          \item s:       Survival time of pooled data by KM method;
#'          \item s.til1:  s.tilda1 defined as, s.til = max(s, s.tau1);
#'          \item s.til2:  s.tilda2 defined as, s.til = max(s, s.tau2);
#'          \item w1:      Weight function w1(t) at each of unique event times;
#'          \item w2:      Weight function w2(t) at each of unique event times;
#'          \item V11:     Variance statistic at each of unique event times for w1;
#'          \item V12:     Covariance statistic at each of unique event times for w1 and w2;
#'          \item V22:     Variance statistic at each of unique event times for w2;
#'          }
#' \item{corr}{Correlation between two weigthed log-rank score statistics U1 and U2, equivalent to the correlation between two normalized weigthed log-rank statistics Z1 and Z2, where Zi = Ui/sqrt(var(Ui))}
#' \item{cov}{Covariance between two weighted log-rank score statistics, U1 and U2.}
#' }
#' @examples
#' wlr.cov0(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1), group=c(0,1,0,1,0,1,0,1), rho1=0, gamma1=0, tau1 = NULL, s.tau1=0.5,rho2=0, gamma2=1, tau2 = NULL, s.tau2=0.5, f.ws1=NULL, f.ws2=NULL)
#' @export  
#' 
wlr.cov0 = function(time=c(5,7,10,12,12,15,20,20), event=c(1,0,0,1,1,0,1,1),
                    group=c(0,1,0,1,0,1,0,1), rho1=0, gamma1=0, tau1 = NULL, s.tau1=0.5,
                    rho2=0, gamma2=1, tau2 = NULL, s.tau2=0.5,
                    f.ws1=NULL, f.ws2=NULL) {

  u.eTime = unique(time[event==1]) #vector of unique event times
  u.Ne = length(u.eTime) #number of unique event times
  if (u.Ne == 0) {return(NULL); stop("No events")}
  
  u.eTime = sort(u.eTime) #sort by increasing order of unique event times
  
  Y0 = Y1 = Y = dN = dN0 = dN1 = rep(NA, u.Ne) #risk-set, death-set
  
  for (i in 1:u.Ne) {
    Y0[i] = sum(time>=u.eTime[i] & group == 0)
    Y1[i] = sum(time>=u.eTime[i] & group == 1)
    dN0[i] = sum(time == u.eTime[i] & group == 0)
    dN1[i] = sum(time == u.eTime[i] & group == 1)
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
  
  f.w = function(s=s, f.ws=f.ws1, tau=tau1, s.tau=s.tau1, rho=rho1, gamma=gamma1) {
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
      
      if(gamma == 0){w[1] = 1} else{w[1] = 0}; 
      for (i in 2:u.Ne){ w[i] = s.til[i-1]^rho*(1-s.til[i-1])^gamma }
      #w = s.til^rho*(1-s.til)^gamma
    }
    ot = list()
    ot$w = w; ot$s.til = s.til
    return(ot)
  }
  ow1 = f.w(s=s, f.ws=f.ws1, tau=tau1, s.tau=s.tau1, rho=rho1, gamma=gamma1)
  ow2 = f.w(s=s, f.ws=f.ws2, tau=tau2, s.tau=s.tau2, rho=rho2, gamma=gamma2)
  w1 = ow1$w; s.til1 = ow1$s.til
  w2 = ow2$w; s.til2 = ow2$s.til
  
  V = matrix(NA, nrow=2, ncol=2)  
  v11 = w1*w1*Y0*Y1/(Y^2)*(Y-dN)/(Y-1)*dN
  v12 = w1*w2*Y0*Y1/(Y^2)*(Y-dN)/(Y-1)*dN
  v22 = w2*w2*Y0*Y1/(Y^2)*(Y-dN)/(Y-1)*dN
  
  V[1,1]= sum(v11[!is.nan(v11)])
  V[1,2]=V[2,1]=sum(v12[!is.nan(v12)])
  V[2,2]=sum(v22[!is.nan(v22)])
  
  corr = V[1,2]/(sqrt(V[1,1]*V[2,2]))
  
  #create a dataframe to output the list of unique event times and statistics
  uni.event.time = data.frame(cbind(u.eTime,Y0, Y1, Y, dN0, dN1, dN, s, s.til1, s.til2, w1, w2, v11, v12, v22))
  
  #create a dataframe to output the original data
  data = data.frame(cbind(time, event, group))
  
  o=list()
  o$uni.event.time = uni.event.time
  #o$data = data
  o$corr = corr
  o$cov = V
  
  if(!is.null(f.ws1)){wt1 = f.ws1} else{
    wt1 = data.frame(cbind(rho1, gamma1, tau1, s.tau1))
  }
  if(!is.null(f.ws2)){wt2 = f.ws2} else{
    wt2 = data.frame(cbind(rho2, gamma2, tau2, s.tau2))
  }
  o$wt1 = wt1
  o$wt2 = wt2   
  return(o)
}

