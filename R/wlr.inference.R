#'  Calculation of Actual Rejection Boundary and Statistical Inference for Group Sequential Design Using Weighted Log-rank Tests
#' 
#'  This function calculates the actual rejection boundary for group sequential design 
#'  based on the asymptotic 
#'  distribution of the weighted log-rank test statistic under H1 for the next 
#'  analysis based on the datasets from previous analyses. In standard log-rank test,
#'  the rejection boundaries can be pre-determined according to the alpha spending function
#'  for the information fraction. The boundaries can be pre-determined because 
#'  of the good approximation of information fraction based on total number of events. 
#'  However, using weighted log-rank test, the weight function is usually 
#'  associated with the pooled survival curve from the actual data and the weight
#'  function can significantly impact the Fisher's information. As a result, the
#'  limiting weight function of pooled survival curve based on study design assumptions
#'  often can not provide accurate approximation of actual pooled survival curve using
#'  Kaplan-Meier method. The actual correlation between test statistics at different 
#'  analyses needs to be based on actual data under asymptotic distribution 
#'  derived by Tsiatis (1982).
#'  
#'  
#' @param datasets  The datasets used for the analyses. For the first analysis,
#'                  just provide in the format of list(IA1 = data1). For the ith 
#'                  analysis, needs to provide all datasets up to the ith 
#'                  analysis in the format of list(IA1=data1, IA2=data2, ...,IAi = datai).
#'                  All datasets (data1, data2, ..., datai) must be sorted by
#'                  subject id and having a format of 1 subject 1 record. In addition,
#'                  the same number of subjects is required in order to produce
#'                  consistent results.
#'                  
#'                  Each dataset must include the following variables:
#'                  \itemize{
#'                  \item survTimeCut  Survival time for each analysis
#'                  \item cnsrCut      Censoring status for each analysis (0=event; 1=censored)
#'                  \item group        Treatment group (0 = control, 1 = experimental treatment)
#'                  \item strata1, strata2, strata3: Stratification variables. 
#'                  Only 3 stratification factors are allowed at this moment.
#'                  }
#'                  
#' @param alpha Allocated alpha level in incremental fashion. The sum of alpha 
#' is the total type I error allocated to the group sequential test, 
#' usually 0.025 1-sided.
#' @param side  Side of test 1 or 2. Default side = 1. 
#'              Currently, only 1-sided test is implemented.
#' @param  f.ws  Weight functions of survival rate used for each analysis.
#'               For example (1), if there are 3 analyses planned including IA1, IA2, and FA.
#'               IA1 uses log-rank test, IA2 uses a maxcombo test of (logrank, FH01),
#'               and FA uses a maxcombo test of (logrank, FH01, FH11). Then specify as
#'               f.ws = list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)),
#'               where define lr = function(s){1}; fh01=function(s){1-s}; fh11 = function(s){s*(1-s)};
#'               For example (2), if only fh01 is used for all three analyses IA1, IA2 and FA.
#'               f.ws = list(IA1=list(fh01), IA2=list(fh01), FA1=list(fh01)).
#'               For example (3), if only logrank is used for a single time analysis, then
#'               f.ws = list(IA1=list(lr)).
#'
#' @return An object with dataframes below.
#'  \itemize{
#'  \item  test.results:  The test results for each analysis and statistical 
#'               inference whether it is considered statistically significant.
#'  \item  corr: Correlation matrix among the weighted log-rank tests at different analyses.
#'  \item  wt: Weight functions used for each analysis
#'  \item  setting: alpha and side of test. 
#'  }
#'  
#' @examples 
#' 
#' N=600; m0 = 12; A=21; r=1; hr = 0.65; w = 1.5; dropOff0 = dropOff1 = 0; 
#' targetEvents = c(300, 397, 496); cuts = 6
#' lambda0 = log(2) / m0; lambda1 = c(log(2)/m0, log(2)/m0*hr)
#' 
#' data0 = sim.pwexp(nSim=1, N = N, A = A, w=w, r=r, lambda0=lambda0, 
#'       lambda1=lambda1, cuts=cuts, dropOff0= dropOff0, dropOff1= dropOff1, 
#'       targetEvents = targetEvents)
#'       
#' 
#' data1 = data0[[1]][sim==1,]; data2 = data0[[2]][sim==1,]; data3 = data0[[3]][sim==1,]
#' #Add strata variables       
#' data1$strata1 = data1$strata2 =data1$strata3 =sample(c(1,2), N, replace = TRUE);
#' data2$strata1 = data2$strata2 =data2$strata3 =sample(c(1,2), N, replace = TRUE);
#' data3$strata1 = data3$strata2 =data3$strata3 =sample(c(1,2), N, replace = TRUE);
#' data1$group = as.numeric(data1$treatment == "experimental")
#' data2$group = as.numeric(data2$treatment == "experimental")
#' data3$group = as.numeric(data3$treatment == "experimental")
#' 
#' km.IA1<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data1)
#' plot(km.IA1,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.IA2<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data2)
#' plot(km.IA2,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' km.FA<-survival::survfit(survival::Surv(survTimeCut,1-cnsrCut)~treatment,data=data3)
#' plot(km.FA,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' 
#' #Define weight functions for weighted log-rank tests
#' lr = function(s){1}
#' fh01 = function(s){(1-s)}
#' fh11 = function(s){s*(1-s)}
#' #stabilized FH(0, 1; 0.5)
#' sfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1-s1)} 
#' #modestly log-rank
#' mfh01 = function(s){s1 = apply(cbind(s, 0.5), MARGIN=1,FUN=max); return(1/s1)}
#' 
#' #Example (1). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and FH01);
#' #               FA uses max(log-rank, FH01, FH11).
#' 
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr, fh01), 
#'          FA=list(lr, fh01, fh11)))
#'   
#' #Example (2a). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #             FA uses max(log-rank, FH01). Unstratified analysis
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = NULL, strata2 = NULL, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr, sfh01), 
#'          FA=list(lr, fh01)))
#'          
#' #Example (2b). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #              FA uses max(log-rank, FH01). 
#' #              Stratified analysis of 2 stratification factors
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr, sfh01), 
#'          FA=list(lr, fh01)))
#'          
#' #Example (2c). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #              FA uses max(log-rank, FH01). 
#' #              Stratified analysis of 3 stratification factors
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = data1$strata3,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr, sfh01), 
#'          FA=list(lr, fh01)))
#'          
#' #Example (3). 2 IAs and FA. IA1 uses log-rank test; IA2 uses max(log-rank and sFH01);
#' #               FA uses max(log-rank, FH01). 
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr, sfh01), 
#'          FA=list(lr, fh01)))
#'          
#' #Example (4). 2 IAs and FA. All analyses use max(logrank, sfh01).
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr, sfh01), 
#'          IA2=list(lr, sfh01), 
#'          FA=list(lr, sfh01)))
#'          
#'  #For each analysis, if the wlr.maxcombo() function can perform
#'  #stratified maxcombo test and produce the equivalent p value.
#'  #However, wlr.maxcombo() function cannot provide statistical inference
#'  #whether it is considered as statistically significant for group
#'  #sequential design.
#'          
#' wlr.maxcombo(time=data1$survTimeCut, event=1-data1$cnsrCut, group=data1$group, 
#'  rho=NULL, gamma=NULL, tau = NULL, s.tau=NULL,
#'  strata1 = data1$strata1, strata2=data1$strata2, strata3=data1$strata3,
#'  f.ws=list(lr, sfh01), side = 1)
#' 
#' #Example (5). 2 IAs and FA. All analyses use logrank.
#'               
#' wlr.inference(datasets=list(IA1=data1, IA2=data2, FA=data3), 
#'     alpha = c(0.01, 0.02, 0.02)/2, side = 1,
#'     strata1 = data1$strata1, strata2 = data1$strata2, strata3 = NULL,
#'     f.ws=list(IA1=list(lr), 
#'          IA2=list(lr), 
#'          FA=list(lr)))
#'         
#'  #For each analysis, wlr() function can perform stratified weighted
#'  #log-rank test for a single weighted log-rank test and produces equivalent p value.
#'  #However, for group sequential design, wlr.maxcombo() function cannot 
#'  #provide statistical inference whether the test result is considered as statistically significant.
#'    
#'  wlr(time=data1$survTimeCut, event=1-data1$cnsrCut, group=data1$group, 
#'  rho=NULL, gamma=NULL, tau = NULL, s.tau=NULL,
#'  strata1 = data1$strata1, strata2=data1$strata2, strata3=data1$strata3,
#'  f.ws=lr, side=1)
#'  
#' @export
#' 
#' 
#'  
wlr.inference = function(datasets=list(IA1=data1, IA2=data2, FA=data3), 
        alpha = c(0.01, 0.02, 0.02)/2, side = 1, 
        strata1 = NULL, strata2 = NULL, strata3 = NULL,
        f.ws=list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11))){

  #K analyses
  K=length(f.ws)
  
  #Number of test components in each analysis
  J = lengths(f.ws)
  
  #Find the complete correlation matrix among all weighted log-rank 
  #tests J[1]+J[2]+...+J[K] dimensions
  
  corr = matrix(1, nrow=sum(J), ncol=sum(J))
  #calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[i]){
      for (ip in i:K){
        for (jp in 1:J[ip]){
          row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the corr matrix
          
          #incremental location pointer for column compared to row of the corr matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            #information matrix for Zij and Zi'j'
            datai = datasets[[i]]; dataip = datasets[[ip]]

            corr[row, col] = wlr.cov2t(time1=datai$survTimeCut, event1=1-datai$cnsrCut, 
                            time2=dataip$survTimeCut, event2=1-dataip$cnsrCut, 
                            group=datai$group, 
                            strata1=strata1, strata2=strata2, strata3=strata3,
                            f.ws1=f.ws[[i]][[j]], f.ws2=f.ws[[ip]][[jp]])$corr
            corr[col, row] = corr[row, col]
          }
        }
      }
    }
  }
  
  #Rejection boundary: recursively solve for the rejection boundary for each analysis
  b = rep(NA, K)
  
  #First Analysis
  #maxcombo has at least 2 components.
  if (J[1] >= 2){
    if(side == 1){
      #1-sided test
      f.b = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
            mean=rep(0, J[1]), corr = corr[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      b[1] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
    } 
  } else{
    if (side == 1) {b[1] = qnorm(1-alpha[1])} else {b[1] = qnorm(1-alpha[1]/2)}
  }
  
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
      if(side == 1){
        #1-sided test
        f.b = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b[j], J[j]))}}
          
          idx1 = 1:sum(J[1:(i-1)])
          if (length(LL1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, mean=rep(0, length(idx1)), 
                                  corr = corr[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #index of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, mean=rep(0, length(idx)), 
                                corr = corr[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        b[i] = uniroot(f=f.b, interval=c(1, 20), tol = 1e-8)$root
      }
    }
  }
  
  #Calculate p value for each analysis
  test.results = NULL
  for (i in 1:K){
    if(side == 1){
      testi = wlr.maxcombo(time=datasets[[i]]$survTimeCut, event=1-datasets[[i]]$cnsrCut,
              group=datasets[[i]]$group, 
              strata1=strata1, strata2=strata2, strata3=strata3, 
              rho = NULL, gamma=NULL, tau = NULL, s.tau=NULL,
              f.ws=f.ws[[i]], side = "one.sided")$test.results
      testi$analysis = i; testi$z.bound = b[i]; 
      
      testi$inference = ifelse(testi$z.max > b[i], "Positive", "Negative")
      test.results = rbind(test.results, testi)
    }
  }
  
  o = list()
  o$corr = corr
  o$test.results = test.results 
  o$wt = f.ws
  setting = list(alpha=alpha, side=side)
  o$setting = setting
  return(o)
}
