#'
#' Extract The Design Parameters From nphDesign Object
#' Extract the design parameters from nphDesign object generated from 
#' finalize.nphDesign or explore.nphDesign or wlr.power.maxcombo functions. This
#' is a utility function for running trial simulations and display the nphDesign.
#'
#' @param  nphDesign nphDesign object generated from finalize.nphDesign or explore.nphDesign or wlr.power.maxcombo functions.
#' 
#' @return An object with design parameters.
#' 
#' nph = D0.A18X1.W1
#' 
#' f.extract(nphDesign=H0.A18X1.W1)
#' 
#' #@export
#' 
f.extract = function(nphDesign){
  tt = 1:100;
  #Find hazard rate lambda
  f.lam.cut = function(h){
    LAM = rep(NA, 100); t = 1:100
    for(i in 1:100){LAM[i] = h(i)}
    lam = LAM[!duplicated(LAM)];
    cut = t[!duplicated(LAM)];
    ans = list(); ans$lam = lam; ans$cut = cut;
    return(ans)
  }
  o=list()
  
  #Check the format of nphDesign object.
  if(is.null(nphDesign$design)){
    #nphDesign object is created by wlr.power.maxcombo() function
    acc.tt = nphDesign$setting$F.entry(tt)
    A = tt[length(unique(acc.tt))]
    w = log(nphDesign$setting$F.entry(A/2)) / log(0.5)
    N = nphDesign$n; 
    r = nphDesign$setting$r
    m0 = nphDesign$medians$m0
    m1 = nphDesign$medians$m1

    lam.cut0 = f.lam.cut(nphDesign$setting$h0)
    lam.cut1 = f.lam.cut(nphDesign$setting$h1)
    
    o$targetEvents = nphDesign$events$events
    o$events = nphDesign$events
    o$maturity = nphDesign$maturity

    if(!is.null(nphDesign$Critical.Values.HR.H1)){
      o$CV.HR = nphDesign$Critical.Values.HR.H1
      o$CV.med = nphDesign$Critical.Values.Medians.H1
    } else {
      o$CV.HR = o$CV.med = NULL
    }
    
    o$alpha = nphDesign$setting$alpha; 
    o$overall.alpha = sum(nphDesign$setting$alpha)
    o$m0 = m0; o$m1 = m1
    o$power = nphDesign$power; o$overall.power = nphDesign$overall.power
    
    o$side=nphDesign$setting$side
    o$fws = nphDesign$wt
    o$DCO = nphDesign$T
    o$b = nphDesign$bounds$b.Hp
    o$AHR = nphDesign$AHR
  }else{
    #nphDesign object is created by finalize.nphDesign() / explore.nphDesign() functions.
    #Find the accural period from F.entry function
    acc.tt = nphDesign$design$setting$F.entry(tt)
    A = tt[length(unique(acc.tt))]
    w = log(nphDesign$design$setting$F.entry(A/2)) / log(0.5)
    
    N = nphDesign$design$n; 
    r = nphDesign$design$setting$r
    
    m0 = nphDesign$design$medians$m0
    m1 = nphDesign$design$medians$m1

    lam.cut0 = f.lam.cut(nphDesign$design$setting$h0)
    lam.cut1 = f.lam.cut(nphDesign$design$setting$h1)
    
    o$targetEvents = nphDesign$design$events$events
    o$events = nphDesign$design$events
    o$maturity = nphDesign$design$maturity
    
    if(!is.null(nphDesign$design$Critical.Values.HR.H1)){
      o$CV.HR = nphDesign$design$Critical.Values.HR.H1
      o$CV.med = nphDesign$design$Critical.Values.Medians.H1
    } else {
      o$CV.HR = o$CV.med = NULL
    }
    
    o$alpha = nphDesign$design$setting$alpha; 
    o$overall.alpha = sum(nphDesign$design$setting$alpha)
    o$m0 = m0; o$m1 = m1
    o$power = nphDesign$design$power; 
    o$overall.power = nphDesign$design$overall.power
    o$side=nphDesign$design$setting$side
    o$fws = nphDesign$design$wt
    o$DCO = nphDesign$design$T
    o$b = nphDesign$design$bounds$b.Hp
    o$AHR = nphDesign$design$AHR
  }

  lam0 = lam.cut0$lam; cut0 = lam.cut0$cut
  lam1 = lam.cut1$lam; cut1 = lam.cut1$cut
  
  if(length(cut0)==1 && cut0 == 1){cut0=NULL} else {cut0 = cut0[cut0!=1]}; 
  if(length(cut1)==1 && cut1 == 1){cut1=NULL} else {cut1 = cut1[cut1!=1]}; 
  cuts = sort(unique(c(cut0, cut1)))
  
  #Number of subintervals
  K = length(cuts) + 1
  cuts.star = c(0, cuts, Inf)
  c0.star = c(0, cut0, Inf)
  c1.star = c(0, cut1, Inf)
  L0 = length(c0.star); L1 = length(c1.star)
  
  if(K > 1){
    lambda0 = lambda1 = rep(NA, K)
    for (i in 1:K){
      if(length(cut0) == 0){lambda0[i] = lam0[1]} else{
        ix0 = cuts.star[i+1] <= c0.star[2:L0] & cuts.star[i+1] > c0.star[1:(L0-1)]
        lambda0[i] = lam0[ix0]
      }
      if(length(cut1) == 0){lambda1[i] = lam1[1]} else{
        ix1 = cuts.star[i+1] <= c1.star[2:L1] & cuts.star[i+1] > c1.star[1:(L1-1)]
        lambda1[i] = lam1[ix1]
      }  
    }
  } else {
    lambda1 = lam1; lambda0 = lam0
  }
  
  o$N = N; o$A=A; o$w=w; o$r=r;
  o$lam0 = lambda0; o$lam1=lambda1;
  o$cuts = cuts;
  
  return(o)
}

