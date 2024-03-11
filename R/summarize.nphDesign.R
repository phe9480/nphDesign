#'
#' Summary of Study Design Object nphDesign Produced by finalize.nphDesign or wlr.power.maxcombo function
#' 
#' Extract the design parameters from nphDesign object generated from 
#' finalize.nphDesign or explore.nphDesign or wlr.power.maxcombo functions. This
#' is a utility function for running trial simulations and display the nphDesign.
#'
#' @param  nphDesign nphDesign object generated from finalize.nphDesign or explore.nphDesign or wlr.power.maxcombo functions.
#' 
#' @return A nice formatted study design summary ready for formal documentation
#' 
#' 
#' @export
#' 
summarize.nphDesign = function(nphDesign = PFS.design, endpoint="Endpoint", caption="Summary of Statistical Design"){
  nph = f.extract(nphDesign)
  n = nph$N; A=nph$A;  w=nph$w; r=nph$r; 
  lambda0=nph$lam0; lambda1 = nph$lam1;
  cuts = nph$cuts; targetEvents = nph$targetEvents; 
  DCO = nph$DCO
  maturity = nph$maturity
  CV.HR = nph$CV.HR; CV.med = nph$CV.med
  b = nph$b
  
  overall.alpha = nph$overall.alpha;
  overall.power = nph$overall.power
  power = nph$power
  m0 = nph$m0; m1 = nph$m1
  
  #delayed effect
  if(!is.null(nph$cuts)) {delay = cuts[1]; AHR = nph$AHR; HR = NULL} else{
    #Prop. Hazards
    delay = NULL; HR = nph$AHR[1]; AHR = NULL
  }
  if(!is.null(nph$cuts)){delay = cuts[1]}
  
  #Ignore delayed effect in display
  delay = NULL
  side = nph$side; alpha = nph$alpha; f.ws = nph$fws
  
  #Number of analyses
  K = length(f.ws)
  J = lengths(f.ws) 
  timing = targetEvents/targetEvents[K]
  
  row.endpoint = c("Endpoint", endpoint)
  row.n = c("N", n)
  row.r = c("Randomization Ratio", paste(r, ": 1"))
  row.A = c("Accrual", paste(A, "months"))
  row.overall.alpha = c("Overall Alpha", 
                        paste(round(overall.alpha,4), "(", side, " sided)", sep=""))
  row.overall.power = c("Overall Power", paste(round(overall.power*100,1), "%", sep=""))  
  row.median = c("Median (Control vs. Exp. Arm)", 
                 paste(round(m0,1), "vs", round(m1,1)))
  
  
  if(!is.null(delay)){
    row.delay = c("Delayed Effect", delay)
    row.target.HR = NULL
  }else{
    row.delay = NULL
    row.target.HR = c("Target HR", HR)
  }
  
  #dataframe for overall study setting
  tab = data.frame(
    rbind(
      row.endpoint, 
      row.n,
      row.r,
      row.A,
      row.overall.alpha,
      row.overall.power,
      row.median,
      row.delay,
      row.target.HR
    )
  )
  
  if(K > 1) {tab$Analysis = ""}
  
  #summary of study design per analysis
  if(K > 1){
    #group sequential design
    for (i in 1:K){
      row.DCO = c("Expected Time of Analysis", round(DCO[i], 1))
      row.IF = c("Information Fraction", round(timing[i],2))
      row.events = c("Target Events (% maturity) for Analysis",
                     paste(round(targetEvents[i]), 
                           "(", round(maturity$overall.maturity[i]*100,1), "%)", sep=""))
      
      if(!is.null(AHR)){row.AHR = c("Expected Average HR", 
                                    round(AHR[i], 3))}else{row.AHR=NULL}
      
      if(!is.null(CV.HR)){
        med.imp = CV.med[i] - m0
        CV.HR.med = paste(round(CV.HR[i], 3), "(", 
                          round(med.imp, 1), ")", sep="")
        row.CV = c("Critical Value in HR (median imp.)", CV.HR.med)
      } else {row.CV = NULL}
      
      #alpha and power
      if(side == 1) {
        row.alpha = c("p value bound (alpha)", round(1-pnorm(b[i]), 4))
      }else{
        row.alpha = c("p value bound (alpha)", round(2*(1-pnorm(b[i])), 4))
      }
      if (round(power[i]*100, 1)==100){power.txt = "> 99%"}else{
        power.txt = paste(round(power[i]*100, 1), "%", sep="")
      }
      row.power = c("Power", power.txt)
      
      tabi = data.frame(rbind(
        row.DCO,
        row.IF,
        row.events,
        row.AHR,
        row.CV,
        row.alpha,
        row.power
      ))
      if(K == 2 && i == 1){tabi$Analysis = "Interim Analysis"}
      if(K == 2 && i == 2){tabi$Analysis = "Final Analysis"}
      if(K > 2 && i < K){tabi$Analysis = paste("Interim Analysis", i)}
      if(K > 2 && i == K){tabi$Analysis = "Final Analysis"}
      
      tab = rbind(tab, tabi)
    }
  } else{
    #Only 1 analysis
    row.DCO = c("Expected Time of Analysis", round(DCO, 1))
    row.events = c("Number of Target Events (% maturity) for Analysis",
                   paste(round(targetEvents), 
                         "(", round(maturity$overall.maturity*100,1), "%)", sep=""))
    if(!is.null(CV.HR)){
      med.imp = CV.med - m0
      CV.HR.med = paste(round(CV.HR, 3), "(", 
                        round(med.imp, 1), ")", sep="")
      row.CV = c("Critical Value in HR (median imp.)", CV.HR.med)
    } else {row.CV = NULL}
    
    #alpha and power
    if(side == 1) {
      row.alpha = c("p value bound (alpha)", round(1-pnorm(b), 4))
    }else{
      row.alpha = c("p value bound (alpha)", round(2*(1-pnorm(b)), 4))
    }
    if (round(power*100, 1)==100){power.txt = "> 99%"} else{
      power.txt = paste(round(power*100, 1), "%", sep="")
    }
    row.power = c("Power", power.txt)
    
    tabi = data.frame(rbind(
      row.DCO,
      row.events,
      row.CV,
      row.alpha,
      row.power
    ))
    tab = rbind(tab, tabi)
  } 
  tab = dplyr::rename(tab, Description = X1, Design = X2)
  
  if(K > 1){    
    out = tab %>% group_by(Analysis) %>% gt()
  } else{
    out = tab %>% gt()
  }
  
  out = out %>% tab_header(
    title = md(caption), 
  ) %>% tab_source_note(
    source_note = md("Time unit is month.")
  )
  return(out)
}

