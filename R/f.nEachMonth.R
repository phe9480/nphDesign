#' Determine number of subjects randomized for every month according to specified weight function
#' 
#' @param N Total number of subjects in a trial
#' @param A Duration of accrual period in months
#' @param w Weigth parameter in the non-uniform enrollment pattern, culmulative proportion of enrollment = (t/A)^w
#' @param r Randomization ratio r:1
#' @return n0 and n1 for number of subjects for each month in two groups (n0: control; n1: experimental)
#' @examples 
#' #Enroll 600 patients in 24 months with accrual weight 2 and 2:1 randomization
#' f.nEachMonth(N=600, A=24, w=2, r=2)
#' @export
f.nEachMonth = function (N=600, A=24, w=2, r=2) {
  
  N1 = N * (r/(r+1))
  N0 = N - N1
  
  #When r > 1, the control arm has smaller number of pts. 
  #Just need to determine enrollment for control arm per month, 
  #then to obtain enrollment for experimental arm by n1i = n0i * r.
  
  n1 = n0 = rep(NA, A) #enrollment by month
  randdt1 = rep(NA, N1) #randomization date
  randdt0 = rep(NA, N0)
  
  #Determine number of pts per month for control arm
  #(i-1)th month cumulative enrolled pts
  cLastN0 = 0
  for (i in 1:A) {
    #ith month: cumulative #pts
    cN0i = max(round((i/A)^w * N0), 1)
    
    n0[i] = max(cN0i - cLastN0, 1)
    if (i == A) {n0[i] = N0 - sum(n0[1:(A-1)]) }
    cLastN0 = cN0i  
  }
  n1 = n0 * r
  
  o = list()
  o$n0 = n0
  o$n1 = n1
  return(o)
}
