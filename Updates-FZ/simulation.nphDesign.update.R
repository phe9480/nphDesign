################## Updated version of simulation.nphDesign to include parallelisation
## Author: Michael Sweeting
## Number of cores to use is an argument of function (n.cores). Defaults to 10
simulation.nphDesign.update = function(nSim=10000, N = 672, A = 21, w=1.5, r=1, lambda0=log(2)/11.7, lambda1=log(2)/11.7*0.745, 
                                        cuts=NULL, dropOff0=0, dropOff1=0, targetEvents = c(290, 397, 496), 
                                        sf = "LDOF", overall.alpha = 0.025, side = 1, alpha = NULL,
                                        logrank="N", fws.options=NULL, H0 = "N", nphDesign = NULL, n.cores = 10) {
  
  ### Parallelisation of simulation.nphDesign
  library(parallel)   #For PSOCK parallel backend
  library(foreach)    #For 'foreach' loop
  library(doParallel) #For using PSOCK cluster with operator '%dopar%'
  n.cores <- n.cores
  my.cluster <- makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  clusterEvalQ(my.cluster, .libPaths(.libPaths())
  )
  registerDoParallel(cl = my.cluster)
  clusterEvalQ(my.cluster, library(nphDesign))
  
  ### create a big list saving simulated data
  sim.data<-vector("list",length=nSim)
  
  ##############
  #Extract the setting parameters from nphDesign object if provided
  ##############
  if(!is.null(nphDesign)){
    nph = f.extract(nphDesign)
    N = nph$N; A=nph$A;  w=nph$w; r=nph$r; lambda0=nph$lam0; lambda1 = nph$lam1;
    cuts = nph$cuts; targetEvents = ceiling(nph$targetEvents); overall.alpha = nph$overall.alpha;
    side = nph$side; alpha = nph$alpha;
    if(is.null(fws.options)){fws.options = list(nph$fws)}
  }
  
  #Simulation for checking type I error
  if (H0 == "Y"){lambda1 = lambda0} 
  
  ##############################
  #M Options of test strategies
  M = length(fws.options)
  
  #K analyses
  K=length(fws.options[[1]])
  
  timing = targetEvents / targetEvents[K]
  
  #if alpha is not provided, use sf to derive alpha. 
  #if alpha is provided, then sf is ignored.
  if(is.null(alpha) && !is.null(overall.alpha)){
    ld.obf = function(s){
      if (side == 1){a = 2*(1 - pnorm(qnorm(1-overall.alpha/2)/sqrt(s)))}
      if (side == 2){a = 2*2*(1 - pnorm(qnorm(1-overall.alpha/4)/sqrt(s)))}
      return(a)
    }
    ld.pk = function(s){overall.alpha * log(1 + (exp(1)-1)*s)}
    
    
    if (sf == "LDOF"){
      gs.alpha = ld.obf(s = timing)
    }
    if (sf == "LDPK") {
      gs.alpha = ld.pk(s = timing)
    }
    if (K == 1){alpha = overall.alpha} else{
      alpha[1] = gs.alpha[1]
      for(i in 2:K){alpha[i] = gs.alpha[i] - gs.alpha[i-1]}
    }
  }
  
  wlr.sim = array(NA, dim=c(nSim, M, K, 5))
  if (logrank == "Y") {lr.sim = array(NA, dim=c(nSim, K, 5))}
  
  for (i in 1:nSim) {
    #(1). Generate data
    dati = sim.pwexp(nSim=1, N = N, A = A, w=w, r=r, lambda0=lambda0, lambda1=lambda1, 
                     cuts=cuts, dropOff0=dropOff0, dropOff1=dropOff1, targetEvents = targetEvents)
    #(2). Add group variable
    for (j in 1:K) {
      dati[[j]]$group = as.numeric(dati[[j]]$treatment == "experimental")
    }  
    
    sim.data[[i]]<-dati
  }
  
  #(3). Testing strategy m
  for (m in 1:M){    
    #Perform weighted log-rank test for each analysis in strategy m
    ## This next line takes a long-time to compute, making simulation difficult. Use parallel computing instead
    wlr <- foreach(i = 1:nSim) %dopar% {
      wlr.inference(datasets=sim.data[[i]], alpha = alpha, side = side, f.ws=fws.options[[m]])$test.results
    }
    stopCluster(cl = my.cluster) ## Stop cluster
    for (i in 1:nSim) {
      wlri = wlr[[i]][!duplicated(wlr[[i]]$analysis),]
      wlri$result = as.numeric(wlri$inference=="Positive")
      wlr.sim[i, m, , ] = as.matrix(wlri[,c(2,3,5,6,8)])
    }
  }    
  
  #(4). Standard log-rank test for all analyses if requested
  if (logrank=="Y"){
    if(K > 1){
      #GSD boundary for each analysis
      if(side == 1) {
        z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha,timing=timing[1:(K-1)], 
                                   sfu=gsDesign::sfLDOF)$upper$bound 
      } else{
        z.bd <- gsDesign::gsDesign(k=K,  alpha=overall.alpha/2,timing=timing[1:(K-1)], 
                                   sfu=gsDesign::sfLDOF)$upper$bound 
      }
    } else {
      if(side == 1) {z.bd = qnorm(1-overall.alpha)} else{z.bd = qnorm(1-overall.alpha/2)}
    }
    for (j in 1:K){
      for (i in 1:nSim) {
        lr.test = survival::survdiff(survival::Surv(survTimeCut, 1-cnsrCut) ~ group, data = sim.data[[i]][[j]])
        
        #convert to z value in correct direction: z>0 means better experimental arm.
        better = as.numeric(lr.test$obs[1] > lr.test$obs[2])
        sign = 2*better - 1
        z = sqrt(lr.test$chisq) * sign
        
        #count power
        lr.sim[i, j, 1] = z
        if (side == 1){ p = 1 - pnorm(z)} else {p = 2*(1 - pnorm(z))}
        
        lr.sim[i, j, 2] = p
        lr.sim[i, j, 3] = j
        lr.sim[i, j, 4] = z.bd[j]
        lr.sim[i, j, 5] = as.numeric(z > z.bd[j])
      }
    }
  }
  
  pow = matrix(NA, nrow=M, ncol=K)
  overall.pow = rep(0, M)
  
  for(m in 1:M){for(j in 1:K){pow[m, j] = sum(wlr.sim[,m,j,5])/nSim}}
  for(m in 1:M){for(i in 1:nSim){
    overall.pow[m] =  overall.pow[m] + as.numeric(sum(wlr.sim[i,m,,5])>0)
  }}
  overall.pow = overall.pow / nSim
  
  o=list()
  o$power = pow; o$overall.power = overall.pow
  o$wlr.simulations = wlr.sim
  if(logrank=="Y"){
    lr.pow = rep(NA, K)
    for (j in 1:K) {lr.pow[j] = sum(lr.sim[,j,5])/nSim}
    
    lr.overall.pow = 0
    for (i in 1:nSim){lr.overall.pow = lr.overall.pow + as.numeric(sum(lr.sim[i,,5])>0)}
    lr.overall.pow = lr.overall.pow / nSim
    
    o$lr.overall.power = lr.overall.pow
    o$lr.power = lr.pow
    o$lr.simulations = lr.sim
  }
  o$sim.data<-sim.data
  return(o)
}