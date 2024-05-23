## FZ: Expand the upper bound for root finding in lines 321 and 322
wlr.power.maxcombo.update = function(T = c(24, 36, 48), events = NULL, 
                                     alpha=c(0.01, 0.02, 0.02)/2, power = 0.9, side = 1, r = 1, n = NULL, 
                                     h0 = function(t){log(2)/12}, S0= function(t){exp(-log(2)/12 * t)},
                                     h1 = function(t){log(2)/12*0.70}, S1= function(t){exp(-log(2)/12 * 0.7 * t)}, 
                                     f.logHR = function(t){log(as.numeric(t<6) + as.numeric(t>= 6)*0.65)},
                                     f.ws = list(IA1=list(lr), IA2=list(lr, fh01), FA=list(lr, fh01, fh11)),
                                     F.entry = function(t){(t/18)*as.numeric(t <= 18) + as.numeric(t > 18)}, 
                                     G.ltfu = function(t){0}, non.centrality = "Heetal2021", show.setting="N"){
  
  #Re-parameterization as consistent with the manuscript; r1 is proportion of experimental arm subjects.
  r1 = r / (r + 1); r0 = 1 - r1 
  
  #K analyses
  K=length(f.ws)
  
  #Number of test components in each analysis
  J = lengths(f.ws)
  
  #If n not provided but T provided. length(T) must be K.
  #Determine an initial n using log-rank test at FA
  if(is.null(n)){
    if(non.centrality != "Schoenfeld"){
      R = wlr.mu(T = T[K], r = r, n = NULL, h0 = h0, S0=S0,
                 h1 = h1, S1=S1, f.logHR = f.logHR,
                 rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
                 F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }else{
      R = wlr.mu.schoenfeld(T = T[K], r = r, n = NULL, h0 = h0, S0=S0,
                            h1 = h1, S1=S1, f.logHR = f.logHR,
                            rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws = NULL,
                            F.entry = F.entry, G.ltfu = G.ltfu)$R2
    }
    if(side == 1){za = qnorm(1-sum(alpha))} else {za = qnorm(1-sum(alpha)/2)}
    n = ceiling(((za + qnorm(power))/R)^2 )
  }
  
  #If T is not provided but the events are provided. length(events) must be K.
  #Determine the calendar times for the required number of events
  if(!is.null(events)){
    T = rep(NA, length(events))
    for (i in 1:length(events)){
      f.root = function(t){
        events[i] - f.nEvents(T = t, r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events$n.events.total
      }
      T[i] = uniroot(f.root, interval= c(1, 1000), tol = 1e-8)$root
    }
  }
  
  #If events are not provided, but T is provided. length(T) must be K.
  #Calculate the events according to T.
  events0 = events1 = rep(NA, K)
  if (is.null(events)){events = rep(NA, K)}
  
  for (i in 1:K){
    n.tmp = f.nEvents(T = T[i], r = r, h0 = h0, S0 = S0, h1 = h1, S1 = S1, F.entry = F.entry, G.ltfu = G.ltfu, n = n)$n.events
    
    if (is.na(events[i])){events[i] = n.tmp$n.events.total}
    events0[i] = n.tmp$n.events0
    events1[i] = n.tmp$n.events1
  }
  
  #Non-centrality parameters for all WLRTs in all analyses
  mu = matrix(NA, nrow=K, ncol=max(J))
  
  for(i in 1:K){
    for(j in 1:J[[i]]){
      if (non.centrality != "Schoenfeld") {
        non.centrality = "Heetal2021"
        mu[i, j] <- wlr.mu(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
                           rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[i]][[j]],
                           F.entry = F.entry, G.ltfu = G.ltfu)$mu2
      }else{
        tt = wlr.mu.schoenfeld(T = T[i], r = r, n = n, h0 = h0, S0=S0, h1 = h1, S1=S1, f.logHR = f.logHR,
                               rho = NULL, gamma = NULL, tau = NULL, s.tau = NULL, f.ws = f.ws[[i]][[j]],
                               F.entry = F.entry, G.ltfu = G.ltfu)
        mu[i, j] = c(tt$mu2)
      }
    }
  }
  
  #Complete correlation matrix structure: J[1]+J[2]+...+J[K] dimensions
  corr.Hp = corr.strict.H0 =matrix(1, nrow=sum(J), ncol=sum(J))
  #calculate the correlation between Zij and Z_i'j'
  for (i in 1:K){
    for (j in 1:J[[i]]){
      for (ip in i:K){
        for (jp in 1:J[ip]){
          row = as.numeric(i>=2)*sum(J[1:(i-1)])+j #row location of the corr matrix
          
          #incremental location pointer for column compared to row of the corr matrix
          incr = (as.numeric(ip>=2)*sum(J[1:(ip-1)])+jp)-(as.numeric(i>=2)*sum(J[1:(i-1)])+j)
          col = row + incr
          #incr controls the computation only limited to upper right corner
          if(incr > 0){
            #information matrix for Zij and Zi'j'
            info = wlr.info(T = c(T[i], T[ip]), r = r, n = NULL, h0 = h0, S0= S0, h1 = h1, S1=S1, 
                            rho=NULL, gamma=NULL, tau=NULL, s.tau=NULL,
                            f.ws=c(f.ws[[i]][[j]],f.ws[[ip]][[jp]]), F.entry = F.entry, G.ltfu = G.ltfu)
            corr.Hp[row, col] = info$corr.Hp
            corr.strict.H0[row, col] = info$corr.H0
            corr.Hp[col, row]  = corr.Hp[row, col]
            corr.strict.H0[col, row] = corr.strict.H0[row, col]
          }
        }
      }
    }
  }
  
  #Rejection boundary: recursively solve for the rejection boundary for each analysis
  b.Hp = b.H0 = rep(NA, K)
  
  #First Analysis
  #maxcombo has at least 2 components.
  if (J[[1]] >= 2){
    if(side == 1){
      #1-sided test
      f.b.Hp = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.Hp[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      f.b.H0 = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-Inf,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.strict.H0[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }      
      b.Hp[1] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
      b.H0[1] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
    } else {
      #2-sided test
      f.b.Hp = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.Hp[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }
      f.b.H0 = function(x){
        1-mvtnorm::pmvnorm(lower=rep(-x,J[1]),upper=rep(x, J[1]), 
                           mean=0, corr = corr.strict.H0[1:J[1],1:J[1]], abseps = 1e-8, maxpts=100000)[1] - alpha[1]
      }      
      b.Hp[1] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
      b.H0[1] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
    }
  } else{
    if (side == 1) {b.Hp[1] = b.H0[1] = qnorm(1-alpha[1])} else {b.Hp[1] = b.H0[1] = qnorm(1-alpha[1]/2)}
  }
  
  #Recursively solve other boundaries from 2nd analysis to Kth analysis
  if(K > 1){
    for(i in 2:K){
      if(side == 1){
        #1-sided test
        f.b.Hp = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b.Hp[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}
          
          idx1 = 1:sum(J[1:(i-1)])
          if (length(UU1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                                  mean=rep(0, length(idx1)),corr = corr.Hp[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                                mean=rep(0, length(idx)), corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        f.b.H0 = function(x){
          LL1 = rep(-Inf, sum(J[1:(i-1)]))
          UU1 = rep(b.Hp[1], J[1])
          if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}
          
          idx1 = 1:sum(J[1:(i-1)])
          if (length(LL1) == 1) {
            P1 = pnorm(UU1)
          } else {
            P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, 
                                  mean=rep(0, length(idx1)), 
                                  corr = corr.strict.H0[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
          }
          
          LL2 = rep(-Inf, sum(J[1:i]))
          UU2 = c(UU1, rep(x, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, 
                                mean=rep(0, length(idx)), 
                                corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(P1 - P2 - alpha[i])
        }
        b.Hp[i] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0[i] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      } else {
        #2-sided test
        f.b.Hp = function(x){
          LL = rep(-b.Hp[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.Hp[j], J[j]))}}
          LL = c(LL, rep(x, J[i]))
          
          UU = rep(b.Hp[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.Hp[j], J[j]))}}
          UU = c(UU, rep(Inf, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          p1 = mvtnorm::pmvnorm(lower=LL, upper=UU, mean=rep(0, length(idx)), 
                                corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          LL2 = LL; UU2 = UU
          LL2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-Inf, J[i])
          UU2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-x, J[i])
          p2 = mvtnorm::pmvnorm(lower=LL2, upper=U2, mean=rep(0, length(idx)), 
                                corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(p1+p2-alpha[i])
        }
        f.b.H0 = function(x){
          LL = rep(-b.H0[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.H0[j], J[j]))}}
          LL = c(LL, rep(x, J[i]))
          
          UU = rep(b.H0[1], J[1])
          if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.H0[j], J[j]))}}
          UU = c(UU, rep(Inf, J[i]))
          
          idx = 1:sum(J[1:i]) #indices of the corr matrix
          p1 = mvtnorm::pmvnorm(lower=LL, upper=UU, mean=rep(0, length(idx)), corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          LL2 = LL; UU2 = UU
          LL2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-Inf, J[i])
          UU2[(sum(J[1:i])-sum(J[1:(i-1)])):sum(J[1:i])] = rep(-x, J[i])
          p2 = mvtnorm::pmvnorm(lower=LL2, upper=U2, mean=rep(0, length(idx)), corr = corr.strict.H0[idx, idx], abseps = 1e-8, maxpts=100000)[1]
          return(p1+p2-alpha[i])
        }
        b.Hp[i] = uniroot(f=f.b.Hp, interval=c(1, 20), tol = 1e-8)$root
        b.H0[i] = uniroot(f=f.b.H0, interval=c(1, 20), tol = 1e-8)$root
      }
    }
  }
  
  #Power calculation: Only consider superiority of 1-sided fashion
  power = rep(NA, K); #Power for each analysis
  
  #piecewise power for each weighted logrank test by each analysis using maxcombo rejection boundary
  power.piece = matrix(NA, nrow=K, ncol=max(J)); 
  
  for (i in 1:K){
    if(J[i] >= 2){
      mui = c(mu[i,]); mui = mui[!is.na(mui)]
      ix = (as.numeric(i>1)*sum(J[1:(i-1)])+1) : sum(J[1:i])
      power[i] = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,J[i]),upper=rep(b.Hp[i], J[i]), 
                                      mean=mui, corr = corr.Hp[ix, ix], abseps = 1e-8, maxpts=100000)[1]
    } else {power[i] = 1 - pnorm(b.Hp[i], mean=mu[i, 1])}  
    
    for (j in 1:J[i]){
      power.piece[i, j] = 1-pnorm(b.Hp[i], mean=mu[i, j])
    }
  }
  
  #Overall power and incremental power
  overall.power = power[1]; incr.power = rep(0, K)
  
  if (K > 1){
    if (side == 1){
      for (i in 2:K){
        LL1 = rep(-Inf, sum(J[1:(i-1)]))
        UU1 = rep(b.Hp[1], J[1])
        if(i > 2){for (j in 2:(i-1)){UU1 = c(UU1, rep(b.Hp[j], J[j]))}}
        mu.i1 = mu[1,] 
        if (i > 2) {for (j in 2:(i-1)){mu.i1 = c(mu.i1, mu[j,])}}
        mu.i1 = mu.i1[!is.na(mu.i1)] #mean vector up to (i-1)th analysis
        
        idx1 = 1:sum(J[1:(i-1)])
        if (length(LL1) == 1) {
          P1 = pnorm(UU1, mean=mu.i1)
        } else {
          P1 = mvtnorm::pmvnorm(lower = LL1, upper = UU1, mean=mu.i1, 
                                corr = corr.Hp[idx1, idx1], abseps = 1e-8, maxpts=100000)[1]
        }
        
        LL2 = rep(-Inf, sum(J[1:i]))
        UU2 = rep(b.Hp[1], J[1])
        for (j in 2:i){UU2 = c(UU2, rep(b.Hp[j], J[j]))}
        
        mu.i = mu[1,] 
        for (j in 2:i){mu.i = c(mu.i, mu[j,])}
        mu.i = mu.i[!is.na(mu.i)] #mean vector up to ith analysis
        
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        P2 = mvtnorm::pmvnorm(lower = LL2, upper = UU2, mean=mu.i, 
                              corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
        
        incr.power[i] = P1 - P2
        overall.power = overall.power + incr.power[i]
      }
    }else {
      for (i in 2:K){
        LL = rep(-b.Hp[1], J[1])
        if(i >= 3){for (j in 2:(i-1)){LL = c(LL, rep(-b.Hp[j], J[j]))}}
        LL = c(LL, rep(b.Hp[i], J[i]))
        
        UU = rep(b.Hp[1], J[1])
        if(i >= 3){for (j in 2:(i-1)){UU = c(UU, rep(b.Hp[j], J[j]))}}
        UU = c(UU, rep(Inf, J[i]))
        
        mu.uptoi = mu[1,] 
        for (j in 2:i){mu.uptoi = c(mu.uptoi, mu[j,])}
        mu.uptoi = mu.uptoi[!is.na(mu.uptoi)] #mean vector up to ith analysis
        
        idx = 1:sum(J[1:i]) #indices of the corr matrix
        incr.power[i] = mvtnorm::pmvnorm(lower = LL, upper = UU, mean=mu.uptoi, 
                                         corr = corr.Hp[idx, idx], abseps = 1e-8, maxpts=100000)[1]
        overall.power = overall.power + incr.power[i]
      }
    }
  }
  
  #Calculate the medians
  f.m0 = function(t){S0(t) - 0.5}
  f.m1 = function(t){S1(t) - 0.5}
  
  m0 = uniroot(f.m0, interval= c(1, 1000), tol = 1e-8)$root ##FZ - increase upper bound from 100 to 1000 as discussed. 05/07/2024
  m1 = uniroot(f.m1, interval= c(1, 1000), tol = 1e-8)$root
  
  all.events = data.frame(cbind(events0, events1, events))
  
  maturity0 = events0/(r0*n);
  maturity1 = events1/(r1*n);
  overall.maturity = events/n
  
  maturity = data.frame(cbind(maturity0, maturity1, overall.maturity))
  
  o = list()
  o$mu = mu
  o$events = all.events
  o$maturity = maturity
  
  o$T = T
  o$power = power
  
  #######Critical Values#########
  CV.HR.H0 = CV.HR.H1 = rep(NULL, K)
  if (sum(J) == K){
    lr0 = 0
    for (i in 1:K){for (j in 1:J[i]){lr0 = f.ws[[i]][[j]](rnorm(1)) + lr0}}
    if (lr0 == K){
      #logrank test in all analyses
      
      for (i in 1:K){
        CV.HR.H0[i] = exp(-b.Hp[i]/sqrt(F.entry(T[i])*r0*r1*events[i]))
        r0.H1 = events0[i] / events[i]
        r1.H1 = 1 - r0.H1
        CV.HR.H1[i] = exp(-b.Hp[i]/sqrt(F.entry(T[i])*r0.H1*r1.H1*events[i]))
      }
    }
  }
  
  #only applicable to proportional hazards and using log-rank test
  ph = (f.logHR(1) == f.logHR(10) && f.logHR(10) == f.logHR(100))
  
  if (ph){
    CV.median.H0 = m0/CV.HR.H0
    CV.median.H1 = m0/CV.HR.H1
  } else{CV.median.H0 = CV.median.H1 = NULL}
  
  o$power.piece = power.piece
  o$overall.power = overall.power
  o$incremental.power = incr.power
  o$n = n
  
  o$medians = data.frame(round(cbind(m0, m1),2))
  o$Critical.Values.Medians.H0 = CV.median.H0
  o$Critical.Values.Medians.H1 = CV.median.H1
  
  o$Critical.Values.HR.H0 = CV.HR.H0
  o$Critical.Values.HR.H1 = CV.HR.H1
  
  o$corr.Hp = corr.Hp
  o$corr.strict.H0 = corr.strict.H0
  o$wt = f.ws
  o$bounds = data.frame(cbind(b.Hp, b.H0))
  #H1 setting of distributions;
  setting = list(r=r, alpha=alpha, side=side, h0=h0, S0=S0, h1=h1, S1=S1, log.HR = f.logHR, F.entry=F.entry, G.ltfu=G.ltfu, non.centrality=non.centrality)
  if (show.setting!="N"){o$setting = setting}
  
  #AHR using un-weighted Cox regression
  ahr = rep(NULL, K)
  for (i in 1:K){
    ahr[i] = wlr.AHR(T=T[i], r=r, n = n, h0=h0, S0=S0,h1=h1, S1 = S1,f.logHR = f.logHR,
                     rho = 0, gamma = 0, tau = NULL, s.tau = 0, f.ws=NULL,
                     F.entry = F.entry, G.ltfu = G.ltfu)$AHR
  }
  
  o$AHR = ahr
  return(o)
}