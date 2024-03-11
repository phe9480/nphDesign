#' Display of Simulation Results Produced From simulation.nphDesign Function
#' 
#' This function plots the density curve of z values and power curve by analysis
#' 
#' @param s.nphDesign An object produced from simulation.nphDesign Function
#'
#' @return Display of the graphs
#'  
#' 
#' @export
#' 
display.simulation = function(s.nphDesign){

  K = length(s.nphDesign$power)
  nSim = dim(s.nphDesign$wlr.simulations)[1]
  
  #Tabulate power
  colnm = NULL
  for (i in 1:K){
    ai = paste("Analysis", i)
    colnm = c(colnm, ai)
  }
    
  knitr::kable(round(s.nphDesign$power*100,2), col.names = colnm, 
               caption="Power for Each Analysis (%)", 
               format = "markdown")
  knitr::kable(round(s.nphDesign$overall.power*100,2), col.names ="Overall Power", 
               caption="Overall Study Power", format = "markdown")
  
  #Power curve over analysis
  pow = s.nphDesign$power; 
  overall.pow = s.nphDesign$overall.power
  plot(1:K, pow, ylim=c(0, 1), 
       type="n", xlim=c(0, K+1), xlab="Analysis", ylab="Power",
       main=paste("Power by simulations: ", nSim, "trials"))
  for  (i in 1:K){
    points(i, pow[i], pch=i, col=i, lwd=3, cex=2)
  }  
  points(K+1, overall.pow, pch=3, col="maroon", lwd=4, cex=2)
  
  legend(0, 0.3, c("Overall Power"), pch=3, col="maroon", bty="n")
  
  #Density curve of z values
  for (i in 1:K){
    fit = density(s.nphDesign$wlr.simulations[,1,i,1])
    plot(fit, ylab="z Values", 
         main = paste("z based on", nSim, "simulations for Analysis", i))
    abline(v=s.nphDesign$wlr.simulations[1,1,i,4], col="green", lwd=3)
    legend(min(fit$x), max(fit$y), c("Bound"), bty="n", col="green", lwd=3, lty=1, cex=0.8)
  }
  
  #If logrank test is requested, then additionally display
  if(!is.null(s.nphDesign$lr.power)){
    lr.pow = s.nphDesign$lr.power; 
    lr.overall.pow = s.nphDesign$lr.overall.power
    plot(1:K, lr.pow, ylim=c(0, 1), 
         type="n", xlim=c(0, K+1), xlab="Analysis", ylab="Power",
         main=paste("Logrank Power by simulations: ", nSim, "trials"))
    for  (i in 1:K){
      points(i, lr.pow[i], pch=i, col=i, lwd=3, cex=2)
    }  
    points(K+1, lr.overall.pow, pch=3, col="maroon", lwd=4, cex=2)
    
    legend(0, 0.3, c("Overall Logrank Power"), pch=3, col="maroon", bty="n")
    
    #Density curve of z values
    for (i in 1:K){
      fit = density(s.nphDesign$lr.simulations[,i,1])
      plot(fit, ylab="z Values", 
           main = paste("Logrank z based on", nSim, "simulations for Analysis", i))
      abline(v=s.nphDesign$lr.simulations[1,i,4], col="green", lwd=3)
      legend(min(fit$x), max(fit$y), c("Bound"), bty="n", col="green", lwd=3, lty=1, cex=0.8)
    }
    
  }
}

