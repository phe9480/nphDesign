#'
#' Summary of Study Design Assumptions
#' 
#' Display the study design assumptions in a nice table format. Require R packages gt and dplyr
#'
#' @param  dist Distribution specification for each arm
#' @param  accr Accrual curve specification for the study
#' @param alphabeta Type I and type II errors
#' @param test The specification of statistical test
#' 
#' @return A nice formatted summary table for study deisgn assumptions
#' 
#' 
#' @export
#' 
#Display design parameters
display.stat.setting = function(dist = PFS.dist, accr = study.accr, 
                                alphabeta = PFS.alphabeta, test = PFS.test,
                                caption="Summary of Accrual and Statistical Setting"){
  d = data.frame(cbind(dist))
  d[row.names(d)=="control", 1] = paste(dist$control$dist, "distribution with median", dist$control$median)
  if (dist$exp$delay==0){exp.d = paste("Exponential distribution with HR", dist$exp$HR)} else{exp.d = paste("Delayed effect: ", dist$exp$delay, "and the HR after delay:", dist$exp$HR)}
  d[row.names(d)=="exp", 1] = exp.d
  
  a = data.frame(cbind(accr))
  ab = data.frame(cbind(alphabeta))
  t = data.frame(cbind(test))
  
  d2 = rename(d, Description = dist)
  a2 = rename(a, Description = accr)
  ab2 = rename(ab, Description = alphabeta)
  t2 = rename(t, Description = test)
  
  d2$stat = "Statistical Distribution"
  a2$stat = "Assumed Study Accrual"
  ab2$stat = "Type I (Alpha) and Type II Error (Beta)"
  t2$stat = "Statistical Test"
  
  out = data.frame(rbind(a2, d2, ab2, t2))
  out$label = row.names(out)
  
  tab = out  %>% group_by(stat)  %>% gt()
  tab = tab  %>% cols_move(columns = c(Description), after = label)
  tab = tab %>% tab_header(
    title = md(caption), 
  ) %>% tab_footnote(
    footnote = "The unit for median and accrual is month.",
    locations = cells_column_labels(columns = Description)
  )
  return(tab)
}
