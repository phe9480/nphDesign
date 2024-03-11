#' Generate a dataset by administrative data cut at a pre-specified number of target events
#'
#' Cut data by according to calendar time when number of Target Events is reached.
#' 
#' @param data Dataset must include variables: 
#'  \describe{
#'    \item{enterTime}{The start time for survival calculation, eg, randomization time}
#'    \item{calendarTime}{The time when event/censoring occurred}
#'    \item{survTime}{Survival time for analysis}
#'    \item{cnsr}{Indicator of event(0) or censor (1)}
#'    }
#' @param targetEvents The number of target events for analysis.
#' @return A dataframe including variables:
#'    calendarCutoff
#'    survTimeCut
#'    cnsrCut
#'    
#' @export
f.dataCut = function(data, targetEvents = 397) {
  data0 = data
  data0.order <- data0[order(data0$calendarTime), ] #order by calendar time
  data.event <- data0.order[data0.order$cnsr == 0,] #Events Only
  
  data.event$event.seq <- seq.int(nrow(data.event)) #event sequence number
  data0$calendarCutoff = data.event$calendarTime[data.event$event.seq == targetEvents] #Data cutoff in calendar time added to the original dataframe as a variable
  data0$survTimeCut = ifelse(data0$calendarTime <= data0$calendarCutoff, data0$survTime, data0$calendarCutoff-data0$enterTime)
  
  data0$cnsrCut = ifelse(data0$calendarTime <= data0$calendarCutoff, data0$cnsr, 1)
  
  return(data0)
}

