% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.stat.setting.R
\name{display.stat.setting}
\alias{display.stat.setting}
\title{Summary of Study Design Assumptions}
\usage{
display.stat.setting(
  dist = PFS.dist,
  accr = study.accr,
  alphabeta = PFS.alphabeta,
  test = PFS.test,
  caption = "Summary of Accrual and Statistical Setting"
)
}
\arguments{
\item{dist}{Distribution specification for each arm}

\item{accr}{Accrual curve specification for the study}

\item{alphabeta}{Type I and type II errors}

\item{test}{The specification of statistical test}
}
\value{
A nice formatted summary table for study deisgn assumptions
}
\description{
Display the study design assumptions in a nice table format. Require R packages gt and dplyr
}
