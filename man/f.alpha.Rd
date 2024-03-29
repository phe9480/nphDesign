% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f.alpha.R
\name{f.alpha}
\alias{f.alpha}
\title{Calculation of Incremental Alpha for Group Sequential Design}
\usage{
f.alpha(overall.alpha = 0.025, side = 1, sf = "LDOF", timing = c(0.75, 1))
}
\arguments{
\item{overall.alpha}{Allocated overall alpha for the GSD}

\item{side}{Side of test}

\item{sf}{Spending function. "LDOF", "LDPK"}

\item{timing}{Timing of IA and FA}
}
\description{
The function includes the incremental alpha for group sequential design.
The incremental alpha is a vector with sum to the overall alpha. It is may be needed for
calling the group sequential design power functions.
}
