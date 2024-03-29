% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scenario_5_prevalence.R
\name{scenario_5_prevalence}
\alias{scenario_5_prevalence}
\title{Prevalence estimation before inspection when lot with heterogeneous contamination and contamination levels fluctuate from lot to lot.}
\usage{
scenario_5_prevalence(mu, sd, m, l, type, n_sim)
}
\arguments{
\item{mu}{the the mean concentration (\eqn{\mu}).}

\item{sd}{standard deviation on the log10 scale (default value 0.8).}

\item{m}{the vector of incremental samples (with equal/unequal weights).}

\item{l}{the number of lots in the production process.}

\item{type}{what type of the results you would like to consider such as "theory" or "simulation".}

\item{n_sim}{number of simulations (large simulations provide more precise estimation).}
}
\value{
Prevalence estimation before inspection when lot with heterogeneous contamination and contamination levels fluctuate from lot to lot by using theoretical or simulation-based results.
}
\description{
\code{\link{scenario_5_prevalence}} provides a prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study.
}
\details{
\code{\link{scenario_5_prevalence}} provides a prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
}
\examples{
mu <- -3
sd <- 0.8
m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
l <- 5000
scenario_5_prevalence(mu, sd, m, l, type = "theory")
}
