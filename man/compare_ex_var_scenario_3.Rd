% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_ex_var_scenario_3.R
\name{compare_ex_var_scenario_3}
\alias{compare_ex_var_scenario_3}
\title{Graphical displays of expectation or variance of different sampling schemes under scenario 3.}
\usage{
compare_ex_var_scenario_3(mulow, muhigh, sd, m1, m2, K, measure)
}
\arguments{
\item{mulow}{the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.}

\item{muhigh}{the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.}

\item{sd}{standard deviation on the log10 scale (default value 0.8).}

\item{m1}{the vector of the first set of incremental samples (with equal/unequal weights).}

\item{m2}{the vector of the second set of incremental samples (with equal/unequal weights).}

\item{K}{shape parameter (default value 0.25).}

\item{measure}{what type of measure you would like to consider for the graph, such as "expectation" or "variance".}
}
\value{
Graphical displays based on expectation or variance when lot with heterogeneous and low-level contamination.
}
\description{
\code{\link{compare_ex_var_scenario_3}} provides graphical displays based on expectation or variance under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
}
\details{
\code{\link{compare_ex_var_scenario_3}} provides a probability of acceptance under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
Under this scenario (a lot with homogeneous contaminations), we employed Poisson gamma distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
}
\examples{
mulow <- 0
muhigh <- 2
m1 <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
10,10,10,10,10,10,10,10,10,10)
m2 <- c(15,5,5,5,10,5,10,5,15,10,5,10,5,25,10,5,10,5,5,10,5,15,10,
5,5,20,5,10,5,10,20,5,10,30,5,20,5,10,5,10,20,15,10,15,10,
10,5,10,15,5)
K <- 0.05
compare_ex_var_scenario_3(mulow, muhigh, sd = 0.8, m1, m2, K, measure = "variance")
compare_ex_var_scenario_3(mulow, muhigh, sd = 0.8, m1, m2, K, measure =  "expectation")
}
