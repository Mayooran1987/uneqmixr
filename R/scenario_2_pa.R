##' \code{\link{scenario_2_pa}} provides a probability of acceptance under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of acceptance estimation when lot with heterogeneous and high-level contamination.
##' @param c acceptance number
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param n number of aggregate samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_2_pa}} provides a probability of acceptance under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, high-level contamination), we employed Poisson lognormal distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of acceptance when lot with heterogeneous and high-level contamination.
##' @examples
##' c <- 0
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,10,15,5,10,10,5)
##' n <- 10
##' n_sim <- 100000
##' scenario_2_pa(c, mu, sd, m, n, n_sim)
##' @usage  scenario_2_pa(c, mu, sd, m, n, n_sim)
##' @export
scenario_2_pa <- function(c, mu, sd = 0.8, m, n, n_sim){
  pd <- scenario_2_pd(mu, sd, m, n_sim)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}

