##' \code{\link{scenario_3_pa}} provides a probability of acceptance  under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of acceptance  estimation when lot with heterogeneous and low-level contamination.
##' @param c acceptance number
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param K shape parameter (default value 0.25).
##' @param n number of aggregate samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation".But,not yet established for this scenario at this time.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_3_pa}} provides a probability of acceptance  under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, low-level contamination), we employed Poisson gamma distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of acceptance  when lot with heterogeneous and low-level contamination.
##' @examples
##' c <- 0
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' K <- 0.25
##' n <- 10
##' scenario_3_pa(c, mu, sd = 0.8, m, K, n, type = "theory")
##' scenario_3_pa(c, mu, sd = 0.8, m, K, n, type = "simulation", n_sim = 2000000)
##' @usage  scenario_3_pa(c, mu, sd, m, K,  n, type, n_sim)
##' @export
scenario_3_pa <- function(c, mu, sd = 0.8, m, K, n, type,  n_sim = NA){
  if (type == "theory") {
     pd <- scenario_3_pd(mu, sd, m, K, type)
     pa <- stats::pbinom(c, n, pd)
     return(pa)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      pd <- scenario_3_pd(mu, sd, m, K, type, n_sim)}
    pa <- stats::pbinom(c, n, pd)
    return(pa)
  } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}

