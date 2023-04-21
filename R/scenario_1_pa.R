##' \code{\link{scenario_1_pa}} provides a probability of acceptance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of acceptance estimation when lot with homogeneous contaminations.
##' @param c acceptance number
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param n number of aggregate samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_1_pa}} provides a probability of acceptance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of acceptance under lot with homogeneous contaminations.
##' @examples
##' c <- 0
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' n <- 10
##' scenario_1_pa(c, mu, sd = 0.8, m, n, type = "theory")
##' scenario_1_pa(c, mu, sd = 0.8, m, n, type = "simulation", n_sim = 1000000)
##' @usage  scenario_1_pa(c, mu, sd, m, n, type, n_sim)
##' @export
scenario_1_pa <- function(c, mu, sd = 0.8, m, n, type, n_sim = NA) {
  if (type == "theory") {
    pd <- scenario_1_pd(mu, sd, m, type, n_sim)
    pa <- stats::pbinom(c, n, pd)
    return(pa)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      pd <- scenario_1_pd(mu, sd, m, type, n_sim)
      pa <- stats::pbinom(c, n, pd)
      return(pa)
    }
  } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}
