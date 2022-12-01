##' \code{\link{scenario_5_prevalence}} provides a prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study.
##' @title Prevalence estimation before inspection when lot with heterogeneous contamination and contamination levels fluctuate from lot to lot.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param l the number of lots in the production process.
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_5_prevalence}} provides a prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return Prevalence estimation before inspection when lot with heterogeneous contamination and contamination levels fluctuate from lot to lot by using theoretical or simulation-based results.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' l <- 5000
##' scenario_5_prevalence(mu, sd, m, l, type = "theory")
##' @usage  scenario_5_prevalence(mu, sd, m, l, type, n_sim)
##' @export
scenario_5_prevalence <- function(mu, sd =0.8, m, l, type, n_sim = NA ){
  if (type == "theory") {
    # print("Not yet established, please use simulation-based results")
    mu_0 <- stats::rnorm(l, mu, sd)
    sim2 <- matrix(NA, nrow = 1, ncol = l)
    for (j in 1:l) {
      sim2[,j] <-   scenario_2_pd(mu_0[j], sd, m, type = "theory")
    }
    prev <- sum(sim2)/l
    return(prev)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      mu_0 <- stats::rnorm(l, mu, sd)
      sim2 <- matrix(NA, nrow = 1, ncol = l)
      for (j in 1:l) {
        sim2[,j] <-   scenario_2_pd(mu_0[j], sd, m, type = "simulation", n_sim)
      }
      prev <- sum(sim2)/l
      return(prev)
    }
  }else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}

