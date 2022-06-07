##' \code{\link{scenario_1B_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with homogeneous contaminations based on simulations.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_1B_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with homogeneous contaminants by simulations results.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,10,15,5,10,10,5)
##' n_sim <- 100000
##' scenario_1B_pd(mu, sd, m, n_sim)
##' @usage  scenario_1B_pd(mu, sd, m, n_sim)
##' @export
scenario_1B_pd <- function(mu, sd = 0.8, m, n_sim){
  k <- length(m)
  w <- m/sum(m)
  lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  sim1 <- matrix(NA, nrow = n_sim, ncol = k)
  for (j in 1:k) {
    sim1[,j] <-   stats::rpois( n_sim, lambda)*w[j]
  }
  sim <- apply(sim1, 1, sum)
  pd <- length(which(sim > 0))/n_sim
  return(pd)
}


