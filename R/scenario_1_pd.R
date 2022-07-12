##' \code{\link{scenario_1_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with homogeneous contaminations.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_1_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with homogeneous contaminants by using theoretical or simulation-based results.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' scenario_1_pd(mu, sd = 0.8, m, type = "theory")
##' scenario_1_pd(mu, sd = 0.8, m, type = "simulation", n_sim = 1000000)
##' @usage  scenario_1_pd(mu, sd, m, type, n_sim)
##' @export
scenario_1_pd <- function(mu, sd = 0.8, m, type, n_sim = NA){
  if (type == "theory") {
    lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
    M <- sum(m)
    m_0 <- min(m)
    pd <- 1 - exp(lambda_0*(-M/m_0))
    return(pd)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      k <- length(m)
      w <- m/sum(m)
      lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
      sim1 <- matrix(NA, nrow = n_sim, ncol = k)
      for (j in 1:k) {
        sim1[,j] <-   stats::rpois( n_sim, lambda_0*m[j]/min(m))*w[j]
      }
      sim <- apply(sim1, 1, sum)
      pd <- length(which(sim > 0))/n_sim
      # warning("Please note that you can get more accurate results if you use a large number of simulations")
      return(pd)
    }
  }else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}

