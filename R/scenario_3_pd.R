##' \code{\link{scenario_3_pd}} provides a probability of detection under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with heterogeneous and low-level contamination.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples(with equal/unequal weights).
##' @param K shape parameter (default value 0.25).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_3_pd}} provides a probability of detection under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, low-level contamination), we employed Poisson gamma distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with heterogeneous and low-level contamination by using theoretical or simulations based results.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' K <- 0.05
##' scenario_3_pd(mu, sd = 0.8, m, K, type = "theory")
##' scenario_3_pd(mu, sd = 0.8, m, K, type = "simulation", n_sim = 2000000)
##' @usage  scenario_3_pd(mu, sd, m, K, type, n_sim)
##' @export
scenario_3_pd <- function(mu, sd = 0.8, m, K = 0.25, type, n_sim = NA){
    if (type == "theory") {
      lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
      k <- length(m)
      sim1 <- matrix(NA, nrow = k, ncol = 1)
      for (i in 1:k) {
        # sim1[i,1] <-   ((m[i]*lambda_0)/(min(m) + m[i]*lambda_0))^K
        # sim1[i,1] <-   (K*min(m)/(K*min(m) + m[i]*lambda_0))^K
        sim1[i,1] <-   (min(m)/(min(m) + m[i]*lambda_0))^K
      }
      sim <- apply(sim1, 2, prod)
      pd <- 1 - sim
      return(pd)
    } else if (type == "simulation") {
      if (is.na(n_sim) == TRUE) {
        stop("please set the number of simualtions")
      } else {
        k <- length(m)
        lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        w <- m/sum(m)
        sim1 <- matrix(NA, nrow = n_sim, ncol = k)
        for (j in 1:k) {
          # sim1[,j] <-  extraDistr::rgpois(n_sim, shape = K, rate = lambda)*w[j]
          sim1[,j] <-  extraDistr::rgpois(n_sim, shape = K, rate =  1/(m[j]*lambda_0/min(m)))*w[j]
        }
        sim <- apply(sim1, 1, sum)
        pd <- length(which(sim > 0))/n_sim
        return(pd)
      }
    }else {
      print("please include what type (theory/ simulation) you would like to consider")
    }
  }

