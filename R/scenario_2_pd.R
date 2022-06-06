##' \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with heterogeneous and high-level contamination.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param n_sim number of simulations (large simulations provide more precious estimation).
##' @details \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, high-level contamination), we employed Poisson lognormal distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with heterogeneous and high-level contamination.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,10,15,5,10,10,5)
##' n_sim <- 100000
##' scenario_2_pd(mu, sd, m, n_sim)
##' @usage  scenario_2_pd(mu, sd, m, n_sim)
##' @export
scenario_2_pd <- function(mu, sd = 0.8, m, n_sim){
  k <- length(m)
  rpoislog <- function(S, mu, sig, nu = 1, condS = FALSE, keep0 = FALSE){
    sim <- function(nr) {
      lamx <- rnorm(nr)
      x <- rpois(nr, exp(sig * lamx + mu + log(nu)))
      if (!keep0)
        x <- x[x > 0]
      return(x)
    }
    if (S < 1)
      stop("S is not positive")
    if (!is.finite(S))
      stop("S is not finite")
    if ((S/trunc(S)) != 1)
      stop("S is not an integer")
    if (sig < 0)
      stop("sig is not positive")
    if (nu < 0)
      stop("nu is not positive")
    if (condS) {
      simVec <- vector("numeric", 0)
      fac <- 2
      nr <- S
      while (length(simVec) < S) {
        simvals <- sim(nr * fac)
        simVec <- c(simVec, simvals)
        fac <- (1/(length(simvals)/(nr * fac))) * 2
        fac <- ifelse(is.finite(fac), fac, 1000)
        nr <- S - length(simvals)
      }
      simVec <- simVec[1:S]
    }
    else simVec <- sim(S)
    return(simVec)
  }
  w <- m/sum(m)
  sim1 <- matrix(NA, nrow = n_sim, ncol = k)
  for (j in 1:k) {
    sim1[,j] <-  rpoislog( n_sim, mu, sd, keep0 = TRUE)*w[j]
  }
  sim <- apply(sim1, 1, sum)
  pd <- length(which(sim > 0))/n_sim
  return(pd)
}


