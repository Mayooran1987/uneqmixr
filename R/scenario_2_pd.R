##' \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with heterogeneous and high-level contamination.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, high-level contamination), we employed Poisson lognormal distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with heterogeneous and high-level contamination by using theoretical or simulation-based results.
##' @examples
##' mu <- -3
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' scenario_2_pd(mu, sd = 0.8, m, type = "theory")
##' scenario_2_pd(mu, sd = 0.8, m, type = "simulation", n_sim = 2000000)
##' @usage  scenario_2_pd(mu, sd, m, type, n_sim)
##' @export
scenario_2_pd <- function(mu, sd = 0.8, m, type, n_sim = NA){
  if (type == "theory") {
    # print("Not yet established, please use simulation-based results")
    prob <- matrix(NA, nrow = length(m), ncol = 1)
    for (j in 1:length(m)) {
      # some codes used from "poilog" package
      dpoilog <- function(n,mu,sig){
        if (length(mu) > 1 | length(sig) > 1) stop('vectorization of mu and sig is currently not implemented')
        if (any((n[n != 0]/trunc(n[n != 0])) != 1)) stop('all n must be integers')
        if (any(n < 0)) stop('one or several values of n are negative')
        if (!all(is.finite(c(mu,sig)))) stop('all parameters should be finite')
        if (sig <= 0) stop('sig is not larger than 0')
        .C('poilog1',n = as.integer(n),mu = as.double(mu),sig2 = as.double(sig^2),nrN = as.integer(length(n)),val = double(length(n)))$val
      }
      prob[j,] <- dpoilog( 0, log10(m[j]/min(m)) + mu, sd)
    }
    pd <- 1 - apply(prob, 2, prod)
    return(pd)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
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
        # sim1[,j] <-  rpoislog( n_sim, m[j]*mu/min(m), sd, keep0 = TRUE)*w[j]
        sim1[,j] <-  rpoislog( n_sim, (log((m[j]/min(m)),10) + mu), sd, keep0 = TRUE)*w[j]
        # if X follows PLN(mu, sigma) then cX follows PLN(log(c)+mu, sigma)); where c is any positive constant.
      }
      sim <- apply(sim1, 1, sum)
      pd <- length(which(sim > 0))/n_sim
      return(pd)
    }
  }else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}

