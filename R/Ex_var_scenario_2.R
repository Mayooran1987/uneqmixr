##' \code{\link{Ex_var_scenario_2}} provides the expected value or variance of the number of microorganisms in the aggregate sample under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' @title Expected value or variance estimation of microorganism in the aggregate sample under scenario 2.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param measure what type of measure you would like to consider for the graph, such as "expectation" or "variance".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{Ex_var_scenario_2}} provides the expected value or variance of the number of microorganisms in the aggregate sample under scenario 2 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return expected value or variance of the number of microorganisms in the aggregate sample.
##' @examples
##' mu <- 2
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' Ex_var_scenario_2(mu, sd = 0.8, m, type = "theory",  measure = "expectation")
##' Ex_var_scenario_2(mu, sd = 0.8, m, type = "simulation",  measure = "expectation", n_sim = 2000000)
##' Ex_var_scenario_2(mu, sd = 0.8, m, type = "theory",  measure = "variance")
##' Ex_var_scenario_2(mu, sd = 0.8, m, type = "simulation",  measure = "variance", n_sim = 2000000)
##' @usage  Ex_var_scenario_2(mu, sd = 0.8, m, type, measure, n_sim)
##' @export
Ex_var_scenario_2 <- function(mu, sd = 0.8, m, type, measure, n_sim = NA) {
  # lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  if (measure == "expectation") {
    if (type == "theory") {
      expect <- (exp(0.5 * sd^2) / sum(m)) * sum(m * exp(log((m / min(m)), 10) + mu))
      return(expect)
    } else if (type == "simulation") {
      if (is.na(n_sim) == TRUE) {
        stop("please set the number of simualtions")
      } else {
        k <- length(m)
        rpoislog <- function(S, mu, sig, nu = 1, condS = FALSE, keep0 = FALSE) {
          sim <- function(nr) {
            lamx <- rnorm(nr)
            x <- rpois(nr, exp(sig * lamx + mu + log(nu)))
            if (!keep0) {
              x <- x[x > 0]
            }
            return(x)
          }
          if (S < 1) {
            stop("S is not positive")
          }
          if (!is.finite(S)) {
            stop("S is not finite")
          }
          if ((S / trunc(S)) != 1) {
            stop("S is not an integer")
          }
          if (sig < 0) {
            stop("sig is not positive")
          }
          if (nu < 0) {
            stop("nu is not positive")
          }
          if (condS) {
            simVec <- vector("numeric", 0)
            fac <- 2
            nr <- S
            while (length(simVec) < S) {
              simvals <- sim(nr * fac)
              simVec <- c(simVec, simvals)
              fac <- (1 / (length(simvals) / (nr * fac))) * 2
              fac <- ifelse(is.finite(fac), fac, 1000)
              nr <- S - length(simvals)
            }
            simVec <- simVec[1:S]
          } else {
            simVec <- sim(S)
          }
          return(simVec)
        }
        w <- m / sum(m)
        sim1 <- matrix(NA, nrow = n_sim, ncol = k)
        for (j in 1:k) {
          sim1[, j] <- rpoislog(n_sim, (log((m[j] / min(m)), 10) + mu), sd, keep0 = TRUE)
          # if X follows PLN(mu, sigma) then cX follows PLN(log(c)+mu, sigma)); where c is any positive constant.
        }
        sim3 <- apply(sim1, 2, mean)
        sim4 <- matrix(NA, nrow = 1, ncol = k)
        for (j in 1:k) {
          sim4[, j] <- sim3[j] * w[j]
        }
        sim <- apply(sim4, 1, sum)
        result <- sum(sim)
        return(result)
      }
    } else {
      print("please include what type (theory/ simulation) you would like to consider")
    }
  } else if (measure == "variance") {
    if (type == "theory") {
      var <- (exp(0.5 * sd^2) / (sum(m) * sum(m))) * (sum(m * m * exp(log((m / min(m)), 10) + mu))) + (((exp(sd^2) - 1) * exp(sd^2)) / (sum(m) * sum(m))) * sum(m * m * exp(2 * (log((m / min(m)), 10) + mu)))
      return(var)
    } else if (type == "simulation") {
      if (is.na(n_sim) == TRUE) {
        stop("please set the number of simualtions")
      } else {
        k <- length(m)
        rpoislog <- function(S, mu, sig, nu = 1, condS = FALSE, keep0 = FALSE) {
          sim <- function(nr) {
            lamx <- rnorm(nr)
            x <- rpois(nr, exp(sig * lamx + mu + log(nu)))
            if (!keep0) {
              x <- x[x > 0]
            }
            return(x)
          }
          if (S < 1) {
            stop("S is not positive")
          }
          if (!is.finite(S)) {
            stop("S is not finite")
          }
          if ((S / trunc(S)) != 1) {
            stop("S is not an integer")
          }
          if (sig < 0) {
            stop("sig is not positive")
          }
          if (nu < 0) {
            stop("nu is not positive")
          }
          if (condS) {
            simVec <- vector("numeric", 0)
            fac <- 2
            nr <- S
            while (length(simVec) < S) {
              simvals <- sim(nr * fac)
              simVec <- c(simVec, simvals)
              fac <- (1 / (length(simvals) / (nr * fac))) * 2
              fac <- ifelse(is.finite(fac), fac, 1000)
              nr <- S - length(simvals)
            }
            simVec <- simVec[1:S]
          } else {
            simVec <- sim(S)
          }
          return(simVec)
        }
        w <- m / sum(m)
        sim1 <- matrix(NA, nrow = n_sim, ncol = k)
        for (j in 1:k) {
          sim1[, j] <- rpoislog(n_sim, (log((m[j] / min(m)), 10) + mu), sd, keep0 = TRUE)
          # if X follows PLN(mu, sigma) then cX follows PLN(log(c)+mu, sigma)); where c is any positive constant.
        }
        sim3 <- apply(sim1, 2, var)
        sim4 <- matrix(NA, nrow = 1, ncol = k)
        for (j in 1:k) {
          sim4[, j] <- sim3[j] * w[j] * w[j]
        }
        sim <- apply(sim4, 1, sum)
        result <- sim
        return(result)
      }
    } else {
      print("please include what type (theory/ simulation) you would like to consider")
    }
  } else {
    print("Please choose measure as expectation or variance")
  }
  return(result)
}
