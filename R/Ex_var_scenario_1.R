##' \code{\link{Ex_var_scenario_1}} provides the expected value or variance of the number of microorganisms in the aggregate sample under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Expected value or variance estimation of microorganism in the aggregate sample under scenario 1.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param measure what type of measure you would like to consider for the graph, such as "expectation" or "variance".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{Ex_var_scenario_1}} provides the expected value or variance of the number of microorganisms in the aggregate sample under scenario 1 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return expected value or variance of the number of microorganisms in the aggregate sample.
##' @examples
##' mu <- 2
##' sd <- 0.8
##' m <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' Ex_var_scenario_1(mu, sd = 0.8, m, type = "theory",  measure = "expectation")
##' Ex_var_scenario_1(mu, sd = 0.8, m, type = "simulation",  measure = "expectation", n_sim = 2000000)
##' Ex_var_scenario_1(mu, sd = 0.8, m, type = "theory",  measure = "variance")
##' Ex_var_scenario_1(mu, sd = 0.8, m, type = "simulation",  measure = "variance", n_sim = 2000000)
##' @usage  Ex_var_scenario_1(mu, sd = 0.8, m, type, measure, n_sim)
##' @export
Ex_var_scenario_1 <- function(mu, sd = 0.8, m, type, measure, n_sim = NA){
  k <- length(m)
  w <- m/sum(m)
  lambda_0 <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  if (measure == "expectation") {
    if (type == "theory") {
      expect <- lambda_0*sum(m*m/(min(m)*sum(m)))
      return(expect)
    } else if (type == "simulation") {
      if (is.na(n_sim) == TRUE) {
        stop("please set the number of simualtions")
      } else {
        sim1 <- matrix(NA, nrow = n_sim, ncol = k)
        for (j in 1:k) {
          sim1[,j] <- stats::rpois(n_sim, lambda_0*m[j]/min(m))
        }
        sim3 <- apply(sim1, 2, mean)
        sim4 <- matrix(NA, nrow = 1, ncol = k)
        for (j in 1:k) {
          sim4[,j] <- sim3[j]*w[j]
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
        Vari <- lambda_0*sum(m*m*m/(min(m)*sum(m)*sum(m)))
        return(Vari)
      } else if (type == "simulation") {
        if (is.na(n_sim) == TRUE) {
          stop("please set the number of simualtions")
        } else {
          sim1 <- matrix(NA, nrow = n_sim, ncol = k)
          for (j in 1:k) {
            sim1[,j] <- stats::rpois(n_sim, lambda_0*m[j]/min(m))
          }
          sim3 <- apply(sim1, 2, var)
          sim4 <- matrix(NA, nrow = 1, ncol = k)
          for (j in 1:k) {
            sim4[,j] <- sim3[j]*w[j]*w[j]
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


