##' \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with heterogeneous contamination and concentration levels fluctuating from sublot to sublot.
##' @param mu the the mean concentration (\eqn{\mu}).
##' @param sd_b standard deviation of concentration level between sublots on the log10 scale (default value 0.8).
##' @param sd_w standard deviation of concentration level within sublot on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_2_pd}} provides a probability of detection under scenario 2 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with heterogeneous, high-level contamination), we employed Poisson lognormal distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection when lot with heterogeneous and high-level contamination by using theoretical or simulation-based results.
##' @examples
##' m1 <- c(5,5,5,5,5,5,5,5,5,5)
##' m2 <- c(5,5,5,5,5,5,5,5,5,5)
##' m3 <- c(5,5,5,5,5,5,5,5,5,5)
##' m4 <- c(5,5,5,5,5,5,5,5,5,5)
##' mu <- -3
##' m <- list(m1,m2,m3,m4)
##' sd_b <- 0.2
##' sd_w <- 0.8
##' scenario_5_pd(mu, sd_b, sd_w = 0.8, m, type = "theory")
##' scenario_5_pd(mu, sd_b, sd_w = 0.8, m, type = "simulation", n_sim = 20000)
##' @usage  scenario_5_pd(mu, sd_b, sd_w, m, type, n_sim)
##' @export
scenario_5_pd <- function(mu, sd_b, sd_w = 0.8, m, type, n_sim = NA) {
  # warning("\033[1;31m","Please define m as a list form of each set of incremental samples from sublots, for example m <- list(m1,m2,...)", call. = FALSE)
  set.seed(1000)
  mu1 <- rnorm(length(m), mu, sd_b)
  sim2 <- matrix(NA, nrow = 1, ncol = length(m))
  for (j in 1:length(m)) {
    sim2[, j] <- 1 - scenario_2_pd(mu1[j], sd_w, unlist(m[j]), type, n_sim)
  }
  p_d <- 1 - prod(sim2)
  return(p_d)
}
