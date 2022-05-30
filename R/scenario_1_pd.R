##' \code{\link{scenario_1_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of detection estimation when lot with homogeneous contaminations.
##' @param lambda expected cell count in the based sample (for this research, we used a quantity of sample which is to be the minimum quantity of selected incremental samples).
##' @param m the quantity (weight) of the aggregate sample.
##' @param m1 the quantity (weight) of the based sample for the risk assessment.
##' @details \code{\link{scenario_1_pd}} provides a probability of detection under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of detection under lot with homogeneous contaminations.
##' @examples
##' lambda <- 0.05
##' m <- 25
##' m1 <- 5
##' scenario_1_pd(lambda, m, m1)
##' @usage  scenario_1_pd(lambda, m, m1)
##' @export
scenario_1_pd <- function(lambda, m, m1){
  pd <- 1 - exp(lambda*(-m/m1))
  return(pd)
}
