##' \code{\link{scenario_1A_pa}} provides a probability of acceptance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Probability of acceptance estimation when lot with homogeneous contaminations.
##' @param c acceptance number
##' @param lambda expected cell count in the based sample (for this research, we used a quantity of sample which is to be the minimum quantity of selected incremental samples).
##' @param M the quantity (weight) of the aggregate sample.
##' @param m1 the quantity (weight) of the based sample for the risk assessment.
##' @param n number of aggregate samples which are used for inspection.
##' @details \code{\link{scenario_1A_pa}} provides a probability of acceptance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Probability of acceptance under lot with homogeneous contaminations.
##' @examples
##' c <- 0
##' lambda <- 0.05
##' M <- 60
##' m1 <- 5
##' n <- 10
##' scenario_1A_pa(c,lambda, M, m1, n)
##' @usage  scenario_1A_pa(c,lambda, M, m1, n)
##' @export
scenario_1A_pa <- function(c,lambda, M, m1,n){
  pd <- scenario_1A_pd(lambda, M, m1)
  pa <- stats::pbinom(c, n, pd)
  return(pa)
}
