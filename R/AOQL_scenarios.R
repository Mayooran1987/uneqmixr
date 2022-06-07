##' \code{\link{AOQL_scenarios}} provides the Average Outgoing Quality (AOQ) curve and calculates Average Outgoing Quality Level (AOQL) value based on expected microbial counts in each scenario.
##' @title Construction of  AOQ curve and calculate AOQL value based on average microbial counts
##' @param c acceptance number
##' @param llim the upper limit for graphing the arithmetic mean of cell count
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m the vector of incremental samples (with equal/unequal weights).
##' @param n number of aggregate samples which are used for inspection.
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @param scenario what scenario we have considered such as  \code{"1B"}or \code{"2"}or \code{"3"}
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @details  Since \eqn{P_a} is the probability of acceptance, \eqn{\lambda} is the arithmetic mean of cell count and the outgoing contaminated arithmetic mean of cell count of incremental samples is given by \eqn{AOQ} as the product \eqn{\lambda P_a}.
##'           The quantity \eqn{AOQL} is defined as the maximum proportion of outgoing contaminated incremental samples and is given by \deqn{AOQL ={\max_{\lambda \geq 0}}{\lambda P_a}}
##' @return AOQ curve and AOQL value based on expected microbial counts in each scenario.
##' @seealso  \link{scenario_1B_pa}, \link{scenario_2_pa},\link{scenario_3_pa}
##' @examples
##' c <- 0
##' llim <- 0.5
##' sd <- 0.8
##' m1 <- c(10,10,10,10,10,10)
##' m2 <- c(10,12,18,20)
##' n <- 10
##' scenario <- "1B"
##' n_sim <- 100000
##' AOQL_scenarios(c,llim, sd, m=m1, n, scenario, K, n_sim)
##' AOQL_scenarios(c,llim, sd, m=m2, n, scenario, K, n_sim)
##' @usage  AOQL_scenarios(c,llim, sd, m, n, scenario, K, n_sim)
##' @export
AOQL_scenarios <- function(c,llim, sd, m, n, scenario,K =0.25, n_sim){
k <- length(m)
M <- sum(m)
Sampling_scheme <- NULL
f_spr <- function(k,M) {
  sprintf("Scheme(k=%.0f, M=%.0f)", k, M)
}
lambda <- seq(0, llim, by = 1e-04)
mu <- log(lambda, 10) - (sd^2/2)*log(10, exp(1))

# AOQ <- matrix(NA, nrow = length(mu), ncol = 1)
# for (i in 1:length(mu)) {
#   AOQ[i,1] <-  lambda[i]*(1 - scenario_1B_pd(mu[i],sd,m1,n_sim))^n
#   # AOQ[i,2] <-  lambda[i]*(1 - scenario_1B_pd(mu[i],sd,m2,n_sim))^n
# }
if (scenario == "1B") {
  AOQ <- matrix(NA, nrow = length(mu), ncol = 1)
  for (i in 1:length(mu)) {
    # AOQ[i,1] <-  lambda[i]*(1 - scenario_1B_pd(mu[i],sd,m,n_sim))^n
    AOQ[i,1] <-  lambda[i]*scenario_1B_pa(c, mu[i], sd, m, n, n_sim)
  }
} else if (scenario == "2") {
  AOQ <- matrix(NA, nrow = length(mu), ncol = 1)
  for (i in 1:length(mu)) {
    # AOQ[i,1] <-  lambda[i]*(1 - scenario_2_pd(mu[i],sd,m,n_sim))^n
    AOQ[i,1] <-  lambda[i]*scenario_2_pa(c, mu[i], sd, m, n, n_sim)
  }
} else if (scenario == "3") {
  AOQ <- matrix(NA, nrow = length(mu), ncol = 1)
  for (i in 1:length(mu)) {
    # AOQ[i,1] <-  lambda[i]*(1 - scenario_3_pd(mu[i],sd,m,K = 0.25,n_sim))^n
    AOQ[i,1] <-  lambda[i]*scenario_3_pa(c, mu[i], sd, m, n, K, n_sim)
  }
} else {
  print("please include what scenario you would like to consider")
}

Prob <- data.frame(lambda, AOQ )
colnames(Prob ) <- c("lambda", f_spr(k,M))

melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "AOQ")
ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = AOQ, group = Sampling_scheme, colour = Sampling_scheme)) +
  ggplot2::theme_classic() + ggplot2::ylab(expression(AOQ)) + ggthemes::scale_colour_colorblind() +
  ggplot2::xlab(expression("arithmetic mean cell count (" ~ lambda*~")")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.50)) +
ggplot2::geom_hline(yintercept = AOQ[which.max(AOQ)],linetype = "dashed") +
  ggplot2::annotate("text", x = 4*lambda[which.max(AOQ)], y = AOQ[which.max(AOQ)], label = sprintf("\n AOQL = %0.4f", round(AOQ[which.max(AOQ)], digits = 4)), size = 3)
}

