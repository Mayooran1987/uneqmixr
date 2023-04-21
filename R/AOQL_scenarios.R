##' \code{\link{AOQL_scenarios}} provides the Average Outgoing Quality (AOQ) curve and calculates Average Outgoing Quality Level (AOQL) value based on expected microbial counts in each scenario.
##' @title Construction of  AOQ curve and calculate AOQL value based on average microbial counts
##' @param c acceptance number
##' @param llim the upper limit for graphing the arithmetic mean of cell count
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m1 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m2 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param n number of aggregate samples which are used for inspection.
##' @param scenario what scenario we have considered such as  \code{"1"}or \code{"2"}or \code{"3"}
##' @param n number of aggregate samples which are used for inspection.
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details  Since \eqn{p_a} is the probability of acceptance, \eqn{\lambda} is the arithmetic mean of cell count and the outgoing contaminated arithmetic mean of cell count of incremental samples is given by \eqn{AOQ} as the product \eqn{\lambda p_a}.
##'           The quantity \eqn{AOQL} is defined as the maximum proportion of outgoing contaminated incremental samples and is given by \deqn{AOQL ={\max_{\lambda \geq 0}}{\lambda p_a}}
##' @return AOQ curve and AOQL value based on expected microbial counts in each scenario.
##' @seealso  \link{scenario_1_pa}, \link{scenario_2_pa}, \link{scenario_3_pa}
##' @examples
##' c <- 0
##' llim <- 0.02
##' sd <- 0.8
##' m1 <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10)
##' m2 <- c(15,5,5,5,10,5,10,5,15,10,5,10,5,25,10,5,10,5,5,10,5,15,10,
##' 5,5,20,5,10,5,10,20,5,10,30,5,20,5,10,5,10,20,15,10,15,10,
##' 10,5,10,15,5)
##' n <- 10
##' AOQL_scenarios(c,llim, sd, m1, m2, scenario = "1", n, type = "theory")
##' @usage  AOQL_scenarios(c,llim, sd, m1, m2, scenario, n, type, K, n_sim)
##' @export
AOQL_scenarios <- function(c, llim, sd, m1, m2, scenario, n, type, K = NA, n_sim = NA) {
  Sampling_scheme <- NULL
  # f_spr <- function(k,M) {
  #   sprintf("Scheme(k=%.0f, M=%.0f)", k, M)
  # }
  f_spr <- function(m) {
    M <- sum(m)
    k <- length(m)
    if (var(m) == 0) {
      sprintf("Scheme (equal, k=%.0f, M=%.0f)", k, M)
    } else {
      sprintf("Scheme (un-equal, k=%.0f, M=%.0f)", k, M)
    }
  }
  lambda <- seq(0, llim, by = 1e-04)
  mu <- log(lambda, 10) - (sd^2 / 2) * log(10, exp(1))
  if (type == "theory") {
    if (scenario == "1") {
      AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
      for (i in 1:length(mu)) {
        AOQ[i, 1] <- lambda[i] * scenario_1_pa(c, mu[i], sd, m1, n, type = "theory")
        AOQ[i, 2] <- lambda[i] * scenario_1_pa(c, mu[i], sd, m2, n, type = "theory")
      }
    } else if (scenario == "2") {
      AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
      for (i in 1:length(mu)) {
        AOQ[i, 1] <- lambda[i] * scenario_2_pa(c, mu[i], sd, m1, n, type = "theory")
        AOQ[i, 2] <- lambda[i] * scenario_2_pa(c, mu[i], sd, m2, n, type = "theory")
      }
    } else if (scenario == "3") {
      if (is.na(K) == TRUE) {
        stop("please specify the value of the shape parameter (K) of the Poisson gamma distribution")
      } else {
        AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
        for (i in 1:length(mu)) {
          AOQ[i, 1] <- lambda[i] * scenario_3_pa(c, mu[i], sd, m1, K, n, type = "theory")
          AOQ[i, 2] <- lambda[i] * scenario_3_pa(c, mu[i], sd, m2, K, n, type = "theory")
        }
      }
    } else {
      print("please include what scenario you would like to consider")
    }
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      if (scenario == "1") {
        AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
        for (i in 1:length(mu)) {
          AOQ[i, 1] <- lambda[i] * scenario_1_pa(c, mu[i], sd, m1, n, type = "simulation", n_sim)
          AOQ[i, 2] <- lambda[i] * scenario_1_pa(c, mu[i], sd, m2, n, type = "simulation", n_sim)
        }
      } else if (scenario == "2") {
        AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
        for (i in 1:length(mu)) {
          AOQ[i, 1] <- lambda[i] * scenario_2_pa(c, mu[i], sd, m1, n, type = "simulation", n_sim)
          AOQ[i, 2] <- lambda[i] * scenario_2_pa(c, mu[i], sd, m2, n, type = "simulation", n_sim)
        }
      } else if (scenario == "3") {
        if (is.na(K) == TRUE) {
          stop("please specify the value of the shape parameter(K) of the Poisson gamma distribution")
        } else {
          AOQ <- matrix(NA, nrow = length(mu), ncol = 2)
          for (i in 1:length(mu)) {
            AOQ[i, 1] <- lambda[i] * scenario_3_pa(c, mu[i], sd, m1, K, n, type = "simulation", n_sim)
            AOQ[i, 2] <- lambda[i] * scenario_3_pa(c, mu[i], sd, m2, K, n, type = "simulation", n_sim)
          }
        }
      } else {
        print("please include what scenario you would like to consider")
      }
    }
  }
  Prob <- data.frame(lambda, AOQ)
  # colnames(Prob ) <- c("lambda", f_spr(k,M))
  colnames(Prob) <- c("lambda", f_spr(m1), f_spr(m2))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "AOQ")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = lambda, y = AOQ, group = Sampling_scheme, colour = Sampling_scheme)) +
    ggplot2::theme_classic() +
    ggplot2::ylab(expression(AOQ)) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::xlab(expression("arithmetic mean cell count (" ~ lambda[0] * ~")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.75)) +
    ggplot2::geom_hline(yintercept = c(AOQ[, 1][which.max(AOQ[, 1])], AOQ[, 2][which.max(AOQ[, 2])]), linetype = "dashed") +
    ggplot2::annotate("text",
      x = c(4 * lambda[which.max(AOQ[, 1])], 4 * lambda[which.max(AOQ[, 2])]),
      y = c(AOQ[, 1][which.max(AOQ[, 1])], AOQ[, 2][which.max(AOQ[, 2])]), label = sprintf("\n AOQL = %0.4f", round(c(AOQ[, 1][which.max(AOQ[, 1])], AOQ[, 2][which.max(AOQ[, 2])]), digits = 4)), size = 3
    )

  # plot_sam
  return(plot_sam)
}
