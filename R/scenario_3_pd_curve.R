##' \code{\link{scenario_3_pd_curve}} provides the probability of detection curves under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' @title Construction of  Operating Characteristic (OC) curve under lot with heterogeneous and low-level contamination.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m1 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m2 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param K shape parameter (default value 0.25).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_3_pd_curve}} provides the probability of detection curves under scenario 3 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return probability of detection curves when lot with heterogeneous and low-level contamination.
##' @seealso  \link{scenario_3_pd}
##' @examples
##' mulow <- -5
##' muhigh <- 2
##' sd <- 0.8
##' m1 <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
##' 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
##' m2 <- c(10,15,15,20,20,15,15,15,15,20,10,10,20,20,10,20,10,15,15,20,20,15,
##' 15,15,15,20,10,10,20,20,10,20)
##' K <- 0.05
##' scenario_3_pd_curve(mulow, muhigh, sd, m1, m2, K, type = "theory")
##' @usage  scenario_3_pd_curve(mulow, muhigh, sd, m1, m2, K, type, n_sim)
##' @export
scenario_3_pd_curve <- function(mulow, muhigh, sd = 0.8, m1, m2, K = 0.25, type, n_sim = NA) {
  p_d <- NULL
  Sampling_scheme <- NULL
  # M <- sum(m1)
  # if (sum(m1) != sum(m2))
  #   stop(" weights of aggregate samples are not equal, Please check it")
  # k <- c(length(m1),length(m2))
  # f_spr <- function(k,M) {
  #   sprintf("Scheme(k=%.0f, M=%.0f)", k, M)
  # }
  f_spr <- function(number, m) {
    M <- sum(m)
    k <- length(m)
    if (var(m) == 0) {
      sprintf("Scheme %.0f(equal, k=%.0f, M=%.0f)", number, k, M)
    } else {
      sprintf("Scheme %.0f(un-equal, k=%.0f, M=%.0f)", number, k, M)
    }
  }
  mu <- seq(mulow, muhigh, 0.01)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))

  if (type == "theory") {
    Pd <- matrix(NA, nrow = length(mu), ncol = 2)
    for (i in 1:length(mu)) {
      # for (j in 1:2) {
      Pd[i, 1] <- scenario_3_pd(mu[i], sd = 0.8, m1, K, type = "theory")
      Pd[i, 2] <- scenario_3_pd(mu[i], sd = 0.8, m2, K, type = "theory")
      # }
    }
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      Pd <- matrix(NA, nrow = length(mu), ncol = 2)
      for (i in 1:length(mu)) {
        # for (j in 1:2) {
        Pd[i, 1] <- scenario_3_pd(mu[i], sd = 0.8, m1, K, type = "simulation", n_sim)
        Pd[i, 2] <- scenario_3_pd(mu[i], sd = 0.8, m2, K, type = "simulation", n_sim)
        # }
      }
    }
  } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
  Prob <- data.frame(mu, Pd)
  # colnames(Prob ) <- c("mu", f_spr(k,M))
  colnames(Prob) <- c("mu", f_spr(1, m1), f_spr(2, m2))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("log mean concentration  (" ~ mu[0] * ~")")) +
    ggplot2::ylab(expression(p[d])) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    ) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
      name = "expected cell counts (cfu/g)", breaks = seq(min(mu), max(mu), 1),
      labels = c(sprintf("%f", 10^(seq(min(mu), max(mu), 1) + (sd^2 / 2) * log(10, exp(1)))))
    ))
  # plot_sam
  return(plot_sam)
}
