##' \code{\link{scenario_5_pd_curve}} provides the probability of detection curves under scenario 5 of modelling the quantity of material sampled in the risk assessment study.
##' @title Construction of   probability of detection curves under lot with heterogeneous contamination and concentration levels fluctuating from sublot to sublot.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd_b standard deviation of concentration level between sublots on the log10 scale (default value 0.8).
##' @param sd_w standard deviation of concentration level within sublot on the log10 scale (default value 0.8).
##' @param m11 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m22 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{scenario_5_pd_curve}} provides the probability of detection curves under scenario 5 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return probability of detection curves when lot with heterogeneous and high-level contamination.
##' @seealso   \link{scenario_5_pd}
##' @examples
##' m1 <- c(10,10,10,10,10,10,10,10,10,10)
##' m2 <- c(10,10,10,10,10,10,10,10,10,10)
##' m3 <- c(10,10,10,10,10,10,10,10,10,10)
##' m4 <- c(10,10,10,10,10,10,10,10,10,10)
##' m5 <- c(10,10,10,10,10,10,10,10,10,10)
##' m11 <- list(m1,m2,m3,m4,m5)
##' m_1 <- c(15, 5, 5, 5, 10, 5, 10, 5, 15, 10)
##' m_2 <- c(5, 10, 5, 25, 10, 5, 10, 5, 5, 10)
##' m_3 <- c(5, 15, 10, 5, 5, 20, 5, 10, 5, 10)
##' m_4 <- c(20, 5, 10, 30, 5, 20, 5, 10, 5, 10)
##' m_5 <- c(20, 15, 10, 15, 10, 10, 5, 10, 15, 5)
##' m22 <- list(m_1,m_2,m_3,m_4,m_5)
##' mulow <- -8
##' muhigh <- -1
##' sd_b <- 0.2
##' sd_w <- 0.8
##' scenario_5_pd_curve(mulow, muhigh, sd_b, sd_w, m11, m22, type = "theory")
##' scenario_5_pd_curve(mulow, muhigh, sd_b, sd_w, m11, m22, type = "simulation", n_sim = 1000000)
##' @usage  scenario_5_pd_curve(mulow, muhigh, sd_b, sd_w, m11, m22, type, n_sim)
##' @export
scenario_5_pd_curve <- function(mulow, muhigh, sd_b, sd_w = 0.8, m11, m22, type, n_sim = NA) {
  p_d <- NULL
  Sampling_scheme <- NULL
  f_spr <- function(number, m, m1) {
    M <- sum(m)
    k <- length(m)
    l <- length(m1)
    if (var(m) == 0) {
      sprintf("Scheme %.0f(equal, k=%.0f, M=%.0f, l=%.0f)", number, k, M, l)
    } else {
      sprintf("Scheme %.0f(un-equal, k=%.0f, M=%.0f, l=%.0f)", number, k, M, l)
    }
  }
  mu <- seq(mulow, muhigh, 0.1)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))

  if (type == "theory") {
    Pd <- matrix(NA, nrow = length(mu), ncol = 2)
    for (i in 1:length(mu)) {
      # for (j in 1:2) {
      Pd[i, 1] <- scenario_5_pd(mu[i], sd_b, sd_w, m11, type = "theory")
      Pd[i, 2] <- scenario_5_pd(mu[i], sd_b, sd_w, m22, type = "theory")
      # }
    }
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      Pd <- matrix(NA, nrow = length(mu), ncol = 2)
      for (i in 1:length(mu)) {
        # for (j in 1:2) {
        Pd[i, 1] <- scenario_4_pd(mu[i], sd_b, sd_w, m11, type = "simulation", n_sim)
        Pd[i, 2] <- scenario_4_pd(mu[i], sd_b, sd_w, m22, type = "simulation", n_sim)
        # }
      }
    }
  } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
  Prob <- data.frame(mu, Pd)
  # colnames(Prob ) <- c("mu", f_spr(k,M))
  colnames(Prob) <- c("mu", f_spr(1, unlist(m11), m11), f_spr(2, unlist(m22), m22))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "p_d")
  plot_sam <- ggplot2::ggplot(melten.Prob) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = p_d, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("log mean concentration  (" ~ mu * ~")")) +
    ggplot2::ylab(expression(P[d])) +
    ggthemes::scale_colour_colorblind() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75, 0.25), axis.line.x.top = ggplot2::element_line(color = "red"),
      axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
    ) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
      name = "expected cell counts (cfu/g)", breaks = seq(min(mu), max(mu), 1),
      labels = c(sprintf("%f", 10^(seq(min(mu), max(mu), 1) + (sd^2 / 2) * log(10, exp(1)))))
    ))
  # plot_sam
  return(plot_sam)
}
