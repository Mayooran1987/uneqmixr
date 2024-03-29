##' \code{\link{compare_ex_var_scenario_3}} provides graphical displays based on expectation or variance under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' @title Graphical displays of expectation or variance of different sampling schemes under scenario 3.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m1 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m2 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param K shape parameter (default value 0.25).
##' @param measure what type of measure you would like to consider for the graph, such as "expectation" or "variance".
##' @details \code{\link{compare_ex_var_scenario_3}} provides a probability of acceptance under scenario 3 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with homogeneous contaminations), we employed Poisson gamma distribution to the model number of micro-organisms in the incremental samples. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Graphical displays based on expectation or variance when lot with heterogeneous and low-level contamination.
##' @examples
##' mulow <- 0
##' muhigh <- 2
##' m1 <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10)
##' m2 <- c(15,5,5,5,10,5,10,5,15,10,5,10,5,25,10,5,10,5,5,10,5,15,10,
##' 5,5,20,5,10,5,10,20,5,10,30,5,20,5,10,5,10,20,15,10,15,10,
##' 10,5,10,15,5)
##' K <- 0.05
##' compare_ex_var_scenario_3(mulow, muhigh, sd = 0.8, m1, m2, K, measure = "variance")
##' compare_ex_var_scenario_3(mulow, muhigh, sd = 0.8, m1, m2, K, measure =  "expectation")
##' @usage  compare_ex_var_scenario_3(mulow, muhigh, sd, m1, m2, K, measure)
##' @export
compare_ex_var_scenario_3 <- function(mulow, muhigh, sd = 0.8, m1, m2, K, measure = "variance") {
  variance <- NULL
  expectation <- NULL
  Sampling_scheme <- NULL
  # if (sum(m1) != sum(m2))
  #   stop(" weights of aggregate samples are not equal, Please check it")
  # k <- c(length(m1),length(m2))
  f_spr <- function(m) {
    M <- sum(m)
    k <- length(m)
    if (var(m) == 0) {
      sprintf("Scheme (equal, k=%.0f, M=%.0f)", k, M)
    } else {
      sprintf("Scheme (un-equal, k=%.0f, M=%.0f)", k, M)
    }
  }
  # f_spr(1,m2)
  mu <- seq(mulow, muhigh, 0.01)
  lambda <- 10^(mu + (sd^2 / 2) * log(10, exp(1)))

  scenario_3_expectation <- function(lambda, K, m) {
    expect <- (K / sum(m)) * sum(m * (m / min(m)) * lambda)
    return(expect)
  }
  scenario_3_variance <- function(lambda, K, m) {
    var <- (K / (sum(m) * sum(m))) * ((sum(m * m * m) * lambda / min(m)) + (sum(m * m * m * m) * lambda * lambda / (min(m) * min(m))))
    return(var)
  }
  # mu <- 2
  # lambda <- 0.54
  # K <- 0.05
  # m <- c(10,20,20,10,10,10,20,10,10,10,20,20,10,10,10,20,10,10,10)
  # scenario_3_expectation(lambda,K,m)
  # scenario_3_variance(ambda,K,m)

  Ex <- matrix(NA, nrow = length(lambda), ncol = 2)
  for (i in 1:length(lambda)) {
    # for (j in 1:2) {
    Ex[i, 1] <- scenario_3_expectation(lambda[i], K, m1)
    Ex[i, 2] <- scenario_3_expectation(lambda[i], K, m2)
    # }
  }
  Va <- matrix(NA, nrow = length(lambda), ncol = 2)
  for (i in 1:length(lambda)) {
    # for (j in 1:2) {
    Va[i, 1] <- scenario_3_variance(lambda[i], K, m1)
    Va[i, 2] <- scenario_3_variance(lambda[i], K, m2)
    # }
  }
  if (measure == "variance") {
    Prob <- data.frame(mu, Va)
    colnames(Prob) <- c("mu", f_spr(m1), f_spr(m2))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "variance")
    plot_sam <- ggplot2::ggplot(melten.Prob) +
      ggplot2::geom_line(ggplot2::aes(x = mu, y = variance, group = Sampling_scheme, colour = Sampling_scheme)) +
      # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression("log mean concentration  (" ~ mu[0] * ~")")) +
      ggplot2::ylab(expression(Variance)) +
      ggthemes::scale_colour_colorblind() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
        axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
      ) +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
        name = "expected cell counts (cfu/g)", breaks = seq(min(mu), max(mu), 1),
        labels = c(sprintf("%f", 10^(seq(min(mu), max(mu), 1) + (sd^2 / 2) * log(10, exp(1)))))
      ))
  } else if (measure == "expectation") {
    Prob <- data.frame(mu, Ex)
    colnames(Prob) <- c("mu", f_spr(m1), f_spr(m2))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "expectation")
    plot_sam <- ggplot2::ggplot(melten.Prob) +
      ggplot2::geom_line(ggplot2::aes(x = mu, y = expectation, group = Sampling_scheme, colour = Sampling_scheme)) +
      # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression("log mean concentration  (" ~ mu[0] * ~")")) +
      ggplot2::ylab(expression(Expectation)) +
      ggthemes::scale_colour_colorblind() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
        axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")
      ) +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~.,
        name = "expected cell counts (cfu/g)", breaks = seq(min(mu), max(mu), 1),
        labels = c(sprintf("%f", 10^(seq(min(mu), max(mu), 1) + (sd^2 / 2) * log(10, exp(1)))))
      ))
  } else {
    print("Please choose measure as expectation or variance")
  }
  return(plot_sam)
}
