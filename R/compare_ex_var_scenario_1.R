##' \code{\link{compare_ex_var_scenario_1}} provides graphical displays based on expectation or variance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Graphical displays of expectation or variance of different sampling schemes under scenario 1.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m1 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m2 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param measure what type of measure you would like to consider for the graph, such as "expectation" or "variance".
##' @details \code{\link{compare_ex_var_scenario_1}} provides graphical displays based on expectation or variance under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (a lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples. (this section will be updated later on)
##' @return Graphical displays based on expectation or variance under lot with homogeneous contaminations.
##' @examples
##' mulow <- -1
##' muhigh <- 1
##' sd <- 0.8
##' m1 <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10)
##' m2 <- c(15,5,5,5,10,5,10,5,15,10,5,10,5,25,10,5,10,5,5,10,5,15,10,
##' 5,5,20,5,10,5,10,20,5,10,30,5,20,5,10,5,10,20,15,10,15,10,
##' 10,5,10,15,5)
##' compare_ex_var_scenario_1(mulow, muhigh, sd = 0.8, m1, m2, measure = "variance")
##' compare_ex_var_scenario_1(mulow, muhigh, sd = 0.8, m1, m2, measure =  "expectation")
##' @usage  compare_ex_var_scenario_1(mulow, muhigh, sd, m1, m2, measure)
##' @export
compare_ex_var_scenario_1 <- function(mulow, muhigh, sd = 0.8, m1, m2, measure = "variance"){
  variance <- NULL
  expectation <- NULL
  Sampling_scheme <- NULL
  # M1 <- sum(m1)
  # M2 <- sum(m2)
  # if (sum(m1) != sum(m2))
  #   stop(" weights of aggregate samples are not equal, Please check it")
  # k <- c(length(m1),length(m2))
  f_spr <- function(m) {
    M <- sum(m)
    k <- length(m)
    if (var(m) == 0) {
      sprintf("Scheme (equal, k=%.0f, M=%.0f)", k, M)
    } else{
      sprintf("Scheme (un-equal, k=%.0f, M=%.0f)", k, M)
    }
  }
  # f_spr(1,m2)
  mu <- seq(mulow, muhigh, 0.01)
  lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))

  Ex <- matrix(NA, nrow = length(mu), ncol = 2)
  for (i in 1:length(mu)) {
    # for (j in 1:2) {
    Ex[i,1] <-  lambda[i]*sum(m1*m1/(min(m1)*sum(m1)))
    Ex[i,2] <-  lambda[i]*sum(m2*m2/(min(m2)*sum(m2)))
    # }
  }

  Va <- matrix(NA, nrow = length(mu), ncol = 2)
  for (i in 1:length(mu)) {
    # for (j in 1:2) {
    Va[i,1] <-  lambda[i]*sum(m1*m1*m1/(min(m1)*sum(m1)*sum(m1)))
    Va[i,2] <-  lambda[i]*sum(m2*m2*m2/(min(m2)*sum(m2)*sum(m2)))
    # }
  }
  if (measure == "variance") {
    Prob <- data.frame(mu, Va)
    colnames(Prob ) <- c("mu", f_spr(m1), f_spr(m2))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "variance")
    plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = variance, group = Sampling_scheme, colour = Sampling_scheme)) +
      # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
      ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu[0]*~")")) + ggplot2::ylab(expression(Variance)) + ggthemes::scale_colour_colorblind() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                     axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),

                                                               labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  } else if (measure == "expectation") {
    Prob <- data.frame(mu, Ex)
    colnames(Prob ) <- c("mu", f_spr(m1), f_spr(m2))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "expectation")
    plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = expectation, group = Sampling_scheme, colour = Sampling_scheme)) +
      # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
      ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu[0]*~")")) + ggplot2::ylab(expression(Expectation)) + ggthemes::scale_colour_colorblind() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                     axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),

                                                               labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  } else{
    print("Please choose measure as expectation or variance")
  }
  return(plot_sam)

}

