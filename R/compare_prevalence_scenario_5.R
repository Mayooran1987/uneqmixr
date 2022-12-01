##' \code{\link{compare_prevalence_scenario_5}} provides graphical displays based on prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study.
##' @title Graphical displays of prevalence before inspection of different sampling schemes under scenario 5.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param m1 the vector of the first set of incremental samples (with equal/unequal weights).
##' @param m2 the vector of the second set of incremental samples (with equal/unequal weights).
##' @param l the number of lots in the production process.
##' @param type what type of the results you would like to consider such as "theory" or "simulation".
##' @param n_sim number of simulations (large simulations provide more precise estimation).
##' @details \code{\link{compare_prevalence_scenario_5}} provides graphical displays based on prevalence before inspection under scenario 5 of modelling the quantity of material sampled in the risk assessment study. (this section will be updated later on)
##' @return Graphical displays based on prevalence before inspection when lot with heterogeneous contamination and contamination levels fluctuate from lot to lot by using theoretical or simulation-based results.
##' @examples
##' mulow <- -6
##' muhigh <- 1
##' m1 <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
##' 10,10,10,10,10,10,10,10,10,10)
##' m2 <- c(15,5,5,5,10,5,10,5,15,10,5,10,5,25,10,5,10,5,5,10,5,15,10,
##' 5,5,20,5,10,5,10,20,5,10,30,5,20,5,10,5,10,20,15,10,15,10,
##' 10,5,10,15,5)
##' l <- 2000
##' compare_prevalence_scenario_5(mulow, muhigh, sd = 0.8, m1, m2, l,
##' type = "theory")
##' @usage  compare_prevalence_scenario_5(mulow, muhigh, sd, m1, m2, l, type, n_sim)
##' @export
compare_prevalence_scenario_5 <- function(mulow, muhigh, sd = 0.8, m1, m2, l, type, n_sim = NA){
  Sampling_scheme <- NULL
  Prev <- NULL
  f_spr <- function(m) {
    M <- sum(m)
    k <- length(m)
    if (var(m) == 0) {
      sprintf("Scheme (equal, k=%.0f, M=%.0f)",k, M)
    } else{
      sprintf("Scheme (un-equal, k=%.0f, M=%.0f)",k, M)
    }
  }
  mu <- seq(mulow, muhigh, 0.001)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  if (type == "theory") {
    # message("Not yet established, please use simulation-based results")
    Prev <- matrix(NA, nrow = length(mu), ncol = 2)
    for (i in 1:length(mu)) {
      # for (j in 1:2) {
      Prev[i,1] <-  scenario_5_prevalence(mu[i], sd = 0.8, m1,l, type)
      Prev[i,2] <-  scenario_5_prevalence(mu[i], sd = 0.8, m2,l, type)
      # }
    }
    Prob <- data.frame(mu, Prev)
    colnames(Prob ) <- c("mu", f_spr(m1),f_spr(m2))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "Prev")
    plot_sam <- ggplot2::ggplot(melten.Prob,ggplot2::aes(x = mu, y = Prev, group = Sampling_scheme, colour = Sampling_scheme)) +
      ggplot2::geom_line() +
      # ggplot2::geom_smooth() +
      # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
      ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(Prevalance)) + ggthemes::scale_colour_colorblind() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                     axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
      ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
                                                               labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
    return(plot_sam)
  } else if (type == "simulation") {
    if (is.na(n_sim) == TRUE) {
      stop("please set the number of simualtions")
    } else {
      Prev <- matrix(NA, nrow = length(mu), ncol = 2)
      for (i in 1:length(mu)) {
        # for (j in 1:2) {
        Prev[i,1] <-  scenario_5_prevalence(mu[i], sd = 0.8, m1,l, type, n_sim)
        Prev[i,2] <-  scenario_5_prevalence(mu[i], sd = 0.8, m2,l, type, n_sim)
        # }
      }
      Prob <- data.frame(mu, Prev)
      colnames(Prob ) <- c("mu", f_spr(m1),f_spr(m2))
      melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "Prev")
      plot_sam <- ggplot2::ggplot(melten.Prob,ggplot2::aes(x = mu, y = Prev, group = Sampling_scheme, colour = Sampling_scheme)) +
        ggplot2::geom_line() +
        # ggplot2::geom_smooth() +
        # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
        ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(Prevalance)) + ggthemes::scale_colour_colorblind() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.25, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                       axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
        ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
                                                                 labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
      return(plot_sam)
    }
    } else {
    print("please include what type (theory/ simulation) you would like to consider")
  }
}

