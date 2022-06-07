##' \code{\link{scenario_1A_OC}} provides the Operating Characteristic (OC) curves under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' @title Construction of  Operating Characteristic (OC) curve under lot with homogeneous contaminations.
##' @param c acceptance number
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation on the log10 scale (default value 0.8).
##' @param M the quantity (weight) of the aggregate sample.
##' @param m1 the quantity (weight) of the based sample for the risk assessment (for this research, we used a quantity of sample which is to be the minimum quantity of selected incremental samples).
##' @param n number of aggregate samples which are used for inspection.
##' @details \code{\link{scenario_1A_OC}} provides the Operating Characteristic (OC) curves under scenario 1 of modelling the quantity of material sampled in the risk assessment study.
##' Under this scenario (lot with homogeneous contaminations), we employed Poisson distribution to the model number of micro-organisms in the incremental samples.
##' The purpose of this function used for compares different sets of sampling schemes when lot with heterogeneous and high-level contamination.
##' Under this scenario expected cell count in each incremental sample can be written in terms of the based incremental sample's expected cell count (for this study, the based sample is a sample that is to be the minimum quantity).
##' The probability of acceptance is plotted against mean log10 concentration and expected cell counts. Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated later on)
##' @return Operating Characteristic (OC) curves under lot with homogeneous contaminations.
##' @seealso  \link{scenario_1A_pd}
##' @examples
##' c <- 0
##' mulow <- -6
##' muhigh <- 0
##' sd <- 0.8
##' M <- 60
##' m1 <- c(5,10,15)
##' n <- 10
##' scenario_1A_OC(c, mulow, muhigh, sd, M, m1, n)
##' @usage  scenario_1A_OC(c, mulow, muhigh, sd, M, m1, n)
##' @export
scenario_1A_OC <- function(c, mulow, muhigh, sd = 0.8, M, m1, n){
  P_a <- NULL
  Sampling_scheme <- NULL
  f_spr <- function(m1,m) {
    sprintf("Scheme(m_1=%.0f, M=%.0f)", m1, M)
  }
  mu <- seq(mulow, muhigh, 0.01)
  lambda <- 10^(mu + ((sd^2)/2) * log(10, exp(1)))
  # pd <- 1-exp(lambda*(-m/m1))
  # pd <- pd1(lambda,m,m1)


  Pa <- matrix(NA, nrow = length(mu), ncol = length(m1))
  for (j in 1:length(m1)) {
    for (i in 1:length(mu)) {
      # Pa[i,j] <-  (1 - scenario_1A_pd(lambda[i],M,m1[j]))^n
      Pa[i,j] <-  scenario_1A_pa(c, lambda[i], M, m1[j], n)
    }
  }

  # pa <- (1-pd)^n
  Prob <- data.frame(mu, Pa)
  colnames(Prob ) <- c("mu", f_spr(m1,M))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "expected cell counts (cfu/g)", breaks = seq(min(mu),max(mu),1),
                                                             labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  return(plot_sam)
}
