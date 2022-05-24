##' \code{\link{scenarious_1_OC}} provides the Operating Characteristic (OC) curve for known microbiological distribution such as Poisson.
##' The probability of acceptance is plotted against mean log10 concentration and expected cell counts.
##' @title Construction of  Operating Characteristic (OC) curve under lot with homogeneous contaminations.
##' @param mulow the lower value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param muhigh the upper value of the mean concentration (\eqn{\mu}) for use in the graphical display's x-axis.
##' @param sd standard deviation of the lognormal and Poisson-lognormal distributions on the log10 scale (default value 0.8)
##' @param m the quantity (weight) of the aggregate sample
##' @param m1 the quantity (weight) of the based sample for the risk assessment (for this research, we used a quantity of sample which is to be the minimum quantity of selected incremental samples)
##' @param n number of aggregate samples which are used for inspection
##' @details Based on the food safety literature, the expected cell count is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}. (this section will be updated soon)
##' @return Operating Characteristic (OC) curves under lot with homogeneous contaminations
##' @examples
##' mulow <- -6
##' muhigh <- 0
##' sd <- 0.8
##' m <- 25
##' m1 <- c(5,10,15)
##' n <- 20
##' scenarious_1_OC(mulow,muhigh,sd,m,m1,n)
##' @usage  scenarious_1_OC(mulow,muhigh,sd,m,m1,n)
##' @export
scenarious_1_OC <- function(mulow,muhigh,sd,m,m1,n){
  P_a <- NULL
  Sampling_scheme <- NULL
  f_spr <- function(m1,m) {
    sprintf("Scheme(m_1=%.0f, m=%.0f)", m1, m)
  }
  mu <- seq(mulow, muhigh, 0.01)
  lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  # pd <- 1-exp(lambda*(-m/m1))
  # pd <- pd1(lambda,m,m1)

  pd1 <- function(lambda,m,m1){
    pd <- 1 - exp(lambda*(-m/m1))
    return(pd)
  }
  Pa <- matrix(NA, nrow = length(mu), ncol = length(m1))
  for (j in 1:length(m1)) {
    Pa[,j] <-  (1 - pd1(lambda,m,m1[j]))^n
  }

  # pa <- (1-pd)^n
  Prob <- data.frame(mu, Pa)
  colnames(Prob ) <- c("mu", f_spr(m1,m))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
  plot_sam <- ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("OC curve based on Lognormal distribution") +
    ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.75), axis.line.x.top = ggplot2::element_line(color = "red"),
                   axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean concentration (cfu/g)", breaks = seq(min(mu),max(mu),1),
                                                             labels = c(sprintf("%f", 10^(seq(min(mu),max(mu),1) + (sd^2/2) * log(10, exp(1)))))))
  return(plot_sam)
}

