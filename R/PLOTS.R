# PLOTS.R
# This file gathers all plotting functions that do not involve mapping
# SOME CALLS TO GOF.R FUNCTIONS SHOULD BE UPDATED TO USE THE PLOT_ALL FUNCTION


#' plot_all
#' @description This function creates a 2 by 2 plot by calling all common plotting functions
#' It plots for any distribution and method,
#' as well as the distribution parameters and Goodness Of Fit values
#' Returns nothing, saves nothing
#' @param dat
#' @param GOF.list
#' @param param
#' @param distr
#' @param method
#'
#' @return
#' @export
#'
#' @examples
plot_all  <- function(dat, GOF.list, param, distr = "distr", method = "method") {


  windows()
  par(mfrow = c(2, 2))
  plot_density(dat, GOF.list, param, as.character(distr))
  plot_rlevel(dat, param, as.character(distr))

  # Add a table with the goodness of fit estimations
  nbs <- matrix(round(c(GOF.list$CS, GOF.list$KS, GOF.list$AD), 2), ncol = 3)
  rownames(nbs) <- c("Goodness of fit")
  colnames(nbs) <- c("CS", "KS", "AD")
  addtable2plot(0, 0, nbs, bty = "o", bg = "lightgray", display.rownames = TRUE, xpad = 0, ypad = 0)

  # Add a table with the fitting results
  if(distr == 'gumbel' | distr == 'gamma')  {
    nbs <- matrix(round(c(param$estimate[1], param$estimate[2], param$se[1], param$se[2]), 2), ncol = 2)
    rownames(nbs) <- c("Location", "Scale")
  } else {
    nbs <- matrix(round(c(param$estimate[1], param$estimate[2], param$estimate[3], param$se[1], param$se[2], param$se[3]), 2), ncol = 2)
    rownames(nbs) <- c("Location", "Scale", "Shape")
  }
  colnames(nbs) <- c("Estimate", "Std Err")
  xmax <- max(dat)
  addtable2plot(0, xmax, nbs, bty = "o", bg = "lightgray", display.rownames = TRUE, xpad = 0, ypad = 0) #150,570

  plot_ecdf(dat, param, as.character(distr))
  text(0, 1, paste("Distrib=", as.character(distr), "/ Method=", as.character(method)), cex = 1.2, adj = 0)
  plot_qq(dat, param, as.character(distr))
}



#' plot_density
#' @description plot fitted probability density function to estimated empirical pdf
#' Returns nothing, saves nothing
#' @param dat
#' @param GOF.list
#' @param param
#' @param distr
#'
#' @return
#' @export
#'
#' @examples
plot_density  <- function(dat, GOF.list, param, distr = "distr") {

  xmax <- max(dat)*1.2
  x <- seq(0, xmax, xmax / 100)

  ymax <- max(density(dat)$y)*1.2

  # Plotting input dat, this is common to all distributions
  hist(dat, xlab = "Flood discharge (m3/s)",ylab = "Probability density",freq = FALSE,
       breaks = seq(0, xmax, xmax / 15), col = "gray", main = NULL, xlim = c(0, xmax), ylim = c(0, ymax))
  par(new = TRUE)

  # Distribution specific y vector
  if(distr == 'gumbel')   y <- dgumbel(x, param$estimate[1], param$estimate[2])
  if(distr == 'gamma')    y <- dgamma(x, param$estimate[1], param$estimate[2])
  if(distr == 'gev')      y <- evd::dgev(x, param$estimate[1], param$estimate[2], param$estimate[3])  # I should have done that for most functions coming from packages...
  if(distr == 'gl')       y <- f.genlogis(x, param$estimate[1], param$estimate[2], param$estimate[3])
  if(distr == 'pearson')  y <- f.gamma(x, param$estimate[1], param$estimate[2], param$estimate[3])

  plot(x, y, xlim = c(0, xmax), ylim = c(0, ymax), type = "l", lwd = 2, col = "black", xlab = "", ylab = "")
  par(new = TRUE)
  plot(density(dat), main = "Density distribution and data histogramm",
       xlim = c(0, xmax), ylim = c(0, ymax), lty = 1, lwd = 3, col = "blue", xlab = "", ylab = "")

 legend("topright", inset = .05, c("Model","Empirical" ), col = c("black","blue"),lty = c(1, 1),lwd=c(2, 3),
        merge = TRUE, bg = "gray90")
}


#' plot_rlevel
#' @description Plot return levels
#' @param dat
#' @param param
#' @param distr
#'
#' @return Returns nothing, saves nothing
#' @export
#'
#' @examples
plot_rlevel <- function(dat, param, distr = "distr") {

  # Common to all distributions
  xmin <- min(dat)
  xmax <- max(dat)*1.5
  y <- seq(xmin, xmax, length = 100)
  empq <- sort(dat)

  # The x vector is distribution specific
  if(distr == 'gumbel') {
    x <- 1 / (1 - pgumbel(y, param$estimate[1], param$estimate[2]))
    # empT <- 1/(1-(seq(1:length(empq))-0.44)/(length(empq))+0.12) # Gringorten, optimized for the gumbel distribution
    empT <- 1/(1 - (seq(1:length(empq)) - 0.50) / (length(empq)))   # Hazen, a traditional choice
  }
  if(distr == 'gamma') {
    x <- 1 / (1 - pgamma(y, param$estimate[1], param$estimate[2]))
    empT <- 1/(1 - (seq(1:length(empq)) - 0.50) / (length(empq)))   # Hazen, a traditional choice
  }
  if(distr == 'gev')  {
    x <- 1 / (1 - evd::pgev(y, param$estimate[1], param$estimate[2], param$estimate[3]))
    # empT <- 1/(1-(seq(1:length(empq))-0.44)/(length(empq))+0.12) # Gringorten, optimized for the gumbel distribution
    empT <- 1/(1 - (seq(1:length(empq)) - 0.50) / (length(empq)))   # Hazen, a traditional choice
  }
  if(distr == 'gl')   {
    x <- 1 / (1 - F.genlogis(y, param$estimate[1], param$estimate[2], param$estimate[3]))
    # empT <- 1/(1-(seq(1:length(empq))-0.35)/(length(empq)))  # APL
    empT <- 1/(1 - (seq(1:length(empq)) - 0.50) / (length(empq)))   # Hazen, a traditional choice
  }
  if(distr=='pearson') {
    x <- 1/(1-F.gamma(y, param$estimate[1], param$estimate[2], param$estimate[3]))
    empT <- 1/(1 - (seq(1:length(empq)) - 0.50) / (length(empq)))   # Hazen, a traditional choice
  }

  # xaxt="n" is to not plot the x axis ticks, as I specify them later
  plot(log(log(x)), y, xlim = c(0, log(log(1000))), xaxt = "n", ylim = c(0, xmax),
       main = "Return levels", xlab = "Return period (years)", ylab = "Flood discharge (m3/s)",type = "l",lwd = 2)
  tix <- c(5, 10, 20, 50, 100, 200, 500)
  axis(1, at = log(log(tix)), labels = tix)

  # plot empirical dat points
  points(log(log(empT)), empq, pch = 16, col = "blue")
  grid(nx = 7, ny = 10, lwd = 2) # grid only in y-direction

}


#' plot_ecdf
#' @description Plot estimated and empirical cumulative distribution function
#' @param dat
#' @param param
#' @param distr
#'
#' @return Returns nothing, saves nothing
#' @export
#'
#' @examples
plot_ecdf  <- function(dat, param, distr = "distr") {

  xmax <- max(dat)*1.2
  x <- seq(0, xmax, xmax / 100)

  # Distribution specific y vector
  if(distr == 'gumbel') y <- pgumbel(x, param$estimate[1], param$estimate[2])
  if(distr == 'gamma')  y <- pgamma(x, param$estimate[1], param$estimate[2])
  if(distr == 'gev')    y <- evd::pgev(x, param$estimate[1], param$estimate[2], param$estimate[3])
  if(distr == 'gl')     y <- F.genlogis(x, param$estimate[1], param$estimate[2], param$estimate[3])
  if(distr == 'pearson') y <- F.gamma(x, param$estimate[1], param$estimate[2], param$estimate[3])


  plot(ecdf(dat), main = "ECDF", xlim = c(0, xmax), ylim = c(0, 1),
       xlab = "", ylab = "", lty = 21, col = "blue")
  par(new = TRUE)
  plot(x, y, xlim = c(0, xmax), ylim = c(0, 1),
       type = "l",lwd = 2, col = "black", xlab = "Flood discharge (m3/s)", ylab = "Cumulative probability")
}


#' plot_qq
#' @description QQ plot of empiricial against modelled
#' @param dat
#' @param param
#' @param distr
#'
#' @return Returns nothing, saves nothing
#' @export
#'
#' @examples
plot_qq  <- function(dat, param, distr = "distr") {

  # Compute plotting position
  # pvalues <-(seq(1:length(dat))-0.35)/length(dat) # APL
  p.values <- (seq(1:length(dat)) - 0.5) / length(dat)   # Hazen, a traditional choice
  y <- sort(dat)

  if(distr == 'gamma')  x <- sort(rgamma(p.values, param$estimate[1], param$estimate[2]))
  if(distr == 'gumbel') {
    # pvalues <- (seq(1:length(dat))-0.44)/(length(dat)+0.12) # Gringorten, optimized for the gumbel distribution
    x <- sort(rgumbel(p.values, param$estimate[1], param$estimate[2]))
    }
  if(distr == 'gev') {
    # pvalues <- (seq(1:length(dat))-0.44)/(length(dat)+0.12) # Gringorten, optimized for the gumbel distribution
    x <- sort(evd::rgev(p.values, param$estimate[1], param$estimate[2], param$estimate[3]))
  }
  if(distr == 'gl')     x <- invF.genlogis(p.values, param$estimate[1], param$estimate[2], param$estimate[3])
  if(distr == 'pearson') x <- sort(rand.gamma(p.values, param$estimate[1], param$estimate[2], param$estimate[3]))

  plot(x, y, ylab = "Empirical flood dischare (m3/s)", xlab = "Modelled flood dischare (m3/s)",
       main = "Quantile-Quantile Plot",pch = 16, col = "blue")

  abline(0, 1, lwd = 2, col = "black")
}
