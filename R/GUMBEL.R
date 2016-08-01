# GUMBEL.R
# All function related to fitting the Gumbel distribution


#' Fitting the Gumbel distribution with MLE
#' @description Function to fit the Gumbel distribution with the maximum likelihood method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters (2) and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gumbel_mle(XXX)
gumbel_mle <- function(dat) {
  
  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    param <- gum.fit(dat, show = FALSE)
    param$estimate <- param$mle
    invisible(param)

  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse = "", sep = ""))   
    invisible(param)
    }
}


#' Fitting the Gumbel distribution with Lmom
#' @description Function to fit the Gumbel distribution with the linear moment method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters (2) and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gumbel_Lmom(XXX))
gumbel_Lmom <- function(dat) {

  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    dat.mom <- Lmoments(dat)
    param$estimate <- invisible(as.numeric(par.gumb(dat.mom[1],dat.mom[2]))) # Standard error not yet implemented
    invisible(param)
  } else {
    print(paste("Warning: this station has less than ",GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))    
    invisible(param)
  }
}

#' Fitting the Gumbel distribution with mom
#' @description Function to fit the Gumbel distribution with the ordinary moments method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters (2) and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gumbel_mom(XXX))
gumbel_mom <- function(dat) {

  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    sigma <- moments(dat)[2]
    mu <-  moments(dat)[1]
    param$estimate[2] <- sigma / sqrt(1.645)
    param$estimate[1] <- mu - 0.5772 * param$estimate[2]
    # Standard error is not yet implemented
    invisible(param)
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))  
    invisible(param)
  }
}

#' Fitting the Gumbel distribution with Bayesian inference
#' @description Function to fit the Gumbel distribution with BayesianMCMC method
#' We do not provide a prior for the BayesianMCMC function
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gumbel_bayes(XXX)
gumbel_bayes <- function(dat) {

  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  
  if (length(dat) >= GLOBAL_min_years_data) {

    # Should we do gumbel with a prior? I tried, it didn't work
    fail_safe <- failwith(NULL, BayesianMCMC)
    bayes <- fail_safe(dat, nbpas = 5000, nbchaines = 2, 
                       confint = c(0.05, 0.95), dist = "GUMBEL")
    
    if (is.null(bayes) == TRUE) {
      print("Warning: the function BayesianMCMC failed in gumbel_bayes")
      invisible(param)
    } else {
      
      ## Addition to return parameters
      #   # Solution 1
      param$estimate <- bayes$parametersML

#       # Solution 2
#       param$estimate[1] <- mean(as.vector(bayes$parameters[, 1, 1:3]))
#       param$estimate[2] <- mean(as.vector(bayes$parameters[, 2, 1:3]))

      param$se[1] <- sd(as.vector(bayes$parameters[, 1, ]))
      param$se[2] <- sd(as.vector(bayes$parameters[, 2, ]))

      invisible(param)
    }
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                collapse = "", sep = ""))  
    invisible(param)  
  }
}

#' Calculating the posterior predictive distribution for Gumbel
#' @description Function to calculate the posterior predictive distribution after calling gumbel_bayes
#' @param (mmrp, mupars, spars) parameters returned by gumbel_bayes. mupars, spars are the ensemble of param$estimate
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples Needs example
get_posterior_gumbel <- function(mmrp, mupars, spars) {

  qqsample1 <- sapply(seq(length(mupars)), function(st) {
    mean_temp = mupars[st]
    st_temp = spars[st]
    invF.gumb(F = (1 - 1 / mmrp), mean_temp, st_temp)
    }, simplify = "array")
  # 1 is for collums only 0.5 to return the median
  qqr <- apply(qqsample1, 1, quantile, 0.5)
}

