# GEV.R
# All functions related to fitting the Generalized Extreme Value distribution
# All distribution functions names are "distrib_method"

#' Fitting the GEV distribution with MLE
#' @description Function to fit the GEV distribution with the maximum likelihood method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gev_mle(evd::rgev(10000, loc=0, scale=1, shape=0))
gev_mle <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA))    
  if (length(dat) >= GLOBAL_min_years_data) {

    
    fail_safe <- failwith(NULL, fgev)
    fitted.param <- fail_safe(dat)
   
    if (is.null(fitted.param) == TRUE) {
      print("Warning: the function fgev failed in gev_mle")
      invisible(param)
    } else {
    # fitted.param <- fgev(dat)
    param$estimate <- fitted.param$estimate
    param$se <- fitted.param$std.err
    invisible(param)
    }
    } else {
      print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse = "", sep = ""))  
      invisible(param)
    }
}

#' Fitting the GEV distribution with Lmom
#' @description Function to fit the GEV distribution with the linear moment method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se).
#' Standard error is not yet implemented
#' @export
#'
#' @examples gev_Lmom(evd::rgev(10000, loc=0, scale=1, shape=0))
gev_Lmom <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA))   
  if (length(dat) >= GLOBAL_min_years_data) {
 
    dat.Lmom <- Lmoments(dat)
    
    
    fail_safe <- failwith(NULL, par.GEV)
    fitted.param <- fail_safe(dat.Lmom[1], dat.Lmom[2], dat.Lmom[4])
    
    if (is.null(fitted.param) == TRUE) {
      print("Warning: the function par.GEV failed in gev_Lmom")
      invisible(param)
    } else {
    # fitted.param <- as.numeric(par.GEV(dat.mom[1], dat.mom[2], dat.mom[4]))
    # Creating the returning list
    param$estimate <- c(fitted.param$xi, fitted.param$alfa, - fitted.param$k)
    # Standard error is not yet implemented
    invisible(param)
    }
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse = "", sep = ""))  
    invisible(param)  
  }
}

#' Fitting the GEV distribution with mom
#' @description Function to fit the GEV distribution with the ordinary moments method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' Standard error is not yet implemented
#' @export
#'
#' @examples gev_mom(evd::rgev(10000, loc=0, scale=1, shape=0))
gev_mom <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA))  
  if (length(dat) >= GLOBAL_min_years_data) {

  gevmom <- moments(dat)
  # print(gevmom)
  kfunc <- function(k) {
    # print(k)
    if (abs(k) > 0.00000001) {
      sign(k) * ((-1) * gamma(1+3*k) + 3*gamma(1 + k) * gamma(1 + 2*k) - 2*(gamma(1 + k))^3)/
        ((gamma(1 + 2*k) - (gamma(1 + k))^2)^(3/2)) - as.numeric(gevmom[4])
    } else {
      (12*sqrt(6) * zeta(3)) / (pi^3) - as.numeric(gevmom[4])
    }
  }
  # xi <- newtonRaphson(kfunc, -0.1, tol = 0.0001)$root  # FKB: initial code
  fail_safe <- failwith(NULL, newtonRaphson)  # START FKB HACK
  xi <- fail_safe(kfunc, -0.1, tol = 0.0001)
  xi <- xi$root
    
  if (is.null(xi) == TRUE) {
    print("Warning: the function newtonRaphson failed in gev_mom")   # PB WITH locp!!
    invisible(param)
  } else {
    
      if (xi == 0) {   # END FKB HACK
      # xi <- xi$root  
      scp = (as.numeric(gevmom[2]) * sqrt(6)) / (pi)
      locp = as.numeric(gevmom[1]) - as.numeric(gevmom[2]) * 0.57721566490
      } else {
      # xi <- xi$root  
      scp = (as.numeric(gevmom[2])*abs(xi))/sqrt((gamma(1+2*xi)-(gamma(1+xi))^2))
      locp = as.numeric(gevmom[1])-(scp/xi)*(1-gamma(1+xi))
      }
  # z <- list()   # commented FKB
  # z$dat <- dat  # commented FKB
  # z$mle replaced by z$estimate to integrate with the other functions
  param$estimate <- c(locp, scp, (-1)*xi)  # z replaced by param FKB 
  # z$vals <- cbind(locp, scp, xi) # commented FKB
  # class(z) <- "gev.fit"
  # Standard error is not yet implemented
  invisible(param)
  }
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse = "",sep = ""))  
    invisible(param)  
    }
}

#' Fitting the GEV distribution with Bayesian inference
#' @description Function to fit the GEV distribution with BayesianMCMC method
#' WE assume that the shape parameter only has a prior with mean zero and standard deviation 0.2 (dnorm(x[3], 0, 0.2))
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gev_bayes(evd::rgev(10000, loc=0, scale=1, shape=0))
gev_bayes <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    # Prior for Bayes
    myprior <- function (x) {
      # x = vector of parameter values: c(location, scale, shape)
      # I assume the shape parameter only has a prior with mean zero and standard deviation 0.2
      dnorm(x[3], 0, 0.2)
    }
    
    fail_safe <- failwith(NULL, BayesianMCMC)
    bayes <- fail_safe(dat, nbpas = 5000, nbchaines = 2, confint = c(0.05, 0.95), dist = "GEV", apriori = myprior)

    if (is.null(bayes) == TRUE) {
      print("Warning: the function BayesianMCMC failed in gev_bayes")
      invisible(param)
    } else {
      
      ## Addition to return parameters
        # Solution 1
      param$estimate <- bayes$parametersML
      # Sign correction for the shape parameter to be consistent with the other fitting functions
      param$estimate[3] <- - param$estimate[3]

      # Solution 2
#       param$estimate[1] <- mean(as.vector(bayes$parameters[, 1, 1:3]))
#       param$estimate[2] <- mean(as.vector(bayes$parameters[, 2, 1:3]))
#       param$estimate[3] <- mean(as.vector(bayes$parameters[, 3, 1:3]))
      
      param$se[1] <- sd(as.vector(bayes$parameters[ , 1, ]))
      param$se[2] <- sd(as.vector(bayes$parameters[ , 2, ]))
      param$se[3] <- sd(as.vector(bayes$parameters[ , 3, ]))
      
      invisible(param)
    }
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                collapse = "", sep = ""))  
    invisible(param)  
  }
  
}

#' Calculating the posterior predictive distribution for GEV
#' @description Function to calculate the posterior predictive distribution after calling gev_bayes
#' @param (mmrp, mupars, spars, kpars) parameters returned by gev_bayes. mupars, spars, kpars are the ensemble of param$estimate
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples Needs example
get_posterior_gev <- function(mmrp, mupars, spars, kpars) {

  qqsample1 <- sapply(seq(length(mupars)), function(st) {
    mean_temp <- mupars[st]
    st_temp <- spars[st]
    k_temp <- kpars[st]
    # param$estimate <- c(mean_temp, st_temp, k_temp)
    invF.GEV(F = (1 - 1 / mmrp), mean_temp, st_temp, k_temp)
  },simplify = "array")
  # 1 is for collums only 0.5 to return the median
  qqr <- apply(qqsample1, 1, quantile, 0.5)
}
  
