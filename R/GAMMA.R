# GAMMA.R
# All function related to fitting the Gamma distribution

gamma_mle <- function(dat) {
# Fit GAMMA distribution with Maximum Likelihood Estimator
# Returns param as a list($estimate, $se)
  
  param <- list(estimate = c(NA, NA), se = c(NA, NA))    
  if (length(dat) >= GLOBAL_min_years_data) {
    
    fail_safe <- failwith(NULL, fitdistr) 
    temp.param <- fail_safe(dat, "gamma")
    if (!is.null(temp.param)) {
      param$estimate <- temp.param$estimate
      param$se <- temp.param$sd 
      invisible(param)  
    } else {
      print("Warning: the function fitdistr failed in gamma_mle")
      invisible(param)
    }
    
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))  
    invisible(param)
  }
}

gamma_Lmom <- function(dat) {
# Fit GAMMA distribution with Lmoment estimator
# Returns param as a list($estimate, $se)
  
  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    # param <- mmedist(dat,"gamma")
    
    fail_safe <- failwith(NULL, mmedist) 
    temp.param <- fail_safe(dat, "gamma")
    if (!is.null(temp.param)) {
      # Standard error is not yet implemented 
      param$estimate <- temp.param$estimate
      invisible(param)  
    } else {
      print("Warning: The function fitdistr failed in gamma_mle")
      param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
      invisible(param)
    }
    
  } else {
    print(paste("Warning: this station has less than ",GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))  
    invisible(param)
  }
}

gamma_mom <- function(dat) {
# Fit GAMMA distribution with ordinary moment estimator
# Returns param as a list($estimate, $se)

  param <- list(estimate = c(NA, NA), se = c(NA, NA))     
  if (length(dat) >= GLOBAL_min_years_data) {  
    sigma <- moments(dat)[2]
    mu <-  moments(dat)[1]
    param$estimate <- c()
    param$estimate[2] <- mu / sigma^2
    param$estimate[1] <- mu^2 / sigma^2
    # Standard error is not yet implemented
    invisible(param)
  } else {
    print(paste("Warning: this station has less than ",GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))  
    invisible(param)
  }
}

gamma_bayes <- function(dat) {
# Fit GAMMA distribution with the Bayesian method. This is a dummy function because this method has not been implemented yet  
  
  param <- list(estimate = c(NA, NA), se = c(NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {
    print("Returning NAs. This function has not been implemented yet!")
    invisible(param)
  } else {
  print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
              collapse = "", sep = ""))  
  invisible(param)  
  }
}