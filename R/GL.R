# GL.R
# All function related to fitting the Generalized Logistics distribution

#' Fitting the GEV distribution with MLE
#' @description Function to fit the GEV distribution with the maximum likelihood method.
#' This function was copied from Kolbjorn's initial file
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' Standard error is not yet implemented
#' @export
#'
#' @examples gl_mle(XXX))
gl_mle <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 
                    mulink = identity, siglink = identity, shlink = identity, 
                    muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
                    method = "Nelder-Mead", maxit = 10000, ...) {
  
  show="FALSE"  # HACK FLO. Let's make sure we dont print anything 
  z <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA))  # HACK FLO
  
  if (length(xdat) >= GLOBAL_min_years_data) {  
    # z <- list()
    npmu <- length(mul) + 1
    npsc <- length(sigl) + 1
    npsh <- length(shl) + 1
    z$trans <- FALSE
    dat.Lmom <- Lmoments(xdat)
    par.init<-par.genlogis(dat.Lmom[1],dat.Lmom[2], dat.Lmom[4])
    in2 <- par.init[2]
    in1 <- par.init[1]
    
    if (is.null(mul)) {
      mumat <- as.matrix(rep(1, length(xdat)))
      
      if (is.null(muinit)) 
        muinit <- in1
    } else {
      z$trans <- TRUE
      mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
      
      if (is.null(muinit)) 
        muinit <- c(in1, rep(0, length(mul)))
    }
    if (is.null(sigl)) {
      sigmat <- as.matrix(rep(1, length(xdat)))
      
      if (is.null(siginit)) 
        siginit <- in2
    } else {
      z$trans <- TRUE
      sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
      
      if (is.null(siginit)) 
        siginit <- c(in2, rep(0, length(sigl)))
    }
    if (is.null(shl)) {
      shmat <- as.matrix(rep(1, length(xdat)))
      
      if (is.null(shinit)) 
        shinit <- 0.1
    } else {
      z$trans <- TRUE
      shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
      
      if (is.null(shinit)) 
        shinit <- c(0.1, rep(0, length(shl)))
    }
    
    z$model <- list(mul, sigl, shl)
    z$link <- deparse(substitute(c(mulink, siglink, shlink)))
    init <- c(muinit, siginit, shinit)
    gl.lik <- function(a) {
      mu <- mulink(mumat %*% (a[1:npmu]))
      sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
      xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
      y <- (xdat - mu) / sc
      y <- 1 - xi * y
      if (any(y <= 0) || any(sc <= 0)) 
        return(10^6)
      sum(log(sc)) - sum(log(y^((1/xi)-1)))+2*sum(log((1+y^(1/xi)))) 
    }
    
    #   x <- optim(init, gl.lik, hessian = TRUE, method = method,  # Original code commented by FKB
    #       control = list(maxit = maxit, ...))
    
    fail_safe <- failwith(NULL, optim)  # START FKB HACK
    x <- fail_safe(init, gl.lik, hessian = TRUE, method = method, control = list(maxit = maxit, ...))
    
    if (is.null(x) == TRUE) {
      print("Warning: The function optim failed in gl_mle")
      invisible(param)
    } else {
      ## End hack FKB
      
      
      z$conv <- x$convergence
      mu <- mulink(mumat %*% (x$par[1:npmu]))
      sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
      xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
      z$nllh <- x$value
      z$dat <- xdat
      if (z$trans) {
        z$dat <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
      }
      # z$mle replaced by z$estimate to integrate with the other functions
      z$estimate <- x$par
      
      
      
      # z$cov <- solve(x$hessian)  # initially this way
      # z$cov <- tryCatch(solve(x$hessian), finally = "Warning: could not solve Hessian")  # Protection FLO 
      # z$cov <- try(solve(x$hessian))  # TO FIX!!!!!!!!!!!!!!!!
      # z$se <- sqrt(diag(z$cov))  # initially this way
      # z$se <- try(sqrt(diag(z$cov)))  # TO FIX!!!!!!!!!!!!!!!!
      
      
      z$vals <- cbind(mu, sc, xi)
      if (show) {
        if (z$trans) 
          print(z[c(2, 3, 4)])
        else print(z[4])
        if (!z$conv) 
          print(z[c(5, 7, 9)])
      }
      # class(z) <- "gev.fit"  # Commented by FLO
      invisible(z)
    } 
  }
  else {
    print(paste("Warning: this station has less than ",GLOBAL_min_years_data," years of data. Use another method!", 
                collapse = "",sep = ""))  
    invisible(param)  
  }
}

#' Fitting the GL distribution with Lmom
#' @description Function to fit the Generalized Logistics distribution with the linear moments method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gl_Lmom(XXX)
gl_Lmom <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA)) 
  if (length(dat) >= GLOBAL_min_years_data) {

    dat.Lmom <- Lmoments(dat)
    # ADD FAILSAFE?
    fitted.param <- par.genlogis(dat.Lmom[1], dat.Lmom[2], dat.Lmom[4])
    # Creating the returning list. Standard error is not yet implemented
    param$estimate <- c(fitted.param$xi, fitted.param$alfa, fitted.param$k)
    return(param)
    } else {
      print(paste("Warning: This station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse = "",sep = ""))  
      invisible(param)  
    }
}

#' Fitting the GL distribution with mom
#' @description Function to fit the Generalized Logistics distribution with the ordinary moments method
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se).
#' Standard error is not yet implemented
#' @export
#'
#' @examples gl_mom(XXX)
gl_mom <- function(dat) {

  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA))  # HACK FLO
  if (length(dat) >= GLOBAL_min_years_data) {
    alpha1 <-  moments(dat)[1]
    alpha2 <- moments(dat)[2]
    alpha3 <- moments(dat)[4]
    fitted.param <- c()
    # This is an approximation as the initial value of the optim
    x0 <- 2 / (3 * pi) * atan( - 0.59484 * alpha3)
    findshape <- function(k) {
    g1 <- gamma(1 + k) * gamma(1 - k)
    g2 <- gamma(1 + 2*k) * gamma(1 - 2*k)
    g3 <- gamma(1 + 3*k) * gamma(1 - 3*k)  
    sign(k) * (3*g2*g1 - g3 - 2*g1^3) / ((g2 - g1^2)^(1.5)) - alpha3
  }  
   
  # rfind <- newtonRaphson(findshape, x0, maxiter = 100, tol = 0.000001) # Original code
	## Start Hack FKB
  fail_safe <- failwith(NULL, newtonRaphson)  # START FKB HACK
  rfind <- fail_safe(findshape, x0, maxiter = 100, tol = 0.000001)
    
  if (is.null(rfind) == TRUE) {
    print("Warning: The function newtonRaphson failed in gl_mom")
    invisible(param)
  } else {
  ## End hack FKB
    
  fitted.param[3] <-  rfind$root
  g1 <- gamma(1 + fitted.param[3]) * gamma(1 - fitted.param[3])
  g2 <- gamma(1 + 2*fitted.param[3]) * gamma(1 - 2*fitted.param[3])
	
  fitted.param[2] <- abs(fitted.param[3]) * alpha2 / sqrt(g2-g1^2)
  fitted.param[1] <- alpha1 - fitted.param[2] / fitted.param[3] * (1 - g1) 
  
  # Creating the returning list. Standard error is not yet implemented
  param$estimate <- c(fitted.param[1], fitted.param[2], fitted.param[3])
  invisible(param)
  }
  } else {
    print(paste("Warning: this station has less than ",GLOBAL_min_years_data," years of data. Use another method!", 
                  collapse="",sep=""))  
    invisible(param)  
    }
}

#' Fitting the GL distribution with Bayesian inference
#' @description Function to fit the Generalized Logistics distribution with BayesianMCMC method
#' We assume that the shape parameter only has a prior with mean zero and standard deviation 0.2 (dnorm(x[3], 0, 0.2))
#' @param dat the data that needs fitting (i.e. flood data)
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples gl_bayes(XXX)
gl_bayes <- function(dat) {
# Fit GL distribution with the Bayesian method  
  param <- list(estimate = c(NA, NA, NA), se = c(NA, NA, NA)) 
  
  if (length(dat) >= GLOBAL_min_years_data) {
    # Prior for Bayes
    myprior <- function (x) {
      # x = vector of parameter values: c(location, scale, shape)
      # I assume the shape parameter only has a prior with mean zero and standard deviation 0.2
      dnorm(x[3], 0, 0.2)
    }
    
    fail_safe <- failwith(NULL, BayesianMCMC)
    bayes <- fail_safe(dat, nbpas = 5000, nbchaines = 3, 
                       confint = c(0.05, 0.95), dist = "GENLOGIS", apriori = myprior)

    if (is.null(bayes) == TRUE) {
      print("Warning: the function BayesianMCMC failed in gl_bayes")
      invisible(param)
    } else {
      
      ## Addition to return parameters
      #   # Solution 1
      param$estimate <- bayes$parametersML
      
#       # Solution 2
#       param$estimate[1] <- mean(as.vector(bayes$parameters[, 1, 1:3]))
#       param$estimate[2] <- mean(as.vector(bayes$parameters[, 2, 1:3]))
#       param$estimate[3] <- mean(as.vector(bayes$parameters[, 3, 1:3]))
      
      param$se[1] <- sd(as.vector(bayes$parameters[, 1, 1:3]))
      param$se[2] <- sd(as.vector(bayes$parameters[, 2, 1:3]))
      param$se[3] <- sd(as.vector(bayes$parameters[, 3, 1:3]))
      
      invisible(param)
    }
  } else {
    print(paste("Warning: this station has less than ", GLOBAL_min_years_data," years of data. Use another method!", 
                collapse = "", sep = ""))  
    invisible(param)  
  }
}

#' Calculating the posterior predictive distribution for GL
#' @description Function to calculate the posterior predictive distribution after calling gev_bayes
#' @param (mmrp, mupars, spars, kpars) parameters returned by gev_bayes. mupars, spars, kpars are the ensemble of param$estimate
#' @return param Estimated parameters and standard error returned as a list($estimate, $se)
#' @export
#'
#' @examples Needs example
get_posterior_gl <- function(mmrp, mupars, spars, kpars) {

  qqsample1 <- sapply(seq(length(mupars)), function(st) {
    mean_temp <- mupars[st]
    st_temp <-spars[st]
    k_temp <- kpars[st]
    invF.genlogis(F = (1 - 1 / mmrp), mean_temp, st_temp, k_temp)
    }, simplify = "array")
  # 1 is for collums only 0.5 to return the median
  qqr <- apply(qqsample1, 1, quantile, 0.5)
}
