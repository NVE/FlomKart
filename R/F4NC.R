# F4NC.r
# Supporting GOF functions to create the NetCDF bigmama file!

#' gof_cs
#' @description Chi Square calculation. Does not work well
#' @param dat
#' @param param
#' @param distr
#' @importFrom dplyr failwith
#' @return
#' @export
#' @import nsRFA
#' @examples
gof_cs <- function(dat, param, distr = "distr") {

  # Parametrize the number of bins
  nb.bins = trunc(1.88 * length(dat)^(2/5) / 2)   # 1/2 of Formula from http://kb.palisade.com/index.php?pg=kb.page&id=57
  xmin <- min(dat)
  xmax <- max(dat)
  xdiff <- xmax - xmin

  bin.size <- xdiff / nb.bins
  dat.breaks <- seq(xmin, xmax, bin.size)
  observed.binned <- table(cut(dat, breaks = dat.breaks))  # binning dat table
  p <- c()
  CS <- NA

  # This depends on the distribution
  if (distr == 'gumbel') {
    for (i in 1:nb.bins) {
      fail_safe <- failwith(NA, evd::pgumbel)
      temp <- fail_safe(dat.breaks[i+1], param[1], param[2]) - fail_safe(dat.breaks[i], param[1], param[2])
      p <- append(p, temp)
    }
  }
  if (distr == 'gamma') {
    for (i in 1:nb.bins) {
      fail_safe <- failwith(NA, pgamma)
      temp <- fail_safe(dat.breaks[i+1], param[1], rate = param[2]) - fail_safe(dat.breaks[i], param[1], rate = param[2])
      p <- append(p, temp)
    }
  }
  if (distr == 'gev') {
    for (i in 1:nb.bins) {
      fail_safe <- failwith(NA, evd::pgev)  # evd::pgev  # also tried nsRFA::F.GEV
      temp <- fail_safe(dat.breaks[i+1], param[1], param[2], param[3]) - fail_safe(dat.breaks[i], param[1], param[2], param[3])
      p <- append(p, temp)
    }
  }
  if (distr == 'gl') {
    for (i in 1:nb.bins) {
      fail_safe <- failwith(NA, nsRFA::F.genlogis)
      temp <- fail_safe(dat.breaks[i+1], param[1], param[2], param[3]) - fail_safe(dat.breaks[i], param[1], param[2], param[3])
      p <- append(p, temp)
    }
  }
  if (distr == 'pearson') {
    for (i in 1:nb.bins) {
      fail_safe <- failwith(NA, nsRFA::F.gamma)
      temp <- fail_safe(dat.breaks[i+1], param[1], param[2], param[3]) - fail_safe(dat.breaks[i], param[1], param[2], param[3])
      p <- append(p, temp)
    }
  }
  if (any(is.na(p)) == FALSE && any(is.nan(p)) == FALSE && length(p) == nb.bins) {
    # A condition on sum(p) >0.99 could be added
  N <- length(dat)
  temp <- 0

      for (i in 1:nb.bins) {

      if (p[i] < 0.001) {  # This condition is to avoid CS becoming unstable (very high values)
        next(i)
      } else {
        temp = temp + N * p[i] * ( (observed.binned[i] / N - p[i]) / p[i] )^2

        }
    }
    CS <- temp
  }
  # We could add if/else (is.na(CS)) {print("Warning: gof_cs has failed with distr...)}
  invisible(CS)
}


#' gof_ks
#' @description Kolmogorov Smirnov. Gives results with warnings but this is because some flood data are repeated (for some strange reason)
#' @param dat
#' @param param
#' @param distr
#' @param test.stat
#' @param p.value
#'
#' @return
#' @export
#' @import goftest tidyverse evd nsRFA glogis fBasics
#' @examples
gof_ks <- function(dat, param, distr = "distr", test.stat = TRUE , p.value = FALSE) {

  KS <- NA
  fail_safe <- failwith(NA, stats::ks.test)

  if (distr == 'gumbel') {
    temp <- fail_safe(dat, "pgumbel", param[1], param[2])
  }
  if (distr == 'exp') {
    temp <- fail_safe(dat, "F.exp", param[1], param[2])
  }
  if (distr == 'gamma') {
    temp <- fail_safe(dat, "pgamma", param[1], rate = param[2])
  }
  if (distr == 'gev') {
    temp <- fail_safe(dat, "pgev", param[1], param[2], param[3])
  }
  if (distr == 'gl') {
    temp <- fail_safe(dat, "F.genlogis", param[1], param[2], param[3])
  }
  if (distr == 'gp') {
    temp <- fail_safe(dat, "F.genpar", param[1], param[2], param[3])
  }
  if (distr == 'pearson') {
    temp <- fail_safe(dat, "F.gamma", param[1], param[2], param[3])
  }

  if (p.value == TRUE && is.list(temp) == TRUE) {
    KS <- temp$p.value
  }
  else if (test.stat == TRUE && is.list(temp) == TRUE && is.numeric(temp$statistic) == TRUE && !is.infinite(temp$statistic) == TRUE) {
    KS <- temp$statistic
  }
  # We could add if (is.na(KS)) {print("Warning: gof_ks has failed with distr...)}
  invisible(KS)
}


#' gof_ad
#' @description Anderson Darling
#' @param dat
#' @param param
#' @param distr
#' @param test.stat
#' @param p.value
#' @importFrom goftest ad.test
#' @return
#' @export
#' @import goftest tidyverse evd nsRFA glogis fBasics
#' @examples
gof_ad <- function(dat, param, distr = "distr", test.stat = TRUE , p.value = FALSE) {

  AD <- NA
  fail_safe <- failwith(NA, goftest::ad.test)

  if (distr == 'gumbel') {
    temp <- fail_safe(dat, "pgumbel", param[1], param[2])
  }
  if (distr == 'exp') {
    dat <- dat[dat > param[1]]
    temp <- fail_safe(dat, "F.exp", param[1], param[2])
  }
  if (distr == 'gamma') {
    temp <- fail_safe(dat, "pgamma", param[1], rate = param[2])
  }
  if (distr == 'gev') {
    temp <- fail_safe(dat, "pgev", param[1], param[2], param[3])
  }
  if (distr == 'gl') {
    temp <- fail_safe(dat, "F.genlogis", param[1], param[2], param[3])
  }
  if (distr == 'gp') {
    dat <- dat[dat > param[1]]
    temp <- fail_safe(dat, "F.genpar", param[1], param[2], param[3])
  } 
  if (distr == 'pearson') {
    temp <- fail_safe(dat, "F.gamma", param[1], param[2], param[3])
  }

  if (p.value == TRUE && is.list(temp) == TRUE) {
    AD <- temp$p.value
  } else if (test.stat == TRUE && is.list(temp) == TRUE && is.numeric(temp$statistic) == TRUE && !is.infinite(temp$statistic) == TRUE) {
    AD <- temp$statistic
  }
  # We could add if (is.na(AD)) {print("Warning: gof_ad has failed with distr...)}
  invisible(AD)
}

#' QS4NC
#' @description Quantile Score. Only uses dat.sample so the calling script should be sending DAT.SAMPLE!
#' This is implemented according to Equation 2 of the "Cross validation study and verification document"
#' @param dat
#' @param r.levels
#' @param r.periods
#' @import nsRFA
#' @return
#' @export
#'
#' @examples
QS4NC <- function(dat, r.levels, r.periods) {

  r.prob <-  1 - 1 / r.periods
  QS <- array(NA, dim = length(r.levels))

    for (z in seq(along = r.levels)) {

      binary_vector <- as.integer(dat <= r.levels[z])
      QS[z] <- mean( (dat - r.levels[z]) * (r.prob[z] - binary_vector) )
    }

    invisible(QS)
}

#' BS4NC
#' @description Brier Score.
#' This is implemented according to Equation 1 of the "Cross validation study and verification document"
#' @param dat
#' @param threshold
#' @param param
#' @param distr
#' @import nsRFA
#' @return
#' @export
#'
#' @examples
BS4NC <- function(dat, threshold, param, distr = "distr") {

  BS <- array(NA, dim = 6)
  modelled.prob <- array(NA, dim = 6)
  empirical <- array(NA, dim = 6)
  Pu <- array(NA, dim = 6)

#   m <- max(dat.full) * 0.8
#   threshold <- seq(m / 6, m, m / 6)

  if(distr == 'gumbel') {
    modelled.prob <- pgumbel(threshold, param[1], param[2])  # Could be protected with "failwith"
  }
  if(distr == 'gamma') {
    modelled.prob <- pgamma(threshold, param[1], rate = param[2])
  }
  if(distr == 'gev') {
    modelled.prob <- evd::pgev(threshold, param[1], param[2], param[3])
  }
  if(distr == 'gl') {
    modelled.prob <- F.genlogis(threshold, param[1], param[2], param[3])
  }
  if(distr == 'pearson') {
    modelled.prob <- F.gamma(threshold, param[1], param[2], param[3])
  }

  for (z in 1:6) {
    binary_vector <- as.integer(dat >= threshold[z]) # / length(dat) # * 30 otherwise, very small values
    Pu[z] <- (1 - modelled.prob[z])
    BS[z] <- mean( (Pu[z] - binary_vector)^2 )
  }
invisible(BS)

}

############################################################################
# QS4NC_OLD <- function(dat, dat.full, param, distr = "distr") {
#   # Quantile Score. Only uses dat.sample so the calling script should be sending DAT.SAMPLE!
#   # This is implemented according to Equation 2 of the "Cross validation study and verification document"
#
#   QS <- array(NA, dim = 6)
#   modelled.prob <- array(NA, dim = 6)
#
#   m <- max(dat.full) * 0.8
#   threshold <- seq(m / 6, m, m / 6)
#
#   if(distr == 'gumbel') {
#     modelled.prob <- pgumbel(threshold, param[1], param[2])
#   }
#   if(distr == 'gamma') {
#     modelled.prob <- pgamma(threshold, param[1], rate = param[2])
#   }
#   if(distr == 'gev') {
#     modelled.prob <- evd::pgev(threshold, param[1], param[2], param[3])
#   }
#   if(distr == 'gl') {
#     modelled.prob <- F.genlogis(threshold, param[1], param[2], param[3])
#   }
#   if(distr == 'pearson') {
#     modelled.prob <- F.gamma(threshold, param[1], param[2], param[3])
#   }
#
#   for (z in seq(along = threshold)) {
#
#     index1 <- which(dat <= threshold[z])
#     index2 <- which(dat > threshold[z])
#
#     QS[z] <- sum((dat[index1] - threshold[z])) * (modelled.prob[z] - 1) / length(dat)
#     + sum((dat[index2] - threshold[z])) * modelled.prob[z] / length(dat)
# #     QS[z] <- sum((dat[index1] - threshold[z])) * (modelled.prob[z] - 1 / length(dat))
# #     + sum((dat[index2] - threshold[z])) * modelled.prob[z]
#
#   }
#   invisible(QS)
# }

# BS4NC_OLD <- function(dat, param, r.periods, distr = "distr") {
#   # Brier Score.
#   # This is implemented according to Equation 1 of the "Cross validation study and verification document" TO DOUBLE CHECK
#
#   BS <- array(NA, dim = length(r.periods))
#   r.prob <- 1 - 1 / r.periods
#   # empirical.cdf <- c()
#   modelled.cdf <- c()
#
#   if(distr == 'gumbel') {
#     modelled.cdf <- pgumbel(dat, param[1], param[2])  # modelled.cdf.sample?
#   }
#   if(distr == 'gamma') {
#     modelled.cdf <- pgamma(dat, param[1], rate = param[2])
#   }
#   if(distr == 'gev') {
#     modelled.cdf <- evd::pgev(dat, param[1], param[2], param[3])
#   }
#   if(distr == 'gl') {
#     modelled.cdf <- nsRFA::F.genlogis(dat, param[1], param[2], param[3])
#   }
#   if(distr == 'pearson') {
#     modelled.cdf <- nsRFA::F.gamma(dat, param[1], param[2], param[3])
#   }
#   empirical.cdf <- ecdf(dat)
#
#   # Fill vectors with max and min divergence (density function or return level)
#   # Computing return levels for the specified return periods
#   for (z in seq(along = r.periods)) {
#
#     empirical <- length(which(empirical.cdf(dat) >= r.prob[z])) / length(dat) * 30  # * 30 otherwise, very small values
#     modelled <- length(which(modelled.cdf >= r.prob[z])) / length(dat) * 30
#     BS[z] <- (modelled - empirical)^2
#   }
#   invisible(BS)
#
# }
##############################################################################


#' RLEVELS4NC
#'
#' @param param.estimate
#' @param return.periods
#' @param distr
#'
#' @return
#' @export
#'
#' @examples
RLEVELS4NC <- function(param.estimate, return.periods, distr = "distr") {

  r.levels <- array(NA, dim = length(return.periods))
  fail_vector <- array(NA, dim = length(return.periods))
  if(distr == 'gumbel') {
    fail_safe <- failwith(fail_vector, invF.gumb)
    r.levels <- fail_safe(F = (1 - 1 / return.periods), param.estimate[ 1], param.estimate[2])
    if (is.na(r.levels[1]) == TRUE) { print("Warning: the function invF.gumb failed in RLEVELS4NC") }
  }
  if(distr == 'gamma') {
    fail_safe <- failwith(fail_vector, qgamma)
    r.levels <- fail_safe(p = (1 - 1 / return.periods), shape = param.estimate[1], scale = 1 / param.estimate[2])
    # qgamma comes from Rlab package? Be careful!!!!!!!!
    if (is.na(r.levels[1]) == TRUE) { print("Warning: the function qgamma failed in RLEVELS4NC") }
  }
  if(distr == 'gev') {
    fail_safe <- failwith(fail_vector, evd::qgev)
    r.levels <- fail_safe( (1 - 1 / return.periods),
                          param.estimate[1], param.estimate[2], param.estimate[3])
    if (is.na(r.levels[1]) == TRUE) { print("Warning: the function evd::qgev failed in RLEVELS4NC") }
  }
  if(distr == 'gl') {
    fail_safe <- failwith(fail_vector, nsRFA::invF.genlogis)
    r.levels <- fail_safe(F = (1 - 1 / return.periods),
                          param.estimate[1], param.estimate[2], param.estimate[3])
    if (is.na(r.levels[1]) == TRUE) { print("Warning: the function invF.genlogis failed in RLEVELS4NC") }
  }
  if(distr == 'pearson') {
    fail_safe <- failwith(fail_vector, nsRFA::invF.gamma)
    r.levels <- fail_safe(F = (1 - 1 / return.periods),
                          param.estimate[1], param.estimate[2], param.estimate[3])
    if (is.na(r.levels[1]) == TRUE) { print("Warning: the function invF.gamma failed in RLEVELS4NC") }
  }
  invisible(r.levels)
}



#' gof_nt
#' @description Function calculating the NT reliability measure  TO FINISH
#' @param dat
#' @param rperiods.nt
#' @param param
#' @param distr
#'
#' @return
#' @export
#'
#' @examples
gof_nt <- function(dat, rperiods.nt, param, distr = "distr") {
# st is for station, rs is for subsample, threshold[k] is for quantile (could be indexed)

  ## Calculation of the predicted quantile thresholds.
   if(distr == 'gumbel')  thresholds <- evd::qgumbel(1 - 1 / rperiods.nt, param[1], param[2])
   if(distr == 'gamma') thresholds <- qgamma(1 - 1 / rperiods.nt, param[1], scale = 1 / param.estimate[2])

     if(distr == 'gev') thresholds <- evd::qgev(1 - 1 / rperiods.nt, param[1], param[2], param[3])

     if(distr == 'gl')     thresholds <- invF.genlogis(1 - 1 / rperiods.nt, param[1], param[2], param[3])
     if(distr == 'pearson') thresholds <- invF.gamma(1 - 1 / rperiods.nt, param[1], param[2], param[3])

  NT <- array(NA, dim = length(rperiods.nt))

 nn_mean <- 0
   n1 <- 0
   n2 <- 0

  for (k in 1:6) {

    # NT[k] <- sum(dat > thresholds[k])

    nn_mean <- sum(dat > thresholds[k])
     n1 <- pbinom(nn_mean, size = length(dat), prob = 1 / rperiods.nt[k])  # 50 is the number of tries
     if (nn_mean == 0) {
       n2 <- 0

     } else {
       n2 <- pbinom((nn_mean - 1), size = length(dat), prob = 1 / rperiods.nt[k])
     }
     NT[k] <- runif(1, n2, n1)

  }

  invisible(NT)
}






##############################################################################################
#

#' run_all_ff
#' @description FF function copied from run_all.r
#' Test of run_all for ff only.
#' @param dat
#' @param station.nb.vect
#' @param distr
#' @param method
#'
#' @return
#' @export
#'
#' @examples
run_all_ff <- function(dat,station.nb.vect, distr = "distr", method = "method") {

  st.name <- c()
  ff.1 <- c()
  ff.2 <- c()

  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])
    m <- max(sdat$flom_DOGN)

    if (is.null(sdat) == TRUE) {
      print(i)
      station.nb.vect[i] <- NA
      st.name[i] <- "NO DATA"
      ff.1[i] <- NA
      ff.2[i] <- NA

      next(i)
    } else if (length(sdat$flom_DOGN) < min_years_data) {
      st.name[i] <- "not enough data"
      ff.1[i] <- NA
      ff.2[i] <- NA

      next(i)
    } else {
      st.name[i] <- as.character(sdat$name[1])
      if (is.na(st.name[i]) == TRUE) { st.name[i] <- "NO_st.name" }
      # Creation of the split sample
      L <- length(sdat$flom_DOGN)
      #       # Historical split
      #       split <- trunc(length(sdat$flom_DOGN) / 2)
      #       sample.1 <- sdat$flom_DOGN[1:split]
      #       sample.2 <- sdat$flom_DOGN[(split + 1):L]

      # Random split MAYBE PARAMETRIZE THE SAMPLING SIZE
      indices.full <- seq(1, L, 1)
      indices.1 <- sample.int(L, size = L /2, replace = FALSE)
      indices.2 <- setdiff(indices.full, indices.1)
      sample.1 <- sdat$flom_DOGN[indices.1]
      sample.2 <- sdat$flom_DOGN[indices.2]

      #########################
      if ( distr == "gamma") {
        # GAMMA -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("gamma_", as.character(method), collapse = "", sep = ""))
        fail_safe <- failwith(NULL, FUN)
        param.1 <- fail_safe(sample.1)
        param.2 <- fail_safe(sample.2)
        sample.1.max <- max(sample.1)
        sample.2.max <- max(sample.2)

        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgamma)
          ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[1], rate = param.1$estimate[2])^length(sample.2)
        }
        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgamma)
          ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[1], rate = param.2$estimate[2])^length(sample.1)
        }
      }
      #########################
      if ( distr == "gumbel") {
        # GUMBEL -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("gumbel_", as.character(method), collapse = "", sep = ""))
        fail_safe <- failwith(NULL, FUN)
        param.1 <- fail_safe(sample.1)
        param.2 <- fail_safe(sample.2)
        sample.1.max <- max(sample.1)
        sample.2.max <- max(sample.2)

        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgumbel)
          ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[1], param.1$estimate[2])^length(sample.2)
        }
        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgumbel)
          ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[1], param.2$estimate[2])^length(sample.1)
        }
      }
      ############
      if ( distr == "gev") {
        # GEV -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("gev_", as.character(method), collapse = "", sep = ""))
        fail_safe <- failwith(NULL, FUN)
        param.1 <- fail_safe(sample.1)
        param.2 <- fail_safe(sample.2)
        sample.1.max <- max(sample.1)
        sample.2.max <- max(sample.2)

        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgev)
          ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[3], param.1$estimate[1], param.1$estimate[2])^length(sample.2)
        }
        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, pgev)
          ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[3], param.2$estimate[1], param.2$estimate[2])^length(sample.1)
        }
      }
      ############
      if ( distr == "gl") {
        # GL -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("gl_", as.character(method), collapse = "", sep = ""))
        fail_safe <- failwith(NULL, FUN)
        param.1 <- fail_safe(sample.1)
        param.2 <- fail_safe(sample.2)
        sample.1.max <- max(sample.1)
        sample.2.max <- max(sample.2)

        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, F.genlogis)
          ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[1], param.1$estimate[2], param.1$estimate[3])^length(sample.2)
        }
        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, F.genlogis)
          ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[1], param.2$estimate[2], param.2$estimate[3])^length(sample.1)
        }
      }
      ############
      if ( distr == "pearson") {
        # PEARSON -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("pearson_", as.character(method), collapse = "", sep = ""))
        fail_safe <- failwith(NULL, FUN)
        param.1 <- fail_safe(sample.1)
        param.2 <- fail_safe(sample.2)
        sample.1.max <- max(sample.1)
        sample.2.max <- max(sample.2)

        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, F.gamma)
          ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[1], param.1$estimate[2], param.1$estimate[3])^length(sample.2)
        }
        if (is.null(param.1) == FALSE) {
          fail_safe <- failwith(NA, F.gamma)
          ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[1], param.2$estimate[2], param.2$estimate[3])^length(sample.1)
        }
      }
      ##########################
    }  # end of big "else"
  }  # end of for loop on stations

  x <- list(snumber = station.nb.vect, name = st.name, ff1 = ff.1, ff2 = ff.2)
  index.1 <- which(!is.na(x$ff1))
  print(index.1)
  # #   index.2 <- intersect(x, which(!is.na(x$ff2)))
  # #   print(index.2)
  x1 <- x$ff1[index.1]
  x2 <- x$ff2[index.1]

  write.table(x, file = "results_ff.csv", sep = ";", col.names = NA, qmethod = "double")
  invisible(data.frame(ff1 = x1, ff2 = x2))
}


# ### COPIED FROM GOF.R AND MAIN.R TO USE THE FF FUNCTION
#
# distr <- "pearson"
# method <- "mle"
# ff <- run_all_ff(dat, station.nb.vect, distr, method)
# sorted.ff.1 <- sort(ff$ff1)
# sorted.ff.2 <- sort(ff$ff2)
# uniform.1 <- seq(0, 1, length.out = length(sorted.ff.1))
# uniform.2 <- seq(0, 1, length.out = length(sorted.ff.2))
# area.1 <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
# area.2 <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
# area.1.abs  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
# area.2.abs <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
# plot(ecdf(ff$ff1))
# plot(ecdf(ff$ff2))
#
# test <- gof_area(dat, station.nb.vect)


# gof_area <- function(dat, station.nb.vect) {
#   # Calculation of the area between the uniform distribution ECDF and the FF distribution ECDDF
#   # According to Blanchet.
#
#   distributions = c("gumbel","gamma","gev","gl","pearson")
#   fitting_methods = c("mle","Lmom","mom")
#   method = "mom"
#   i <- 1
#   for ( distr in distributions ) {
#     # for (method in fitting_methods) {
#     print( paste( "fitting parameters for ", distr))
#
#     ff <- run_all_ff(dat, station.nb.vect, distr, method)
#     sorted.ff.1 <- sort(ff$ff1)
#     sorted.ff.2 <- sort(ff$ff2)
#     uniform.1 <- seq(0, 1, length.out = length(sorted.ff.1))
#     uniform.2 <- seq(0, 1, length.out = length(sorted.ff.2))
#
#     #       area.1$distr <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
#     #       area.2$distr <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
#     #       area.1.abs$distr  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
#     #       area.2.abs$distr <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
#     area.1[i] <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
#     area.2[i] <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
#     area.1.abs[i]  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
#     area.2.abs[i] <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
#     i <- i + 1
#     # }
#   }
#   results <- list(area1 = area.1, area1abs = area.1.abs,  area2 = area.2, area2abs = area.2.abs)
#   return(results)
# }


# gof_ff <- function(st,rs, k) {
#   # st is for station, rs is for subsample, threshold[k] is for quantile (could be indexed)
#   nn_mean <- sum(evaluationdata > threshold[k])
#   n1 <- pbinom(nn_mean, qll, 1 / MyRP_NT[threshold[k]])
#   if (nn_mean == 0) {
#     n2 <- 0
#
#   } else {
#     n2 <- pbinom((nn_mean - 1), qll, 1 / MyRP_NT[threshold[k]])
#     NTvalue_mean_nx[st, rs, k] <- runif(1, n2, n1)
#
#   }
#   invisible(NTvalue_mean_nx)
# }

