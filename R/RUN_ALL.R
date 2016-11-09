# RUN_ALL.R
# This file gathers all the function that analyze all the stations at once


#' fit_distr
#' @description The function that calls the right fitting function for each distribution and fitting method
#' @param dat
#' @param distr
#' @param method
#'
#' @return
#' @export
#'
#' @examples
fit_distr <- function(dat, distr = "distr", method = "method") {

  # Paste the function call and assign it to a temporary function FUN
  FUN <- match.fun(paste(as.character(distr), "_", as.character(method), collapse = "",sep = ""))
  FUN(dat)
}


#' run_all
#' @description The function to run all distributions and all methods and output csv files of the GOF
#' @param dat
#' @param station.nb.vect
#' @param distr
#' @param method
#'
#' @return
#' @export
#'
#' @examples
run_all <- function(dat, station.nb.vect, distr = "distr", method = "method") {

  st.name <- c()
  p1.val <- c()
  p1.se <- c()
  p2.val <- c()
  p2.se <- c()
  p3.val <- c()
  p3.se <- c()
  cs <- c()
  ks <- c()
  ad <- c()

  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])

    if (is.null(sdat) == TRUE) {
      print(i)
      station.nb.vect[i] <- NA
      next(i)
    } else {
      st.name <- append(st.name, as.character(sdat$name[1]))
      print(sdat$name[1])
      print(i)
      # Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste(as.character(distr), "_", as.character(method), collapse = "", sep = ""))
      # FUN(sdat$flom_DOGN)
      fail_safe <- failwith(NA, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.na(param) == FALSE) {
        p.values <- gof_all(sdat$flom_DOGN, param$estimate, as.character(distr))
        cs <- append(cs, round(p.values$CS, 2))
        ks <- append(ks, round(p.values$KS, 2))
        ad <- append(ad, round(p.values$AD, 2))

        param$estimate <- round(param$estimate, 2)
        param$se <- round(param$se, 2)
        p1.val <- append(p1.val, param$estimate[1])
        p1.se <- append(p1.se, param$se[1])
        p2.val <- append(p2.val, param$estimate[2])
        p2.se <- append(p2.se, param$se[2])
        p3.val <- append(p3.val, param$estimate[3])
        p3.se <- append(p3.se, param$se[3])
        } else {
          p1.val <- append(p1.val, NA)
          p1.se <- append(p1.se, NA)
          p2.val <- append(p2.val, NA)
          p2.se <- append(p2.se, NA)
          p3.val <- append(p3.val, NA)
          p3.se <- append(p3.se, NA)
          cs <- append(cs, NA)
          ks <- append(ks, NA)
          ad <- append(ad, NA)
        }
      }

  }
  x <- list(snumber = na.omit(station.nb.vect), name = st.name, p1_value = p1.val , p1_std.err = p1.se,
                   p2_value = p2.val , p2_std.err = p2.se,
                   p3_value = p3.val , p3_std.err = p3.se,
                   chi_square = cs, kolmogorov = ks, anderson= ad)

  write.table(x, file = paste("results_", as.character(distr), ".csv", collapse = "", sep = ""), sep = ";",
              col.names = NA, qmethod = "double")
  invisible(x)
}


#' run_all_norm
#' @description The function to run all distributions and all methods and output csv files of the GOF
#' @param dat
#' @param station.nb.vect
#' @param distr
#' @param method
#'
#' @return
#' @export
#'
#' @examples
run_all_norm <- function(dat,station.nb.vect, distr = "distr", method = "method") {

# What's the difference with run_all

  st.name <- c()
  p1.val <- c()
  p1.se <- c()
  p2.val <- c()
  p2.se <- c()
  p3.val <- c()
  p3.se <- c()
  cs <- c()
  ks <- c()
  ad <- c()

  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])

    if (is.null(sdat) == TRUE) {
      station.nb.vect[i] <- NA
      next(i)
    } else {
      st.name <- append(st.name, as.character(sdat$name[1]))
      # Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste(as.character(distr), "_", as.character(method), collapse = "", sep = ""))
      FUN(sdat$normQ)

      if (is.na(param) == FALSE) {
        p.values <- gof_all(sdat$normQ, param$estimate, as.character(distr))
        cs <- append(cs, round(p.values$CS, 2))
        ks <- append(ks, round(p.values$KS, 2))  # failwith is missing here, which leads to loop-stoping errors
        ad <- append(ad, round(p.values$AD, 2))  # To think carefully, this whole function could be redesigned...

        param$estimate <- round(param$estimate, 2)
        param$se <- round(param$se, 2)
        p1.val <- append(p1.val, param$estimate[1])
        p1.se <- append(p1.se, param$se[1])
        p2.val <- append(p2.val, param$estimate[2])
        p2.se <- append(p2.se, param$se[2])
        p3.val <- append(p3.val, param$estimate[3])
        p3.se <- append(p3.se, param$se[3])
      } else {
          p1.val <- append(p1.val, NA)
          p1.se <- append(p1.se, NA)
          p2.val <- append(p2.val, NA)
          p2.se <- append(p2.se, NA)
          p3.val <- append(p3.val, NA)
          p3.se <- append(p3.se, NA)
          cs <- append(cs, NA)
          ks <- append(ks, NA)
          ad <- append(ad, NA)
        }
      }

  }
  x <- list(snumber = na.omit(station.nb.vect), name = st.name, p1_value = p1.val , p1_std.err = p1.se,
            p2_value = p2.val , p2_std.err = p2.se,
            p3_value = p3.val , p3_std.err = p3.se,
            chi_square = cs, kolmogorov = ks, anderson= ad)

  write.table(x, file = paste("results_norm_", as.character(distr), ".csv", collapse = "", sep = ""), sep = ";",
              col.names = NA, qmethod = "double")
  invisible(x)
}


#' run_all_ks
#' @description Test of run_all for Ks only.
#' @param dat
#' @param station.nb.vect
#' @param method
#'
#' @return
#' @export
#'
#' @examples
run_all_ks <- function(dat,station.nb.vect, method = "method") {

  st.name <- c()
  ks <- c()
  distrib <- c()
  param <- list()


  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])

    if (is.null(sdat) == TRUE) {
      print(i)
      station.nb.vect[i] <- NA
      st.name[i] <- "NO DATA"
      ks[i] <- NA
      distrib[i] <- NA
      next(i)
    } else if (length(sdat$flom_DOGN) < min_years_3param) {
      st.name[i] <- as.character(sdat$name[1])
      ks[i] <- NA
      distrib[i] <- "not enough data"
      next(i)
      } else {
          st.name[i] <- as.character(sdat$name[1])
          if (is.na(st.name[i]) == TRUE) st.name[i] <- "NO_st.name"

          # GAMMA -> Paste the function call and assign it to a temporary function FUN
        FUN <- match.fun(paste("gamma_", as.character(method), collapse = "", sep = ""))
        # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
        # param <- try(FUN(sdat$flom_DOGN))
        fail_safe <- failwith(NULL, FUN)
        param <- fail_safe(sdat$flom_DOGN)

        if (is.null(param) == FALSE) {
          fail_safe <- failwith(NULL, gof_ks)
          KS <- fail_safe(sdat$flom_DOGN, param$estimate, "gamma")
          print(KS)
          if (is.null(KS) == TRUE) {
            ks[i] <- NA
            distrib[i] <- "KS_GAMMA_FAILED"
          } else if (is.na(KS$p.value) == TRUE) {
              ks[i] <- NA
              distrib[i] <- "KS_GAMMA_FAILED"
            } else {
                fail_safe <- failwith(0, round)
                ks[i] <- fail_safe(KS$p.value,2)
                distrib[i] <- "GAMMA"
            }
        }

## TIDY UP GAMMA COMPLETELY BEFORE THE OTHER DISTRIBS

      # GUMBEL -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gumbel_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ks)
        KS <- fail_safe(sdat$flom_DOGN, param$estimate, "gumbel")
        print(KS)
        if (is.null(KS) == TRUE) {
          ks[i] <- NA
          distrib[i] <- "KS_GUMBEL_FAILED"
        } else if (is.na(KS$p.value) == TRUE) {
          print("jesus...")
          ks[i] <- NA
          distrib[i] <- "KS_GUMBEL_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(KS$p.value,2)
          print("temp")
          print(temp)
          print(ks[i])
          if (temp > ks[i] | is.na(ks[i]) == TRUE) {
            ks[i] <- temp
            distrib[i] <- "GUMBEL"
            print("GUMBEL")
          }
        }
      }
      # GEV -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gev_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ks)
        KS <- fail_safe(sdat$flom_DOGN, param$estimate, "gev")
        if (is.null(KS) == TRUE) {
          ks[i] <- NA
          distrib[i] <- "KS_GEV_FAILED"
        } else if (is.na(KS$p.value) == TRUE) {
          print("jesus...")
          ks[i] <- NA
          distrib[i] <- "KS_GEV_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(KS$p.value,2)
          print("temp")
          print(temp)
          print(ks[i])
          if (temp > ks[i] | is.na(ks[i]) == TRUE) {
            ks[i] <- temp
            distrib[i] <- "GEV"
            print("GEV")
          }
        }
      }
      # GL -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gl_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)
      print(param)
      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ks)
        KS <- fail_safe(sdat$flom_DOGN, param$estimate, "gl")
        print(KS)
        if (is.null(KS) == TRUE) {
          ks[i] <- NA
          distrib[i] <- "KS_GL_FAILED"
        } else if (is.na(KS$p.value) == TRUE) {
          print("jesus...")
          ks[i] <- NA
          distrib[i] <- "KS_GL_FAILED"
        } else {
        fail_safe <- failwith(0, round)
        temp <- fail_safe(KS$p.value,2)
        print("temp")
        print(temp)
        print(ks[i])
          if (temp > ks[i] | is.na(ks[i]) == TRUE) {
          ks[i] <- temp
          distrib[i] <- "GL"
          print("GL")
          }
        }
      }
      # PEARSON -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("pearson_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ks)
        KS <- fail_safe(sdat$flom_DOGN, param$estimate, "pearson")
        print(KS)
        if (is.null(KS) == TRUE) {
          ks[i] <- NA
          distrib[i] <- "KS_PEARSON_FAILED"
        } else if (is.na(KS$p.value) == TRUE) {
          print("jesus...")
          ks[i] <- NA
          distrib[i] <- "KS_PEARSON_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(KS$p.value,2)
          print("temp")
          print(temp)
          print(ks[i])
          if (temp > ks[i] | is.na(ks[i]) == TRUE) {
            ks[i] <- temp
            distrib[i] <- "PEARSON"
            print("PEARSON")
          }
        }
      }


    }
  }

   print(length(na.omit(station.nb.vect)))
   print(length(na.omit(st.name)))
   print(length(na.omit(ks)))
   print(length(na.omit(distrib)))

   print(length(station.nb.vect))
   print(length(st.name))
   print(length(ks))
   print(length(distrib))

  x <- list(snumber = station.nb.vect, name = st.name, gof = ks, distrib = distrib)
  write.table(x, file = "results.csv", sep = ";", col.names = NA, qmethod = "double")
  invisible(x)
}


#' run_all_ad
#' @description TEST OF RUN ALL FOR AD ONLY
#' @param dat
#' @param station.nb.vect
#' @param method
#'
#' @return
#' @export
#'
#' @examples
run_all_ad <- function(dat,station.nb.vect, method = "method") {

  st.name <- c()
  ad <- c()
  distrib <- c()
  param <- list()


  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])

    if (is.null(sdat) == TRUE) {
      print(i)
      station.nb.vect[i] <- NA
      st.name[i] <- "NO DATA"
      ad[i] <- NA
      distrib[i] <- NA
      next(i)
    } else if (length(sdat$flom_DOGN) < min_years_3param) {
      st.name[i] <- as.character(sdat$name[1])
      ad[i] <- NA
      distrib[i] <- "not enough data"
      next(i)
    } else {
      st.name[i] <- as.character(sdat$name[1])
      if (is.na(st.name[i]) == TRUE) st.name[i] <- "NO_st.name"

      print(sdat$name[1])
      print(i)

      # GAMMA -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gamma_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ad)
        AD <- fail_safe(sdat$flom_DOGN, param$estimate, "gamma")
        print(AD)
        if (is.null(AD) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_GAMMA_FAILED"
        } else if (is.na(AD$p.value) == TRUE) {
          print("jesus...")
          ad[i] <- NA
          distrib[i] <- "AD_GAMMA_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          ad[i] <- fail_safe(AD$p.value,2)
          distrib[i] <- "GAMMA"
        }
      } else {
        ad[i] <- NA
        distrib[i] <- "GAMMA_OPTIM_FAILED"
      }

      # GUMBEL -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gumbel_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ad)
        AD <- fail_safe(sdat$flom_DOGN, param$estimate, "gumbel")
        print(AD)
        if (is.null(AD) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_GUMBEL_FAILED"
        } else if (is.na(AD$p.value) == TRUE) {
          print("jesus...")
          ad[i] <- NA
          distrib[i] <- "AD_GUMBEL_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(AD$p.value,2)
          print("temp")
          print(temp)
          print(ad[i])
          if (temp > ad[i] | is.na(ad[i]) == TRUE) {
            ad[i] <- temp
            distrib[i] <- "GUMBEL"
            print("GUMBEL")
          }
        }
      }
      # GEV -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gev_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = param <- NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ad)
        AD <- fail_safe(sdat$flom_DOGN, param$estimate, "gev")
        if (is.null(AD) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_GEV_FAILED"
        } else if (is.na(AD$p.value) == TRUE) {
          print("jesus...")
          ad[i] <- NA
          distrib[i] <- "AD_GEV_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(AD$p.value,2)
          print("temp")
          print(temp)
          print(ad[i])
          if (temp > ad[i] | is.na(ad[i]) == TRUE) {
            ad[i] <- temp
            distrib[i] <- "GEV"
            print("GEV")
          }
        }
      }
      # GL -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("gl_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)
      print(param)
      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ad)
        AD <- fail_safe(sdat$flom_DOGN, param$estimate, "gl")
        print(AD)
        if (is.null(AD) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_GL_FAILED"
        } else if (is.na(AD$p.value) == TRUE) {
          print("jesus...")
          ad[i] <- NA
          distrib[i] <- "AD_GL_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(AD$p.value,2)
          print("temp")
          print(temp)
          print(ad[i])
          if (temp > ad[i] | is.na(ad[i]) == TRUE) {
            ad[i] <- temp
            distrib[i] <- "GL"
            print("GL")
          }
        }
      }
      # PEARSON -> Paste the function call and assign it to a temporary function FUN
      FUN <- match.fun(paste("pearson_", as.character(method), collapse = "", sep = ""))
      # param <- tryCatch(FUN(sdat$flom_DOGN), finally = NA)
      # param <- try(FUN(sdat$flom_DOGN))
      fail_safe <- failwith(NULL, FUN)
      param <- fail_safe(sdat$flom_DOGN)

      if (is.null(param) == FALSE) {
        fail_safe <- failwith(NULL, gof_ad)
        AD <- fail_safe(sdat$flom_DOGN, param$estimate, "pearson")
        print(AD)
        if (is.null(AD) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_PEARSON_FAILED"
        } else if (is.na(AD$p.value) == TRUE) {
          ad[i] <- NA
          distrib[i] <- "AD_PEARSON_FAILED"
        } else {
          fail_safe <- failwith(0, round)
          temp <- fail_safe(AD$p.value,2)
          if (temp > ad[i] | is.na(ad[i]) == TRUE) {
            ad[i] <- temp
            distrib[i] <- "PEARSON"
            print("PEARSON")
          }
        }
      }


    }
  }

  print(length(na.omit(station.nb.vect)))
  print(length(na.omit(st.name)))
  print(length(na.omit(ad)))
  print(length(na.omit(distrib)))

  print(length(station.nb.vect))
  print(length(st.name))
  print(length(ad))
  print(length(distrib))

  x <- list(snumber = station.nb.vect, name = st.name, gof = ad, distrib = distrib)
  write.table(x, file = "results.csv", sep = ";", col.names = NA, qmethod = "double")
  invisible(x)
}


################################################################################################
################################################################################################
################################################################################################



#' run_all_ff2
#' @description Test of run_all for ff only.
#' TO COMPARE WITH run_all_ff in F4NCE.R!!!
#' @param dat
#' @param station.nb.vect
#' @param distr
#' @param method
#' @importFrom dplyr failwith
#' @return
#' @export
#'
#' @examples
run_all_ff2 <- function(dat,station.nb.vect, distr = "distr", method = "method") {

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
      ff.1[i] <- fail_safe(sample.2.max, param.1$estimate[1], param.1$estimate[2])^length(sample.2)
    }
    if (is.null(param.1) == FALSE) {
      fail_safe <- failwith(NA, pgamma)
      ff.2[i] <- fail_safe(sample.1.max, param.2$estimate[1], param.2$estimate[2])^length(sample.1)
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

