# This file creates the NetCDF file containing all GOF analysis

# setwd('C:/Users/flbk/Documents/GitHub/FlomKart/R')  # CHECK DIR

# Load libraries
# library(RNetCDF)
# library(goftest)      # Goodness of fit tests (ad.test)
# library(plyr)         # For the failwith(NA, f) function!

# library(evd)          # Functions for extreme value distributions  library(evir)         # Extreme values in R
# library(fBasics)      # Some basic statistical functions
# library(glogis)       # Generalized logistics
# library(nsRFA)        # invF.gumb()
# library(stats4)       # for mle()
# library(verification) # for the quantile score function
# library(PearsonDS)    # Pearson distribution
# library(ismev)        # An Introduction to Statistical Modeling of Extreme Values

# Source all the required R files
# source('F4NC.R')      # Goodness of fit functions
# source('LOAD.r')

################################### BIG FUNCTION TO PARAMETRIZE

#' update_gofnc
#' @description Function to update the nc database
#' @return
#' @export
#'
#' @examples
update_gofnc <- function(flood_database.nc_path = "output/flood_database.nc") {

## To run at every update --------------
nc <- open.nc(paste(flood_database.nc_path), write = FALSE)  # CHECK DIR
# READ ONLY: we want to load all pre-existing dimensions, parameters and data

# New dimensions
dim.r_periods <- 6  # This is new to GOF (wasn't in NC_DATABASE)
return.periods <- c(10, 20, 50, 100, 200, 500)  # This is new to GOF (wasn't in NC_DATABASE)
rperiods.bs <- c(2,5,10,15,20,30)  # return periods for the Brier score
dim.few_quantiles <- 7  # min/max/q25%/q50%/Q75%/mean/sd

# Create a vector with all the different station numbers, in order to run every station separately in one go
station.nb.vect <- var.get.nc(nc, "station.number")
dim.station <- length(station.nb.vect)

distr.name <- var.get.nc(nc, "distr.name")
method.name <- var.get.nc(nc, "method.name")
dim.distr <- length(distr.name)
dim.method <- length(method.name)
dim.length_rec <- var.get.nc(nc, "dim.length_rec")

dim.param <- 3
dim.characters <- 64

dim.random_runs <- var.get.nc(nc, "dim.random_runs")
min_years_data <- var.get.nc(nc, "min_years_data")
sampling_years <- var.get.nc(nc, "sampling_years")
dim.max_subsample <- max(sampling_years)

}


#' create_empty_gofnc
#' @description Function to create an empty gof_nc database
#' @return
#' @export
#'
#' @examples
create_empty_gofnc <- function(flood_database.nc_path = "output/flood_database.nc") {

  ## To run if creating the nc file from scratch  --------------
  # Get important data and parameters from the already compiled flood_database.nc

  nc <- open.nc(paste(flood_database.nc_path), write = FALSE)
  # READ ONLY: we want to load all pre-existing dimensions, parameters and data

  # New dimensions
  dim.r_periods <- 6  # This is new to GOF (wasn't in NC_DATABASE)
  return.periods <- c(10, 20, 50, 100, 200, 500)  # This is new to GOF (wasn't in NC_DATABASE)
  rperiods.bs <- c(2,5,10,15,20,30)  # return periods for the Brier score
  dim.few_quantiles <- 7  # min/max/q25%/q50%/Q75%/mean/sd

  # Create a vector with all the different station numbers, in order to run every station separately in one go
  station.nb.vect <- var.get.nc(nc, "station.number")
  dim.station <- length(station.nb.vect)

  distr.name <- var.get.nc(nc, "distr.name")
  method.name <- var.get.nc(nc, "method.name")
  dim.distr <- length(distr.name)
  dim.method <- length(method.name)
  dim.length_rec <- var.get.nc(nc, "dim.length_rec")

  dim.param <- 3
  dim.characters <- 64

  dim.random_runs <- var.get.nc(nc, "dim.random_runs")
  min_years_data <- var.get.nc(nc, "min_years_data")
  sampling_years <- var.get.nc(nc, "sampling_years")
  dim.max_subsample <- max(sampling_years)


# Now we create the gof.nc files, define dimensions and start filling it with empty data

gof_nc <- create.nc("output/gof.nc")
att.put.nc(gof_nc, "NC_GLOBAL", "title", "NC_CHAR", "Flood frequency analysis results")
att.put.nc(gof_nc, "NC_GLOBAL", "history", "NC_CHAR", paste("Created on", base::date()))

# Defining dimensions in netCDF
dim.def.nc(gof_nc, "station", dim.station)
dim.def.nc(gof_nc, "distr", dim.distr)
dim.def.nc(gof_nc, "method", dim.method)
dim.def.nc(gof_nc, "param", dim.param)
dim.def.nc(gof_nc, "length.rec", dim.length_rec)
dim.def.nc(gof_nc, "few_quantiles", dim.few_quantiles)
dim.def.nc(gof_nc, "max_string_length", dim.characters)
dim.def.nc(gof_nc, "r.periods", dim.r_periods)

# Variable for the return periods
var.def.nc(gof_nc, varname = "r.periods", vartype = "NC_SHORT",
           dimensions = "r.periods")
att.put.nc(gof_nc, "r.periods", "missing_value", "NC_SHORT", -9999)
att.put.nc(gof_nc, "r.periods", "short_name", "NC_CHAR", "Design return periods")
att.put.nc(gof_nc, "r.periods", "long_name", "NC_CHAR",
           "The following return periods were used: 10, 20, 50, 100, 200 and 500 years")
var.put.nc(gof_nc, "r.periods", return.periods)

var.def.nc(gof_nc, varname = "r.periods.qs", vartype = "NC_SHORT",
           dimensions = "r.periods")
att.put.nc(gof_nc, "r.periods.qs", "missing_value", "NC_SHORT", -9999)
att.put.nc(gof_nc, "r.periods.qs", "short_name", "NC_CHAR", "Design return periods")
att.put.nc(gof_nc, "r.periods.qs", "long_name", "NC_CHAR",
           "The following return periods were used: 2, 5, 10, 15, 20 and 30 years")
var.put.nc(gof_nc, "r.periods.qs", rperiods.bs)

var.def.nc(gof_nc, varname = "r.periods.bs", vartype = "NC_SHORT",
           dimensions = "r.periods")
att.put.nc(gof_nc, "r.periods.bs", "missing_value", "NC_SHORT", -9999)
att.put.nc(gof_nc, "r.periods.bs", "short_name", "NC_CHAR", "Design return periods")
att.put.nc(gof_nc, "r.periods.bs", "long_name", "NC_CHAR",
           "The following return periods were used: 2, 5, 10, 15, 20 and 30 years")
var.put.nc(gof_nc, "r.periods.bs", rperiods.bs)

# Variable for the return levels and various GOF
var.def.nc(gof_nc, varname = "r.levels", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles", "r.periods"))
att.put.nc(gof_nc, "r.levels", "missing_value", "NC_FLOAT", -9999)
att.put.nc(gof_nc, "r.levels", "short_name", "NC_CHAR", "Return levels (m3/s) for different return periods")
att.put.nc(gof_nc, "r.levels", "long_name", "NC_CHAR",
           "The selected return periods for those return levels are: 10, 20, 50, 100, 200 and 500 years")

r.levels <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles, dim.r_periods))
var.put.nc(gof_nc, "r.levels", r.levels)
rm(r.levels)

var.def.nc(gof_nc, varname = "QS", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles", "r.periods"))
att.put.nc(gof_nc, "QS", "missing_value", "NC_FLOAT", -9999)

QS <- array(NA,dim = c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles, dim.r_periods))
var.put.nc(gof_nc, "QS", QS)
rm(QS)

var.def.nc(gof_nc, varname = "BS", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles", "r.periods"))
att.put.nc(gof_nc, "BS", "missing_value", "NC_FLOAT", -9999)

BS <- array(NA,dim = c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles, dim.r_periods))
var.put.nc(gof_nc, "BS", BS)
rm(BS)

# var.def.nc(gof_nc, varname = "NT", vartype = "NC_FLOAT",
#            dimensions = c("station", "distr", "method", "length.rec", "r.periods"))
# att.put.nc(gof_nc, "NT", "missing_value", "NC_FLOAT", -9999)
# NT <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.length_rec, dim.r_periods))
# var.put.nc(gof_nc, "NT", NT)
# rm(NT)

var.def.nc(gof_nc, varname = "CS", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles"))
att.put.nc(gof_nc, "CS", "missing_value", "NC_FLOAT", -9999)

CS <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles))
var.put.nc(gof_nc, "CS", CS)
rm(CS)

var.def.nc(gof_nc, varname = "KS", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles"))
att.put.nc(gof_nc, "KS", "missing_value", "NC_FLOAT", -9999)

KS <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles))
var.put.nc(gof_nc, "KS", KS)
rm(KS)

var.def.nc(gof_nc, varname = "AD", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "length.rec", "few_quantiles"))
att.put.nc(gof_nc, "AD", "missing_value", "NC_FLOAT", -9999)

AD <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.length_rec, dim.few_quantiles))
var.put.nc(gof_nc, "AD", AD)
rm(AD)

sync.nc(gof_nc)

}


#' fillup_gofnc
#' @description Function to fill up the gof_nc database
#' @return
#' @export
#' @import goftest tidyverse fitdistrib evd nsRFA ismev fitdistrplus stats
#' @examples
fillup_gofnc <- function() {

  library(goftest)
  library(tidyverse)
  library(fitdistrib)
  library(evd)
  library(nsRFA)
  library(ismev)
  library(fitdistrplus)
  library(stats)

## To run if updating or creating the nc file from scratch  --------------

nc <- open.nc("output/flood_database.nc", write = FALSE)
gof_nc <- open.nc("output/gof.nc", write = TRUE)  # Put FALSE for read-only
Q <- var.get.nc(nc, "Q")
dim.length_rec <- var.get.nc(nc, "dim.length_rec") # To specify the end of the vector where the full record will be stored

# Create a vector with all the different station numbers, in order to run every station separately in one go
station.nb.vect <- var.get.nc(nc, "station.number")
dim.station <- length(station.nb.vect)

distr.name <- var.get.nc(nc, "distr.name")
method.name <- var.get.nc(nc, "method.name")
dim.distr <- length(distr.name)
dim.method <- length(method.name)

dim.random_runs <- var.get.nc(nc, "dim.random_runs")
min_years_data <- var.get.nc(nc, "min_years_data")
sampling_years <- var.get.nc(nc, "sampling_years")
dim.max_subsample <- max(sampling_years)

## TO TIDY UP so that it is read from the nc databases
# New dimensions
dim.r_periods <- 6  # This is new to GOF (wasn't in NC_DATABASE)
return.periods <- c(10, 20, 50, 100, 200, 500)  # This is new to GOF (wasn't in NC_DATABASE)
rperiods.bs <- c(2,5,10,15,20,30)  # return periods for the Brier score
dim.few_quantiles <- 7  # min/max/q25%/q50%/Q75%/mean/sd

dim.param <- 3
dim.characters <- 64


# Dumping all console output into errorlog.txt
sink("output/errorlog_gofnc.txt")

for (st in seq(along = station.nb.vect)) {
# for (st in 1:3) {
  print(st)
  temp.Q <- as.vector(na.omit(Q[st, ]))
  if (length(temp.Q) <  min_years_data) {  # This is not GLOBAL_min_years to make sure "min_years_data"
    #has been loaded from the nc database
    print("This station number has NULL or not enough data")
    print(station.nb.vect[st])
    next(st)
  } else {

    temp.random_indexes <- var.get.nc(nc, "random_indexes", start=c(st, 1, 1, 1),
               count=c(1, dim.random_runs, length(sampling_years), dim.max_subsample))

    ## Creation of the empirical thresholds. This could be moved to NC_DATABASE.R and saved in a new variable
    ranked.Q <- rank(temp.Q, ties.method = "random")
    empp = (ranked.Q-0.4)/(length(temp.Q)+0.2)
    efunc <- approxfun(empp, temp.Q)
    thresholds <- efunc(1-1/rperiods.bs)
    print(thresholds)


    for (d in 1:5) {
      distr <- distr.name[d]
      print(distr)
      for (m in 1:3) {
        method <- method.name[m]
        print(method)
        print("Computing full dataset")

        # Temporary parameters retrieved from the parameter database. A 1d array
        # reminder: param.estimate[st, d, m, p, dim.length_rec, rs]
        param.estimate <- var.get.nc(nc, "param.estimate",
                                     start = c(st, d, m, 1, dim.length_rec, 1), count = c(1, 1, 1, 3, 1, 1))

        # Create temporary variables with NAs to fill the GOF_NC database
        # reminder r.levels[st, d, m, length_rec, dim.few_quantiles, dim.r_periods]
        temp.r.levels <- rep(NA, dim.r_periods)
        temp.QS <- rep(NA, dim.r_periods)
        temp.BS <- rep(NA, dim.r_periods)
        # temp.NT <- rep(NA, dim.r_periods)


        temp.CS <- NA
        temp.KS <- NA
        temp.AD <- NA
        # Evaluate temporary variables
        if (length(param.estimate) > 1) {
        temp.r.levels <- RLEVELS4NC(param.estimate, return.periods, distr = distr)
        if (!all(is.na(temp.r.levels))) {
          # We use the small return periods for QS and BS (via thresholds for BS)
        temp.QS <- QS4NC(temp.Q, temp.r.levels, rperiods.bs)
        temp.BS <- BS4NC(temp.Q, thresholds, param.estimate, distr = distr)
        # temp.NT <- gof_nt(temp.Q, rperiods.bs, param.estimate, distr = distr)  # TOCHECK

        temp.CS <- gof_cs(temp.Q, param.estimate, as.character(distr))
        temp.KS <- gof_ks(temp.Q, param.estimate, as.character(distr))
        temp.AD <- gof_ad(temp.Q, param.estimate, as.character(distr))
        }
        }
        # WRITE TEMP VARIABLES TO NC on spot dim.length_rec of the length_rec dimension
        # reminder CS[st, d, m, length_rec, dim.few_quantiles]
        var.put.nc(gof_nc, "r.levels", data = temp.r.levels, start = c(st, d, m, dim.length_rec, 1, 1), count = c(1, 1, 1, 1, 1, dim.r_periods))
        var.put.nc(gof_nc, "QS", data = temp.QS, start = c(st, d, m, dim.length_rec, 1, 1), count = c(1, 1, 1, 1, 1, dim.r_periods))
        var.put.nc(gof_nc, "BS", data = temp.BS, start = c(st, d, m, dim.length_rec, 1, 1), count = c(1, 1, 1, 1, 1, dim.r_periods))
        # var.put.nc(gof_nc, "NT", data = temp.NT, start = c(st, d, m, dim.length_rec, 1), count = c(1, 1, 1, 1, dim.r_periods))  # TOCHECK


        var.put.nc(gof_nc, "CS", data = temp.CS, start = c(st, d, m, dim.length_rec, 1), count = c(1, 1, 1, 1, 1))
        var.put.nc(gof_nc, "KS", data = temp.KS, start = c(st, d, m, dim.length_rec, 1), count = c(1, 1, 1, 1, 1))
        var.put.nc(gof_nc, "AD", data = temp.AD, start = c(st, d, m, dim.length_rec, 1), count = c(1, 1, 1, 1, 1))
        # Recreate temporary variables but with 2 extra dimensions for j and rs
        # reminder: param.estimate[st, d, m, p, dim.length_rec, rs]
        rm(param.estimate)
        param.sample.estimate <- var.get.nc(nc, "param.estimate",
                                     start = c(st, d, m, 1, 1, 1),
                                     count = c(1, 1, 1, 3, length(sampling_years), dim.random_runs))

        temp.r.levels <- array(NA,dim=c(length(sampling_years), dim.few_quantiles, dim.r_periods))
        temp.QS <- array(NA,dim=c(length(sampling_years), dim.few_quantiles, dim.r_periods))
        temp.BS <- array(NA,dim=c(length(sampling_years), dim.few_quantiles, dim.r_periods))
        # temp.NT <- array(NA,dim=c(length(sampling_years), dim.r_periods))  # TOCHECK

        temp.CS <- array(NA,dim=c(length(sampling_years), dim.few_quantiles))
        temp.KS <- array(NA,dim=c(length(sampling_years), dim.few_quantiles))
        temp.AD <- array(NA,dim=c(length(sampling_years), dim.few_quantiles))

        # RANDOM START: Fill up with the random runs with sampled data TO DO
        # j is the index from 1 to 23 and sampling_years[j] is 30, 35... 90
        print("Computing sampled dataset")
        if (m < 4) {  # Random sampling of Bayes takes a very long time
        j <- 1
        while (sampling_years[j] <= max(sampling_years)) {
          # sampling_years[j] <= length(temp.Q) && ## This was in the while before which stops it to early for stations with short records
          temp.random.rlevels <- array(NA,dim=c(dim.random_runs, dim.r_periods))
          temp.random.QS <- array(NA,dim=c(dim.random_runs, dim.r_periods))
          temp.random.BS <- array(NA,dim=c(dim.random_runs, dim.r_periods))
          # temp.random.NT <- array(NA,dim=c(dim.random_runs, dim.r_periods))  # TOCHECK

          temp.random.CS <- rep(NA, dim.random_runs)
          temp.random.KS <- rep(NA, dim.random_runs)
          temp.random.AD <- rep(NA, dim.random_runs)

          for (rr in 1:dim.random_runs) {
            # For goodness of fit evaluation, we attribute to sample.Q another random run
            if (rr == dim.random_runs) { sample.q <- temp.Q[temp.random_indexes[1, j, 1:sampling_years[j]]] }

            else { sample.q <- temp.Q[temp.random_indexes[rr + 1, j, 1:sampling_years[j]]] }

            # Could be changed to assess GOF on same random sample
            # sample.q <- temp.Q[temp.random_indexes[rr, j, 1:sampling_years[j]]]


            if (length(param.sample.estimate[1:3, j, rr]) > 1) {
            temp.random.rlevels[rr, 1:dim.r_periods] <- RLEVELS4NC(param.sample.estimate[1:3, j, rr], return.periods, distr = distr)

            # modified by KOE
			temp.random.rlevels.qs <- RLEVELS4NC(param.sample.estimate[1:3, j, rr], rperiods.bs, distr = distr)
#####################################################################################################
            if (!all(is.na(temp.random.rlevels[rr, ]))) {
# modified by KOE
            temp.random.QS[rr, 1:dim.r_periods] <- QS4NC(temp.Q, temp.random.rlevels.qs, rperiods.bs)
#            temp.random.QS[rr, 1:dim.r_periods] <- QS4NC(sample.q, temp.random.rlevels[rr, 1:dim.r_periods], return.periods)
#########################################################################
            temp.random.BS[rr, 1:dim.r_periods] <- BS4NC(temp.Q, thresholds, param.sample.estimate[1:3, j, rr], distr = distr)
#########################################################################################################

            }
            temp.random.CS[rr] <- gof_cs(temp.Q, param.sample.estimate[1:3, j, rr], as.character(distr))  # I think it makes sense this way, but
            temp.random.KS[rr] <- gof_ks(temp.Q, param.sample.estimate[1:3, j, rr], as.character(distr))  # it could be interesting to try sending the sample.q
            temp.random.AD[rr] <- gof_ad(temp.Q, param.sample.estimate[1:3, j, rr], as.character(distr))
            }
          }
          # Calculating and filling the quantiles we want before moving on to the next record length sample
          # We have 7 numbers: mean, sd, min, max, Q25, Q50, Q75
          for (rp in 1:dim.r_periods) {
            if (all(is.na(temp.random.rlevels[1:dim.random_runs, rp])) == FALSE) {
              temp.r.levels[j, 1, rp] <- mean(temp.random.rlevels[1:dim.random_runs, rp], na.rm = TRUE)
              temp.r.levels[j, 2, rp] <- sd(temp.random.rlevels[1:dim.random_runs, rp], na.rm = TRUE)
              temp.r.levels[j, 3, rp] <- min(temp.random.rlevels[1:dim.random_runs, rp], na.rm = TRUE)
              temp.r.levels[j, 4, rp] <- max(temp.random.rlevels[1:dim.random_runs, rp], na.rm = TRUE)
              temp.r.levels[j, 5, rp] <- quantile(temp.random.rlevels[1:dim.random_runs, rp], probs = 0.25, na.rm = TRUE)
              temp.r.levels[j, 6, rp] <- median(temp.random.rlevels[1:dim.random_runs, rp], na.rm = TRUE)
              temp.r.levels[j, 7, rp] <- quantile(temp.random.rlevels[1:dim.random_runs, rp], probs = 0.75, na.rm = TRUE)
            }

            if (all(is.na(temp.random.QS[1:dim.random_runs, rp])) == FALSE) {
              temp.QS[j, 1, rp] <- mean(temp.random.QS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.QS[j, 2, rp] <- sd(temp.random.QS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.QS[j, 3, rp] <- min(temp.random.QS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.QS[j, 4, rp] <- max(temp.random.QS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.QS[j, 5, rp] <- quantile(temp.random.QS[1:dim.random_runs, rp], probs = 0.25, na.rm = TRUE)
              temp.QS[j, 6, rp] <- median(temp.random.QS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.QS[j, 7, rp] <- quantile(temp.random.QS[1:dim.random_runs, rp], probs = 0.75, na.rm = TRUE)
            }

            if (all(is.na(temp.random.BS[1:dim.random_runs, rp])) == FALSE) {
              temp.BS[j, 1, rp] <- mean(temp.random.BS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.BS[j, 2, rp] <- sd(temp.random.BS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.BS[j, 3, rp] <- min(temp.random.BS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.BS[j, 4, rp] <- max(temp.random.BS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.BS[j, 5, rp] <- quantile(temp.random.BS[1:dim.random_runs, rp], probs = 0.25, na.rm = TRUE)
              temp.BS[j, 6, rp] <- median(temp.random.BS[1:dim.random_runs, rp], na.rm = TRUE)
              temp.BS[j, 7, rp] <- quantile(temp.random.BS[1:dim.random_runs, rp], probs = 0.75, na.rm = TRUE)
            }


            ## Process the NT values for each subsample into an absolute area
#             if (all(is.na(temp.random.NT[1:dim.random_runs, rp])) == FALSE) {
#
#
#               sorted.nt <- sort(temp.random.NT[1:dim.random_runs, rp])
#               uniform.seq <- seq(0, 1, length.out = length(sorted.nt))
#               temp.NT[j, rp]  <- sum(abs(sorted.nt - uniform.seq)) / 2 / length(sorted.nt)
#             }

          }

          if (!all(is.na(temp.random.CS))) {
            temp.CS[j, 1] <- mean(na.omit(temp.random.CS))
            temp.CS[j, 2] <- sd(na.omit(temp.random.CS))
            temp.CS[j, 3] <- min(na.omit(temp.random.CS))
            temp.CS[j, 4] <- max(na.omit(temp.random.CS))
            temp.CS[j, 5] <- quantile(na.omit(temp.random.CS), probs = 0.25)
            temp.CS[j, 6] <- median(na.omit(temp.random.CS))
            temp.CS[j, 7] <- quantile(na.omit(temp.random.CS), probs = 0.75)

            }
          if (!all(is.na(temp.random.KS))) {
            temp.KS[j, 1] <- mean(na.omit(temp.random.KS))
            temp.KS[j, 2] <- sd(na.omit(temp.random.KS))
            temp.KS[j, 3] <- min(na.omit(temp.random.KS))
            temp.KS[j, 4] <- max(na.omit(temp.random.KS))
            temp.KS[j, 5] <- quantile(na.omit(temp.random.KS), probs = 0.25)
            temp.KS[j, 6] <- median(na.omit(temp.random.KS))
            temp.KS[j, 7] <- quantile(na.omit(temp.random.KS), probs = 0.75)

            }
          if (!all(is.na(temp.random.AD))) {
            temp.AD[j, 1] <- mean(na.omit(temp.random.AD))
            temp.AD[j, 2] <- sd(na.omit(temp.random.AD))
            temp.AD[j, 3] <- min(na.omit(temp.random.AD))
            temp.AD[j, 4] <- max(na.omit(temp.random.AD))
            temp.AD[j, 5] <- quantile(na.omit(temp.random.AD), probs = 0.25)
            temp.AD[j, 6] <- median(na.omit(temp.random.AD))
            temp.AD[j, 7] <- quantile(na.omit(temp.random.AD), probs = 0.75)

            }

          rm(temp.random.rlevels)
          rm(temp.random.QS)
          rm(temp.random.BS)
          # rm(temp.random.NT)

          rm(temp.random.CS)
          rm(temp.random.KS)
          rm(temp.random.AD)

          if (j == length(sampling_years)) { break }
          else { j <- j + 1 }
        }
      }  # CLOSING IF FOR THE EXCLUSION OF BAYES FOR RANDOM RUNS
        # WRITE TO NC
        # reminder r.levels[st, d, m, length_rec, dim.few_quantiles, dim.r_periods]
        var.put.nc(gof_nc, "r.levels", data = temp.r.levels, start = c(st, d, m, 1, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles, dim.r_periods))
        var.put.nc(gof_nc, "QS", data = temp.QS, start = c(st, d, m, 1, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles, dim.r_periods))
        var.put.nc(gof_nc, "BS", data = temp.BS, start = c(st, d, m, 1, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles, dim.r_periods))
#         var.put.nc(gof_nc, "NT", data = temp.NT, start = c(st, d, m, 1, 1),
#                    count = c(1, 1, 1, length(sampling_years), dim.r_periods))  # TOCHECK

        var.put.nc(gof_nc, "CS", data = temp.CS, start = c(st, d, m, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles))
        var.put.nc(gof_nc, "KS", data = temp.KS, start = c(st, d, m, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles))
        var.put.nc(gof_nc, "AD", data = temp.AD, start = c(st, d, m, 1, 1),
                   count = c(1, 1, 1, length(sampling_years), dim.few_quantiles))

        rm(temp.r.levels)
        rm(temp.QS)
        rm(temp.BS)
        # rm(temp.NT)

        rm(temp.CS)
        rm(temp.KS)
        rm(temp.AD)
      }
    }
  }
  rm(temp.Q)
  rm(temp.random_indexes)
}
sync.nc(gof_nc)
sink()  # Stop printing console onto file


}  # end of fillup_gofnc function

## FINISHED -> CHECK RESULTS ON PANOPLY OR SHINY APP
