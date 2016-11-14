# Creation of the NetCDF file for the whole project
# The design of the NetCDF file is as follows:
# Dimensions: Stations
#             Distributions
#             Estimation methods
#             Length of record
#             Optionally: number of random runs
# Variables:
#   - First the variables depending only on the station:
#       Flom_DOGN
#       Year
#       Time
#       Longitude
#       Latitude
#       Catchment size...
#   - Variables depending on the station, the distribution, the method and the length of record
#       Parameter estimate
#       Standard error of estimate
#       Goodness of fit parameters: CS, KS, AD, BS, QS, FF
#       Return levels for 10, 20,50, 100, 200, 500  year events
#   - Then optionally, the variables depending on the lenght of the record and the number of random runs:
#       index of selected data points

# Clean workspace, load libraries and source scripts  ---------------------------------

## In order to get the fitdistr package
# packages <- c("devtools","quadprog")
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages, rownames(installed.packages())))
# }
# library(devtools)
#
# if (length(setdiff("fitdistrib", rownames(installed.packages()))) > 0) {
#   install_github("fbaffie/fitdistrib")
# }

# rm(list = ls())  # clean out workspace and set working directory

# Load libraries
# library(RNetCDF)      # To work with NetCDF files
# library(fitdistrib) # My own package for fitting distributions
# library(lubridate)
# library(parallel) # to detect number of cores
# library(foreach) # for parallel for loops




################################### BIG FUNCTION TO PARAMETRIZE




#' update_nc
#' @description Function to update the nc database
#' @return
#' @export
#'
#' @examples
update_nc <- function() {


mysystem<-Sys.info()['sysname']
if(mysystem == "Linux") library(doMC) # parallel backend Linux
if(mysystem == "Windows") library(doSNOW) # parallel backend Windows

# Source all the required R files
# source('LOAD.r')      # Required to load Q into database and save station coordinates


## To run at every update --------------
raw_dat <- read.delim("../rawdata/AMS_table_HYDRA_endelig.txt", sep="\t", as.is=TRUE)  # CHECK DIR
raw_dat$year <- year(raw_dat$dt_flood_daily_date)
raw_dat$snumber <- raw_dat$regine_area * 10^5 + raw_dat$main_nr

# Old data for station names only
old_dat <- read.csv("../rawdata/AMS_table_old.csv", sep=";", as.is=TRUE)  # CHECK DIR

for (i in 1:length(raw_dat$year)) {
  if (length(which(old_dat$snumber == raw_dat$snumber[i])) > 0) {
  raw_dat$name[i] <- old_dat$name[which(old_dat$snumber == raw_dat$snumber[i])[1]]
  } else {
    raw_dat$name[i] <- "NO_NAME"
  }
}

# Get rid of the double flood values for 1 year (we take the max by setting the min to NA)
for (i in 1:(length(raw_dat$year) - 1)) {
  if (!is.na(raw_dat$year[i])) {
    if (raw_dat$year[i] == raw_dat$year[i + 1] && any(is.na(c(raw_dat$value_q_daily[i], raw_dat$value_q_daily[i + 1]))) == FALSE ) {
      min_index <- which(c(raw_dat$value_q_daily[i], raw_dat$value_q_daily[i + 1]) == min(c(raw_dat$value_q_daily[i], raw_dat$value_q_daily[i + 1])))
      raw_dat$value_q_daily[i - 1 + min_index] <- NA
    }
  }
}

# Let's first get rid of any data line that has an NA flow and create the dat list
dat <- list()
keep <- which(!is.na(raw_dat$value_q_daily))

dat$regine_area <- raw_dat$regine_area[keep]
dat$main_nr <- raw_dat$main_nr[keep]
dat$snumber <- raw_dat$snumber[keep]
dat$name <- raw_dat$name[keep]

dat$flom_DOGN <- raw_dat$value_q_daily[keep]
dat$date <- raw_dat$dt_flood_daily_date[keep]
dat$year <- raw_dat$year[keep]

dat$fraction_rain <- raw_dat$fraction_rain[keep]

# Making list of snumber and getting rid of stations with not enough data
station.nb.vect.init <- na.omit(unique(dat$snumber))

min_years_data <- 30  # set the minimum number of years to accept a station
max_years_ss <-90 # The maximum record length for subsampling
station.nb.vect <- c()
length_rec <- c()

for (i in seq(along = station.nb.vect.init)) {

  length_record <- length(which(dat$snumber == station.nb.vect.init[i]))
  if (length_record >= min_years_data) {
    station.nb.vect <- c(station.nb.vect, station.nb.vect.init[i])
    length_rec <- c(length_rec, length_record)
  }

}

dim.station <- length(station.nb.vect)
#station.nb.vect <- station.nb.vect[1:dim.station]


# Read station coordinates
utminfo <- read.table("../rawdata/Coordinates_for_R.txt",
                      sep = "\t", header = T)  # CHECK DIR
dat <- save_coordinates(dat, station.nb.vect)  # this function is in LOAD.R

# Read catchments area
catchment.prop <- read.csv("../rawdata/Hydra_FeltparTabell.csv", sep=";")  # CHECK DIR
catchment.prop$COMPOUND_K <- catchment.prop$COMPOUND_K / 1000  # to get the same numbers as in the main data file
catchment.size <- rep("NA", length(station.nb.vect))
catchment.min.height <- rep("NA", length(station.nb.vect))
catchment.max.height <- rep("NA", length(station.nb.vect))

indexes.in.catchmentprop  <- rep("NA", length(station.nb.vect))

for (i in seq(along = station.nb.vect)) {
if (is.numeric(which(catchment.prop$COMPOUND_K == station.nb.vect[i]))) {
  if (as.numeric(length(which(catchment.prop$COMPOUND_K == station.nb.vect[i]))) > 0) {
    indexes.in.catchmentprop[i] <- which(catchment.prop$COMPOUND_K == station.nb.vect[i])[1]
    # [1] because some snumbers occur twice in the Hydra_FeltparTabell.csv
  }
}
}

indexes.in.catchmentprop <- as.numeric(indexes.in.catchmentprop)

catchment.size <- as.numeric(catchment.prop$AREAL_UTM3[as.numeric(indexes.in.catchmentprop)])
catchment.min.height <- catchment.prop$HEIGHT_MIN[as.numeric(indexes.in.catchmentprop)]
catchment.max.height <- catchment.prop$HEIGHT_MAX[as.numeric(indexes.in.catchmentprop)]

################################################
# detect number of cores for parallel computing
ncores<-detectCores(all.tests = FALSE, logical = TRUE)

dim.distr <- 5
dim.method <- 4
dim.param <- 3
dim.length_total_record <- 150
# dim.length_rec <- 30   # This means that we will do subsampling for i1=30,i2=35,.....i23=140 years of record.
dim.random_runs <- 50
dim.characters <- 64
sampling_years <- seq(min_years_data, 5 * trunc(max_years_ss / 5), 5)  # Subsampling from min_years_data until 90 years of data, increments of 5 years
dim.max_subsample <- max(sampling_years)
dim.length_rec <- length(sampling_years)+1   # The last index is for storing the full record.
distr.name <- c("gumbel", "gamma", "gev", "gl", "pearson")
method.name <- c("mle", "Lmom", "mom", "bayes")

}


#' create_nc
#' @description Function to create an empty nc database for FlomKart
#' @return
#' @export
#'
#' @examples
create_nc <- function() {

## To run if creating the nc file from scratch  --------------

# Creation of the CDF dataset and definition of the dimensions
nc <- create.nc("../output/flood_database.nc")  # CHECK DIR
att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", "Flood frequency analysis results")
att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR", paste("Created on", base::date()))

station.name <- rep("NA", dim.station)
station.utmN <- rep(-9999, dim.station)
station.utmE <- rep(-9999, dim.station)
station.long <- rep(-9999, dim.station)
station.lat <- rep(-9999, dim.station)
station.lat <- rep(-9999, dim.station)
station.length_rec <- length_rec


Q <- array(NA,dim=c(dim.station, dim.length_total_record))
years <- array(NA,dim=c(dim.station, dim.length_total_record))
dates <- array(NA,dim=c(dim.station, dim.length_total_record))

# Loop over the stations to fill variables: Q, station.name
for (i in seq(along = station.nb.vect)) {

  indexes <- which(dat$snumber == station.nb.vect[i])

  # Matrices
  Q[i, 1:length_rec[i]] <- dat$flom_DOGN[indexes]
  years[i, 1:length_rec[i]] <- dat$year[indexes]
  dates[i, 1:length_rec[i]] <- dat$date[indexes]

  # Vectors
  station.name[i] <-  as.character(unique(dat$name[indexes]))
  station.utmN[i] <- unique(dat$utmN[indexes])
  station.utmE[i] <-  unique(dat$utmE[indexes])
  station.long[i] <-  unique(dat$long[indexes])
  station.lat[i] <-  unique(dat$lat[indexes])
}

random_indexes <- array(NA, dim = c(dim.station, dim.random_runs, length(sampling_years), dim.max_subsample))

dim.def.nc(nc, "station", dim.station)
dim.def.nc(nc, "distr", dim.distr)
dim.def.nc(nc, "method", dim.method)
dim.def.nc(nc, "param", dim.param)
dim.def.nc(nc, "length.rec", dim.length_rec)
dim.def.nc(nc, "length_total_record", dim.length_total_record)
dim.def.nc(nc, "random_runs", dim.random_runs)
dim.def.nc(nc, "max_string_length", dim.characters)
dim.def.nc(nc, "subsampling", length(sampling_years))
dim.def.nc(nc, "max_subsample", dim.max_subsample)
dim.def.nc(nc, "scalars", 1)
sync.nc(nc)

# Definition of the variables and which dimension they are linked to
# Input data
var.def.nc(nc, varname = "Q", vartype = "NC_FLOAT", dimensions = c("station", "length_total_record"))
att.put.nc(nc, "Q", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "Q", "short_name", "NC_CHAR", "Flood record")
att.put.nc(nc, "Q", "long_name", "NC_CHAR", "Yearly flood records (averaged maximum daily flow in m3/s) for Norway")
var.put.nc(nc, "Q", Q)
sync.nc(nc)

var.def.nc(nc, varname = "years", vartype = "NC_INT", dimensions = c("station", "length_total_record"))
att.put.nc(nc, "years", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "years", "short_name", "NC_CHAR", "Years of record")
att.put.nc(nc, "years", "long_name", "NC_CHAR", "Years of record")
var.put.nc(nc, "years", years)
sync.nc(nc)

var.def.nc(nc, varname = "dates", vartype = "NC_CHAR", dimensions = c("max_string_length", "station", "length_total_record"))
att.put.nc(nc, "dates", "missing_value", "NC_CHAR", "NA")
att.put.nc(nc, "dates", "short_name", "NC_CHAR", "Dates of record")
att.put.nc(nc, "dates", "long_name", "NC_CHAR", "Dates of record")
datesc<-dates
datesc[is.na(datesc)]<-"NA"
var.put.nc(nc, "dates", datesc)
sync.nc(nc)

var.def.nc(nc,  varname = "distr.name", vartype = "NC_CHAR", dimensions = c("max_string_length", "distr"))
att.put.nc(nc, "distr.name", "missing_value", "NC_CHAR", "NA")
att.put.nc(nc, "distr.name", "long_name", "NC_CHAR", "Distributions used: gumbel, gamma, gev, gl, pearson")
var.put.nc(nc, "distr.name", distr.name)
sync.nc(nc)

var.def.nc(nc,  varname = "method.name", vartype = "NC_CHAR", dimensions = c("max_string_length", "method"))
att.put.nc(nc, "method.name", "missing_value", "NC_CHAR", "NA")
att.put.nc(nc, "method.name", "long_name", "NC_CHAR", "Estimation methods used: mle, Lmom, mom, bayes")
var.put.nc(nc, "method.name", method.name)
sync.nc(nc)

var.def.nc(nc,  varname = "station.name", vartype = "NC_CHAR", dimensions = c("max_string_length", "station"))
att.put.nc(nc, "station.name", "missing_value", "NC_CHAR", "NA")
var.put.nc(nc, "station.name", station.name)
sync.nc(nc)

var.def.nc(nc,  varname = "station.utmN", vartype = "NC_INT", dimensions = "station")
att.put.nc(nc, "station.utmN", "missing_value", "NC_INT", -9999)
var.put.nc(nc, "station.utmN", station.utmN)
sync.nc(nc)

var.def.nc(nc,  varname = "station.utmE", vartype = "NC_INT", dimensions = "station")
att.put.nc(nc, "station.utmE", "missing_value", "NC_INT", -9999)
var.put.nc(nc, "station.utmE", station.utmE)
sync.nc(nc)

var.def.nc(nc,  varname = "station.long", vartype = "NC_FLOAT", dimensions = "station")
att.put.nc(nc, "station.long", "missing_value", "NC_FLOAT", -9999)
var.put.nc(nc, "station.long", station.long)

var.def.nc(nc,  varname = "station.lat", vartype = "NC_FLOAT", dimensions = "station")
att.put.nc(nc, "station.lat", "missing_value", "NC_FLOAT", -9999)
var.put.nc(nc, "station.lat", station.lat)

var.def.nc(nc,  varname = "station.number", vartype = "NC_INT", dimensions = "station")
att.put.nc(nc, "station.number", "missing_value", "NC_INT", -9999)
var.put.nc(nc, "station.number", station.nb.vect)

var.def.nc(nc, varname = "catchment.size", vartype = "NC_FLOAT", dimensions = "station")
att.put.nc(nc, "catchment.size", "missing_value", "NC_FLOAT", "NA")
var.put.nc(nc, "catchment.size", catchment.size)

var.def.nc(nc, varname = "catchment.min.height", vartype = "NC_FLOAT", dimensions = "station")
att.put.nc(nc, "catchment.min.height", "missing_value", "NC_CHAR", "NA")
var.put.nc(nc, "catchment.min.height", catchment.min.height)

var.def.nc(nc, varname = "catchment.max.height", vartype = "NC_FLOAT", dimensions = "station")
att.put.nc(nc, "catchment.max.height", "missing_value", "NC_CHAR", "NA")
var.put.nc(nc, "catchment.max.height", catchment.max.height)
sync.nc(nc)

# Saving parameters that will be used in creating the fitting data
# First SCALARS
var.def.nc(nc,  varname = "min_years_data", vartype = "NC_INT", dimensions = NA)
att.put.nc(nc, "min_years_data", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "min_years_data", "short_name", "NC_CHAR", "The minimum number of years the random subsamples have")
var.put.nc(nc, "min_years_data", min_years_data)

var.def.nc(nc,  varname = "dim.random_runs", vartype = "NC_INT", dimensions = NA)
att.put.nc(nc, "dim.random_runs", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "dim.random_runs", "short_name", "NC_CHAR", "The number of random subsamples for every subsampling length")
var.put.nc(nc, "dim.random_runs", dim.random_runs)

var.def.nc(nc,  varname = "dim.length_rec", vartype = "NC_INT", dimensions = NA)
att.put.nc(nc, "dim.length_rec", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "dim.length_rec", "long_name", "NC_CHAR", "The length of the dimension where subsampling results will be stored.")
att.put.nc(nc, "dim.length_rec", "short_name", "NC_CHAR", "This is 30: 23 subsamples up to 140 years of data + margin + i=30 for the full record")
var.put.nc(nc, "dim.length_rec", dim.length_rec)
sync.nc(nc)

# THEN VECTORS AND CO.
var.def.nc(nc,  varname = "sampling_years", vartype = "NC_INT", dimensions = "subsampling")
att.put.nc(nc, "sampling_years", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "sampling_years", "short_name", "NC_CHAR", "The various length of record used for subsampling")
var.put.nc(nc, "sampling_years", sampling_years)

var.def.nc(nc,  varname = "random_indexes", vartype = "NC_INT", dimensions = c("station", "random_runs", "subsampling", "max_subsample"))
att.put.nc(nc, "random_indexes", "missing_value", "NC_INT", -9999)
att.put.nc(nc, "random_indexes", "short_name", "NC_CHAR", "The various indexes used for subsampling")
var.put.nc(nc, "random_indexes", random_indexes)

# Model results (parameter estimates and standard error)
var.def.nc(nc, varname = "param.estimate", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "param", "length.rec", "random_runs"))
att.put.nc(nc, "param.estimate", "missing_value", "NC_FLOAT", -9999)
param.estimate <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.param, dim.length_rec, dim.random_runs))
var.put.nc(nc, "param.estimate", param.estimate)
rm(param.estimate)

var.def.nc(nc, varname = "param.se", vartype = "NC_FLOAT",
           dimensions = c("station", "distr", "method", "param", "length.rec", "random_runs"))
att.put.nc(nc, "param.se", "missing_value", "NC_FLOAT", -9999)
param.se <- array(NA,dim=c(dim.station, dim.distr, dim.method, dim.param, dim.length_rec, dim.random_runs))
var.put.nc(nc, "param.se", param.se)
rm(param.se)

sync.nc(nc)  # Save what we've done so far

# To run if the input data part of the NC file  is not yet defined----------------------------------------------

nc <- open.nc("../output/flood_database.nc", write = TRUE)  # Put FALSE for read-only  # CHECK DIR

var.put.nc(nc, "Q", Q)
var.put.nc(nc, "station.name", station.name)
var.put.nc(nc, "station.utmN", station.utmN)
var.put.nc(nc, "station.utmE", station.utmE)
var.put.nc(nc, "station.long", station.long)
var.put.nc(nc, "station.lat", station.lat)

sync.nc(nc)

}


#' fillup_nc
#' @description Function to fill upthe NetCDF database
#' @return
#' @export
#'
#' @examples
fillup_nc <- function() {

## To run if updating or creating the nc file from scratch  --------------

nc <- open.nc("../output/flood_database.nc", write = TRUE)  # Put FALSE for read-only  # CHECK DIR
Q <- var.get.nc(nc, "Q")

# Dumping all console output into errorlog.txt
sink("../output/errorlog.txt")  # CHECK DIR

if(mysystem=="Windows")cl<-makeCluster(ncores-1) # leave one core free
if(mysystem=="Linux") registerDoMC(ncores-15) # for starting the parallel calculations
if(mysystem=="Windows") registerDoSNOW(cl) # for starting the parallel calculations
# out.temp<-foreach(st=1:10,.packages='fitdistrib') %dopar% {
# foreach starts the paralell loop. In order to gain computing speed when using non-Bayesian methods- the parallel loop should be here
  out.temp<-foreach(st=1:length(station.nb.vect),.packages='fitdistrib') %dopar% {

# for (st in seq(along = station.nb.vect)) {
# for (st in 1:5) {
  print(st)

# Need to store results from each parallel chain in temporary arrays
  par.param.estimate <- array(NA,dim=c(dim.distr,dim.method,dim.param, dim.length_rec,dim.random_runs))
  par.param.se <- array(NA,dim=c(dim.distr,dim.method,dim.param, dim.length_rec,dim.random_runs))


  temp.Q <- as.vector(na.omit(Q[st, ]))

# we might resample with relplacement records that are longer than the original one. No limitation on the original recird length
  temp.random_indexes <- array(NA, dim = c(dim.random_runs, length(sampling_years), dim.max_subsample))
  if (length(temp.Q) <  min_years_data) {
    print("This station number has NULL or not enough data")
    print(station.nb.vect[st])
    next(st)
  } else {

    # Computed random indexes before the for loops on distr and methods
    # We now resample with replacement records of length min
    for (rs in 1:dim.random_runs) {  # sampling.int is by default without replacement
      for (j in 1 : length(sampling_years)){
        # now we sample with replacement. We could also sample with replacement for the full length of record
        # for each station but this requires some playing with indices
        temp.random_indexes[rs, j, 1:sampling_years[j]] <- sample.int(length(temp.Q), sampling_years[j], replace = TRUE)
      } # for j
    } # for rs



    for (d in 1:5) {
      distr <- distr.name[d]
      print(distr)
      for (m in 1:4) {  # k  # no bayes
        # m=4  # k
        method <- method.name[m]
        print(method)
        print("Computing full dataset")

        # Allocation the appropriate function for d and m
        FUN <- match.fun(paste(as.character(distr), "_", as.character(method), collapse = "", sep = ""))
        param <- FUN(temp.Q)


        if (length(na.omit(param$estimate)) == length(param$estimate)) {

          par.param.estimate[d,m,1:2,dim.length_rec,1] <- param$estimate[1:2]
          par.param.se[d,m,1:2,dim.length_rec,1] <- param$se[1:2]

          if (length(param$estimate == 3)) {
            par.param.estimate[d,m,3,dim.length_rec,1] <- param$estimate[3]
            par.param.se[d,m,3,dim.length_rec,1] <- param$se[3]
          } # if
        } # if


        if (m < 4) {  # Random sampling of Bayes takes a very long time  # k
           for (j in 1 : length(sampling_years)){

            for (rs in 1:dim.random_runs) {
              sample_q <- temp.Q[temp.random_indexes[rs, j, 1:sampling_years[j]]]
              param.sample <- FUN(sample_q)

              if (length(na.omit(param.sample$estimate)) == length(param.sample$estimate)) {

                par.param.estimate[d,m,1:2,j,rs] <- param.sample$estimate[1:2]
                par.param.se[d,m,1:2,j,rs] <- param.sample$se[1:2]

                if (length(param.sample$estimate ==3)) {
                  par.param.estimate[d,m,3,j,rs] <- param.sample$estimate[3]
                  par.param.se[d,m,3,j,rs] <- param.sample$se[3]
                } # if

              } # if

            } # for rs

          } # for j

        }  # if m < 4
      } # for m
    } # for d
  } # else
  out<-list(par.param.estimate,par.param.se,temp.random_indexes,temp.Q)
} # dopar
# Write the output from all paralell calculations to the NetCDF file.
# Maybe later this might be improved by letting each paralell loop write to the same NetCDF file.
for(st in 1 : length(station.nb.vect)){
          var.put.nc(nc, "param.estimate", data = out.temp[[st]][[1]], start=c(st, 1, 1, 1, 1, 1),
                     count=c(1, dim.distr, dim.method, dim.param,dim.length_rec, dim.random_runs))
          var.put.nc(nc, "param.se", data = out.temp[[st]][[2]], start=c(st, 1, 1, 1, 1, 1),
                     count=c(1, dim.distr, dim.method, dim.param,dim.length_rec, dim.random_runs))
          var.put.nc(nc, "random_indexes", data = out.temp[[st]][[3]], start=c(st, 1, 1, 1),
                       count=c(1, dim.random_runs, length(sampling_years), dim.max_subsample))


}
rm(out.temp)
stopCluster(cl) # stop the cluster
sink()  # Stop printing console onto file
sync.nc(nc)
close.nc(nc)

}  # end of fillup_nc function

# Names of stations that appear twice. This is not a problem, they are different stations
# because they have a different number
# [1] "Sandvatn"
# [1] "Kjeldstad i Garbergelva"
# [1] "Storvatn"
# [1] "SOndre BjOllaavatn"
# [1] "Russaanes"
# [1] "Stabburselv"
# [1] "Helligskogen"
# [1] "Lillefossen"