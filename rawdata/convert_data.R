dat <- read.table("rawdata/2016-10-14-amsvalues.txt", sep = ";", header = TRUE)
save(dat, file = "data/preprocessed_flood_data.RData")

library(lubridate)
date_vect <- ymd(dat$daily_ams_dates)

postproc_dat <- data.frame(regine = dat$regine, main = dat$main,
                           year = year(date_vect), month = month(date_vect), day = day(date_vect),
                           flow = dat$daily_ams)

save(postproc_dat, file = "data/flood_data.RData")


###########################


# source('LOAD.r')      # Required to load Q into database and save station coordinates


#################################################################################################################
#################################################################################################################

# CLARIFY THIS raw_dat thing!!

# Old data for station names only
old_dat <- read.csv("rawdata/AMS_table_old.csv", sep=";", as.is=TRUE)  # CHECK DIR

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

#################################################################################################################
#################################################################################################################


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

## !!! Station number vector to put into the metadata file
flood_metadata <- list()
flood_metadata$station.nb.vect <- station.nb.vect
flood_metadata$min_years_data <- 30
save(flood_metadata, file = "data/flood_metadata.RData")


# Read station coordinates
utminfo <- read.table("rawdata/Coordinates_for_R.txt",
                      sep = "\t", header = T)  # CHECK DIR
dat <- save_coordinates(dat, station.nb.vect)  # this function is in LOAD.R

# Read catchments area
catchment.prop <- read.csv("rawdata/Hydra_FeltparTabell.csv", sep=";")  # CHECK DIR
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

