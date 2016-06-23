# LOAD.R
# This file gather all the functions related to loading and sampling data

sdat_load <- function(dat, station.nb) {
# Load flood dat for a specific station
  
  # Extracting one station only, and when aarsmax_dogn = TRUE
  x <- which(dat$snumber == station.nb)
  if (length(x) == 0) {
    return(sdat <- NULL)
    stop("This station number does not exist") 
  } else {
    y <- intersect(x, which(dat$aarsmax_dogn == TRUE)) 
  }
  if (length(y) == 0) {
  print(paste("Warning: The station", dat$name[x[1]], "number", dat$snumber[x[1]], "exists but doesn't have data", sep = " ")) 
  return(sdat <- NULL)
  } else {  
      # The list sdat is the selection of the 'y' indexes in dat
    sdat <- dat[y, ]
    sdat$flom_DOGN <- na.omit(sdat$flom_DOGN)  # Need to check if it's really a good idea, otherwise put the na.omit for everyone
    sdat$logQ <- log(sdat$flom_DOGN)
    # sdat$normQ <- (sdat$flom_DOGN-mean(sdat$flom_DOGN)) / sd(sdat$flom_DOGN)
    sdat$normQ <- (sdat$flom_DOGN-min(sdat$flom_DOGN)) / (max(sdat$flom_DOGN)-min(sdat$flom_DOGN))
  }
  invisible(sdat)
}

save_coordinates <- function(dat, station.nb.vect) {
  # This function saves a new csv with the utm coordinates of each station in station.nb.vect
  # NEED TO THINK ABOUT THE OVERALL EFFICIENCY OF THOSE GIS FUNCTIONS
  
  
  utminfo <- read.table("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/Coordinater for kart_for_R.txt",  sep = "\t", header = T)
  

  dat$stn.nve <- paste(dat$regine_area, ".", dat$main_nr,sep="")
  
  for (u in seq(along = dat$stn.nve)) {
    loc.tmp <- which(utminfo$regine_area == dat$regine_area[u] & utminfo$main_no == dat$main_nr[u] ) 
    if(length(loc.tmp) == 0){
      dat$utmN[u] <- NA
      dat$utmE[u] <- NA
      next(i)
    } else {
      loc.tmp <- loc.tmp[1]
      dat$utmN[u] <- utminfo$utm_north_z33[loc.tmp]
      dat$utmE[u] <- utminfo$utm_east_z33[loc.tmp]
      dat$longitude[u] <- utminfo$longitude[loc.tmp]
      dat$latitude[u] <- utminfo$latitude[loc.tmp]
    }
  }
  # You can Write the new data table to a csv file
  write.table(dat, file = "georef_data.csv", sep = ";", 
              col.names = NA, qmethod = "double")
  
  invisible(dat)
}

##############################################################

sdat_plot <- function(dat) {
# Function to plot data histogram, pdf, ecdf and qqplot

  windows()
  par(mfrow = c(3, 1))
  
  xmax <- 1.1 * max(dat)
  xmin <- 0.9 * min(dat)
  hist(dat, xlab = "Flood discharge (m3/s)", freq = FALSE, 
  breaks = seq(xmin, xmax, xmax / 15), col="gray", main = NULL)  
  par(new = TRUE)
  lines(density(dat))
  plot(ecdf(dat), main = "Cumulate distribution ")
  qqnorm(dat, main = "QQ Plot ")  # drawing the QQplot
  abline(0,1)  # drawing a 45-degree reference line
}


load_samples <- function(dat.full, k) {
# Function that creates smaller samples from the original data and return a matrix of the possible samples  

  L <- length(dat.fulll) - 1  # Length of initial dataset
  stab <- list(max = c(), mean = c())
  
  for (j in k:L) {  # i < j so i:j is the NA window we create
    i <- j-k+1
    dat.sample <- dat.fulll
    dat.sample[i:j] <- NA  # Creation of the NA window and na.omit to create a reduced dataset
    dat.sample <- as.vector(na.omit(dat.sample))
  }
}

# Function that returns only 1 sample with k omited values, according to a specific method (random, sorted min and max)
load_1sample <- function(dat.fulll, k, method = "method") {
  
  if (method == "min") {
    dat.new <- sort(dat.fulll)
    dat.new <- dat.new[k:length(dat.new)]
  }
  else if (method == "max") {
    dat.new <- sort(dat.fulll)
    dat.new <- dat.new[1:(length(dat.new) - k)]
  }
}