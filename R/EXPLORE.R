# EXPLORE.R
# EXPERIMENTS WIH USING R TO PLOT MAPS. BASED ON LENA'S SCRIPTS
# Functions related to the initial exploration of the dataset

plot_norway <- function() {
# Function to plot the shape file of Norway.
# This comes from Lena
  
  #load KlausKart in R
  newproj <- "+proj=utm +zone=33 +north +units=m"
  ftext <- "//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/GIS/norge.shp"
  shape <- readShapeSpatial(ftext, proj4string = CRS(newproj), repair = TRUE, force_ring = T, verbose = TRUE)
  
  #ftext<-"../data/33_N2000_Norge"
  ogrInfo(dsn = ftext, layer = "norge")
  temp <- readOGR(ftext, layer = "norge")
  x_min <- min(summary(temp)$bbox[1])
  x_max <- max(summary(temp)$bbox[3])
  y_min <- min(summary(temp)$bbox[2])
  y_max <- max(summary(temp)$bbox[4])
  
  #plot norway map
  plot(temp, border = "grey", axes = FALSE, xlim = c(x_min, x_max), ylim = c(y_min, y_max),
       asp = 1, cex.main = 1.5) #col="",
  
}

plot_all_stations <- function(dat, station.nb.vect, distr = "distr", method = "method") {
# Loop to plot all stations with a dot.
# Calls the function "plot_station" to calculate the size and color of the dot based on a some metric
  
  utminfo <- read.table("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/Coordinater for kart_for_R.txt",  sep = "\t", header = T)
  
  windows()
  titel = "Map of Norway"
  plot_norway()
  
  for (i in seq(along = station.nb.vect)) {
    sdat <- sdat_load(dat, station.nb.vect[i])
 
  if (is.null(sdat) == FALSE)  {
    param <- fit_distr(sdat$normQ, distr = distr, method = method)   
    if (is.na(param) == FALSE) {
    plot_station(sdat, utminfo, param$estimate[1])
    par(new = TRUE)  
  }
  }  
}
}

plot_station <- function(sdat, utminfo, dot.size) {
# Plots each station at its particular location with dot_size a input variable (i.e. goodness of fit) 
# Should probably be replaced completely by plot_sation2
  
  sdat$stn.regine <- trunc(sdat$snumber / 100000)
  sdat$stn.main <- as.numeric(substr(sdat$snumber / 1000, nchar(sdat$snumber / 1000) - 2, 
                                   nchar(sdat$snumber / 1000)))
  sdat$stn.nve <- paste(sdat$stn.regine, ".", sdat$stn.main,sep="")

  # match coordinates into your file
  sdat$utmN <- NA
  sdat$utmE <- NA

  for (u in 1:nrow(sdat)) {
    loc.tmp <- which(utminfo$regine_area == sdat$stn.regine[u] & utminfo$main_no == sdat$stn.main[u] ) 
    if(length(loc.tmp) > 0){
      #sdat$regine[u] <- unique_allstations$regine_enhet[loc.tmp]
    if (length(loc.tmp) > 1){
      loc.tmp <- loc.tmp[1]
    }
    sdat$utmN[u] <-utminfo$utm_north_z33[loc.tmp]
    sdat$utmE[u] <-utminfo$utm_east_z33[loc.tmp]
    }
  }
  points(x = sdat$utmE[1], y = sdat$utmN[1], pch = 19, cex= 0.25 * dot.size, col = "blue") 
}

plot_station2 <- function(sdat, utminfo, dot.size, dot.color) {
# Plots each station at its particular location with dot_size and dot_color as input variables (i.e. goodness of fit)  
  
  sdat$stn.regine <- trunc(sdat$snumber / 100000)
  sdat$stn.main <- as.numeric(substr(sdat$snumber / 1000, nchar(sdat$snumber / 1000) - 2, 
                                     nchar(sdat$snumber / 1000)))
  sdat$stn.nve <- paste(sdat$stn.regine, ".", sdat$stn.main, sep = "")
  
  # match coordinates into your file
  sdat$utmN <- NA
  sdat$utmE <- NA

  loc.tmp <- which(utminfo$regine_area == sdat$stn.regine[1] & utminfo$main_no == sdat$stn.main[1] ) 
  if(length(loc.tmp) > 0){
    loc.tmp <- loc.tmp[1]
    sdat$utmN <-utminfo$utm_north_z33[loc.tmp]
    sdat$utmE <-utminfo$utm_east_z33[loc.tmp]
  }
  
  # plot (labeled) stations (if you want to plot only certain ones, then: 
  # loc.tmp <- which(sdat$yourvariable > threshold),   points(x=sdat$utmE[loc.tmp],....  )
  
  # tiff("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Plots/map.tif",  
  #      width = 14, height = 15,  pointsize = 1/300, units = 'in', res = 250)
  
  #plot your stations
  points(x = sdat$utmE[1], y = sdat$utmN[1], pch = 19, cex= 1.5 * dot.size, 
         col = rainbowPalette(5)[as.numeric(dot.color)])  
  
  #if you want label them
  # text(sdat$utmE[1], sdat$utmN[1],sdat$stn.nve[1], cex=0.55)
  
}

plot_gof <- function(dat, station.nb.vect, method = "MLE", input = input) {
# This is like plot_all_stations but for the variable to plot for each station is given as input
# The structure could be optimized with plot_all_stations
  
  # This is the return of run.all_KS x <- list(snumber = station.nb.vect, name = st_name, ks = ks, distrib = distrib)
  # input <- run.all_KS(dat, station.nb.vect, method = "MLE")
  
  utminfo <- read.table("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/Coordinater for kart_for_R.txt",  sep = "\t", header = T)
  
  windows()
  titel = "Map of Norway"
  plot_norway()
  
  for (i in seq(along = station.nb.vect)) {
    
    if (is.na(input$distrib[i]) == TRUE) next(i)
    
    if (input$distrib[i] == "GAMMA") input$distrib[i] <- as.numeric(1)
    if (input$distrib[i] == "GUMBEL") input$distrib[i] <- as.numeric(2)
    if (input$distrib[i] == "GEV") input$distrib[i] <- as.numeric(3)
    if (input$distrib[i] == "GL") input$distrib[i] <- as.numeric(4)
    if (input$distrib[i] == "PEARSON") input$distrib[i] <- as.numeric(5)
    
    if (!is.null(input$gof[i]) > 0 && !is.na(input$gof[i]))  {
      sdat <- sdat_load(dat, station.nb.vect[i])
      plot_station2(sdat, utminfo, input$gof[i], input$distrib[i])
      par(new = TRUE)  
    }  
  }
}

plot_sample <- function(sdat, station.nb.vect, radius) {
# This function plots all the stations that are located within a particular radius
# of a specified station
  utminfo <- read.table("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/Coordinater for kart_for_R.txt",  sep = "\t", header = T)
  
  sdat$stn.regine <- trunc(sdat$snumber / 100000)
  sdat$stn.main <- as.numeric(substr(sdat$snumber / 1000, nchar(sdat$snumber / 1000) - 2, 
                                     nchar(sdat$snumber / 1000)))
  sdat$stn.nve <- paste(sdat$stn.regine, ".", sdat$stn.main,sep="")
  
  for (u in 1:nrow(sdat)) {
    loc.tmp <- which(utminfo$regine_area == sdat$stn.regine[u] & utminfo$main_no == sdat$stn.main[u] ) 
    if(length(loc.tmp) > 0){
      #sdat$regine[u] <- unique_allstations$regine_enhet[loc.tmp]
      if (length(loc.tmp) > 1){
        loc.tmp <- loc.tmp[1]
      }
      sdat$utmN[u] <-utminfo$utm_north_z33[loc.tmp]
      sdat$utmE[u] <-utminfo$utm_east_z33[loc.tmp]
    }
  }
  points(x = sdat$utmE[1], y = sdat$utmN[1], pch = 19, cex= 0.25 * dot.size, col = "blue") 
  
  
  distances <- sqrt((utminfo$utm_east_z33[2]-utminfo$utm_east_z33[1])^2+(utminfo$utm_north_z33[2]-utminfo$utm_north_z33[1])^2)
  
}


