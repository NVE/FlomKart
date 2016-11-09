###############################


# Distribution of the length data at stations
length_dat <- vector(length = length(station.nb.vect))

for (i in seq(along = station.nb.vect)) {
  sdat <- sdat_load(dat, station.nb.vect[i])
  length_dat[i] <- length(sdat$flom_DOGN)

}
hist(length_dat, breaks = 20 , xlab = "Lenth of flood record", ylab = "Number of stations", main= "", cex.lab = 1.5)
# par(new = TRUE)
# plot(density(length_dat))


#################################
# LMOM PLOT
l1.dat <- c()
l2.dat <- c()
t1.dat <- c()
t2.dat <- c()
lmom.data <- data.frame(l_1 = c(), l_2 = c(), t_3 = c(), t_4 = c())
for (i in seq(along = station.nb.vect)) {
  sdat <- sdat_load(dat, station.nb.vect[i])
  if (!is.null(sdat$flom_DOGN) && length(samlmu(sdat$flom_DOGN)) == 4) {
    l1.dat[i] <- samlmu(sdat$flom_DOGN)[1]
    l2.dat[i] <- samlmu(sdat$flom_DOGN)[2]
    t1.dat[i] <- samlmu(sdat$flom_DOGN)[3]
    t2.dat[i] <- samlmu(sdat$flom_DOGN)[4]
  }
}
lmrd(lmom.data)


windows()
lmom.data <- vector()
for (i in seq(along = station.nb.vect)) {
  sdat <- sdat_load(dat, station.nb.vect[i])
  if (!is.null(sdat$flom_DOGN) && length(samlmu(sdat$flom_DOGN)) == 4) {
    lmom.data <- samlmu(sdat$flom_DOGN)
    lmrd(lmom.data, legend.lmrd = FALSE)
    # par(new = TRUE)
  }
}

sdat <- sdat_load(dat, station.nb.vect[2])

lmom.data <- samlmu(sdat$flom_DOGN)
lmrd(lmom.data)

################### TEST WITH LMOMCO PKG
library(lmomco)
t3 <- vector(mode = "numeric"); t4 <- t3;
for(i in 1:length(station.nb.vect)) {
  sdat <- sdat_load(dat, station.nb.vect[i])
  if (length(sdat$flom_DOGN) > 10) {
    lmom.data <- lmoms(sdat$flom_DOGN)
    t3[i] <- lmom.data$ratios[3]; t4[i] <- lmom.data$ratios[4]
  }
}

# Finally, plot the diagram with a legend at a specified location,
# and "zoom" into the diagram by setting the axis limits.
plotlmrdia(lmrdia(), autolegend=TRUE, xleg=0.1, yleg=.41,
           xlim=c(-.1,.5), ylim=c(-.1,.4), nopoints=TRUE)
points(t3,t4)
points(mean(t3),mean(t4),pch=16,cex=3)
lines(c(T3,T3),c(-1,1),col=8, lty=2)
lines(c(-1,1),c(T4,T4),col=8, lty=2)







#############################
## Test of Kolbjorns code

library(RNetCDF)
library(lmom)

nc <- open.nc("flood_database.nc", write = FALSE)  # Put FALSE for read-only
myq<-var.get.nc(nc,1)
rlength<-apply(!is.na(myq),1,sum)
qlmom<-t(apply(var.get.nc(nc,1),1,samlmu))

lmrd(t(apply(var.get.nc(nc,1),1,samlmu)),pch=16,cex=0.5,xlim=c(0,0.7),ylim=c(0,0.7))
malmom<-get_ma_lmom(qlmom,rlength,40)
points(malmom$lskew,malmom$lcurt,col=2,pch=16,cex=0.6)



get_ma_lmom<-function(elmom,tlength,wsize){

  tlength<-tlength[!is.na(elmom[,4])]

  qskew<-elmom[!is.na(elmom[,4]),3]
  qcurt<-qlmom[!is.na(elmom[,4]),4]
  qskew[rank(qskew)]

  tlentgh<-tlength[order(qskew)]
  qcurt<-qcurt[order(qskew)]
  qskew<-qskew[order(qskew)]

  ns<-length(qskew)
  np2<-ns-40+1
  newcurt<-rep(NA,np2)
  newskew<-rep(NA,np2)
  for(i in 1:np2){
    i2=i+wsize-1
    newcurt[i]=sum(qcurt[i:i2]*rlength[i:i2])/sum(rlength[i:i2])
    newskew[i]=sum(qskew[i:i2]*rlength[i:i2])/sum(rlength[i:i2])
  }
  out<-list()
  out$lskew<-newskew
  out$lcurt<-newcurt
  out
}


source('R/global.R')


# station_group_indexes(gof, distr, method, minmax)
best.ks <- data.frame(method = method.name, gumbel = rep(NA,4), gamma = rep(NA,4), gev = rep(NA,4), gl = rep(NA,4),
                      pearson = rep(NA,4))

worst.ks <- data.frame(method = method.name, gumbel = rep(NA,4), gamma = rep(NA,4), gev = rep(NA,4), gl = rep(NA,4),
                       pearson = rep(NA,4))

for (d in 1:5) {
  for (m in 1:4) {
    best.ks[m, d+1] <- length(station_group_indexes("KS", distr.name[d], method.name[m], "min"))  # length() because we want the number of stations
    # that qualify as best for this method and distr
    worst.ks[m, d+1] <- length(station_group_indexes("KS", distr.name[d], method.name[m], "max"))
  }

}


#### Attempt to plot the average evolution of a GOF as f(record length) with all stations of more than 60 years of data

# d <- 3  # GEV
# m <- 2  # Lmom
#
# d <- 4  # GL
# m <- 1  # mle
#
# d <- 4  # GL
# m <- 3  # mom
#
# d <- 1  # gumbel
# m <- 3  # mom
#
# keep <- c()
# gof = "KS"
#
# for (st in seq(along = station$index)) {
#   if (station$length_rec[st] > 59) {
#         keep <- c(keep, st)
#       }
# }
#
# temp2plot_raw <- array(NA, dim = c(length(keep), length(sampling_years)))
# for (st in seq(along = keep)) {
#   temp2plot_raw[st, ] <- var.get.nc(gof_nc, "KS", start = c(keep[st], d, m, 1, 1),
#                                 count = c(1, 1, 1, length(sampling_years), 1))
# }
#
# temp2plot <- array(NA, dim = length(sampling_years))
# temp2plot <- colSums(temp2plot_raw, na.rm = TRUE) / length(keep)
# plot(temp2plot, col = "blue")
# par(new = TRUE)


