# MAIN.R
# A VARIETY OF FUNCTION CALL EXPERIMENTS
# This file does not contain function definitions. It only calls functions.
# Code to fit various distributions (gamma, gev, Generalized logistics, gumbel) to flood data
# with a choice of estimators (ordinary moments, L-moments, maximum likelihood, Bayesian estimation)
#  - graphs of i vs T for each duration on a multi-panel plot
#  - plot of i vs. D for a range of return periods
#
# 2015-Nov Florian Kobierska

##############################################################################################
rm(list = ls())  # clean out workspace and set working directory
# setwd('//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/R_Scripts')

# Load libraries
library(plotrix)      # For specific plotting functions like the addtable2plot
library(evd)          # Functions for extreme value distributions (GEV and GUMBEL use this package)
library(evir)         # Extreme values in R
library(ismev)        # An Introduction to Statistical Modeling of Extreme Values
library(pracma)       # Practical Numerical Math routines. We use newtonRaphson which actually bugs quite often
library(evdbayes)     # Bayesian exreme value analysis
library(nsRFA)        # Includes l-moment calculations par.gev (GENLOGIS, GENPAR, invF.gumb, invF.gamma for Pearson)
library(fBasics)      # Some basic statistical functions
library(stats4)       # for mle()
library(MASS)         # for fitdistr() 
library(zoo)
library(glogis)       # Generalized logistics
library(fitdistrplus)
library(PearsonDS)    # Pearson distribution TO CHECK WHERE IT\S USED!!!!!
library(goftest)      # For ad.test
# ks.test comes from stats package

# library(rJava)        # Installed with iplots http://www.r-statistics.com/2012/08/how-to-load-the-rjava-package-after-the-error-java_home-cannot-be-determined-from-the-registry/#comment-2532
library(png)          # Installed with iplots
# library(iplots)       # For interactive plots where points are highlighted in several graphs  http://www.r-statistics.com/tag/iplots/
library(plyr)         # For the failwith(NA, f) function!
library(verification) # For goodness of fit tools
library(ggplot2)      # Hadleys plot package
library(magrittr)     # For piping functions

# Load libraries for the map function
library(maptools)     # Required for map
library(rgdal)        # Required for map
library("sp","rgdal") # Required for map
library(RColorBrewer) # Required for map

# To check if really required
require(Hmisc)
require(aqfig)
library(calibrate)

# Prob not required for map
library(timeDate) #
library(data.table) # 
require(Kendall) #
require(hydroGOF) #

##############################################################################################
# global variables for the minimum number of years
GLOBAL_min_years_data <<- 30

# Source all the required R files
source('LOAD.r')
source('RUN_ALL.r')
source('GOF.r')           # Goodness of fit functions
source('PLOTS.r')         # Plotting functions common to all distributions (No maps)
source('EXPERIMENTAL.r')
source('EXPLORE.r')       # Plotting the map of Norway with various parameters
source('PLOTS_4_REPORT.r')  # Where functions for specific plots are defined

# Fitting the distributopns
source('GUMBEL.r')
source('GL.r')
source('GEV.r')
source('GAMMA.r')
source('PEARSON.r')

##############################################################################################
# Load flood data
dat <- read.csv("//nve/fil/h/HM/Interne Prosjekter/Flomkart/Model_fitting/Florian/Data/AMS_table_updated.csv", sep=";")
# Create a vector with all the different station numbers, in order to run every station separately in one go
station.nb.vect <- unique(dat$snumber)

dat <- save_coordinates(dat, station.nb.vect)

# Load data either with the station number vector or directly with the number of the station  
sdat <- sdat_load(dat, station.nb.vect[20]) # 20 is good example
sdat <- sdat_load(dat, 200011)
sdat <- sdat_load(dat, 6200005)             # This is Bulken, a long record
sdat <- sdat_load(dat, 200604)              # This is Elverum, the longest dataset

# Some plots of the raw and normalized data TO INTEGRATE INTO PLOTS:R
sdat_plot(sdat$flom_DOGN)
sdat_plot(sdat$normQ)
plot(sdat$year, sdat$flom_DOGN)  # Plot flood flow evolution over the years

# TO TRY
# library(animation)
# saveGIF({for(i in 1:n) plot(f[[i]])})

# TEST OF LEAFLET MAP
norway_map2(dat,station.nb.vect)


##############################################################################################
# Basic simulation and plot for one specific station
distr <- "gamma"
method <- "mom"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)        
temp.r.levels <- RLEVELS4NC(param$estimate, return.periods, distr = "distr")
test <- as.numeric(gof_cs(sdat$flom_DOGN, param$estimate, distr))
print(test)
GOF.list <- gof_all(sdat$flom_DOGN, param$estimate, distr)
plot_all(sdat$flom_DOGN, GOF.list, param, distr = distr, method = method)

##############################################################################################
## Normalized data DOESNT WORK TO DISCUSS
param <- fit_distr(sdat$normQ, distr = distr, method = method)        
GOF_list <- gof_all(sdat$normQ, param$estimate, distr)
plot_all(sdat$normQ, GOF.list, param, distr, method)

## LOGPEARSON: need to sort out plot_all for this case
param <- pearson_mle(sdat$logQ)
GOF.list <- gof_all(sdat$logQ, param, distr = "pearson")  
# Strange bug in GOF "F.gamma(0.1 * i * xmax, param[1], param[2], param[3]) : (list) object cannot be coerced to type 'double'
# Could be a scale PB as in the general plot
plot_density(sdat$logQ, GOF.list, param, distr = "pearson")  # individually works
plot_ecdf(sdat$logQ, param, distr = "pearson")  # individually works
plot_qq(sdat$logQ, param, distr = "pearson")  # individually works
plot_rlevel(sdat$logQ, param, distr = "pearson")  # individually works

plot_all(sdat$normQ, GOF.list, param, distr = "pearson", mle)  # DOESNT WORK BECAUSE OF SCALE PB

##############################################################################################
# Run all the stations with a specific distribution and export a table with all results
x.gumbel <- run_all(dat, station.nb.vect, distr = "gumbel", method = "mom") 
x.gamma <- run_all(dat, station.nb.vect, distr = "gamma", method = "mle")     # works using Lmom, some pbs with MLE at station 409
x.gev <- run_all(dat, station.nb.vect, distr = "gev", method  = "mom") 
x.gl <- run_all(dat, station.nb.vect, distr = "gl", method ="mle")            # works using mom, some pbs with MLE (station 101) and Lmom (station 89)
x.pearson <- run_all(dat, station.nb.vect, distr = "pearson", method = "mom") # works using mom, some pbs with MLE (station 217) and Lmom (station 53)

x.gumbel.norm <- run_all_norm(dat, station.nb.vect, distr = "gumbel", method = "mle") 

# Run all the stations with all distributions plot the best distribution with a specific GOF
KS <- run_all_ks(dat, station.nb.vect, method = "mle")
plot_gof(dat, station.nb.vect, method = "mle", input = KS)

AD <- run_all_ad(dat, station.nb.vect, method = "mom")
plot_gof(dat, station.nb.vect, method = "mom", input = AD)

AD <- run_all_ad(dat, station.nb.vect, method = "Lmom")
plot_gof(dat, station.nb.vect, method = "Lmom", input = AD)

plot_all_stations(dat,station.nb.vect, distr = "gumbel", method = "mle")

distr <- "pearson"
method <- "mle"
ff <- run_all_ff(dat, station.nb.vect, distr, method)
sorted.ff.1 <- sort(ff$ff1)
sorted.ff.2 <- sort(ff$ff2)
uniform.1 <- seq(0, 1, length.out = length(sorted.ff.1))
uniform.2 <- seq(0, 1, length.out = length(sorted.ff.2))
area.1 <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
area.2 <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
area.1.abs  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
area.2.abs <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
plot(ecdf(ff$ff1))
plot(ecdf(ff$ff2))

test <- gof_area(dat, station.nb.vect)



##### test gives same result for all tests
distr <- "gamma"
method <- "mom"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)      
pred <- dgamma(sdat$flom_DOGN, param$estimate[1], param$estimate[2])
test1 <- quantileScore(sdat$flom_DOGN, pred, 0.9, c(0, 1000, 2000, 3000, 4000))$qs.reliability

distr <- "gumbel"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)    
pred <- dgumbel(sdat$flom_DOGN, param$estimate[1], param$estimate[2])
test2 <- quantileScore(sdat$flom_DOGN, pred, 0.9, c(0, 1000, 2000, 3000, 4000))$qs.reliability

distr <- "gev"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
pred <- dgev(sdat$flom_DOGN, param$estimate[3], param$estimate[1], param$estimate[2])
test3 <- quantileScore(sdat$flom_DOGN, pred, 0.9, c(0, 1000, 2000, 3000, 4000))$qs.reliability
              
##############################################################################################
# TEST stability plots
plot_all_stations(dat, station.nb.vect, distr = "gamma", method = "Lmom")
legend("topright", inset = .05, c("gamma","gl" ), col = c("blue","red"),
       lty = c(1,1),lwd = c(2,3), merge = TRUE, bg = "gray90")

# TEST reliability plots
plot_all_stations(dat, station.nb.vect, distr = "gamma", method = "Lmom")
legend("topright", inset = .05, c("gamma","gl" ), col = c("blue","red"),
       lty = c(1,1),lwd = c(2,3), merge = TRUE, bg = "gray90")

##############################################################################################
# TEST OF SAMPLING FUNCTION
k <- 5
dat.min10 <- load_1sample(sdat$flom_DOGN, k, "min")
param <- fit_distr(dat.min10, distr = distr, method = method)        
GOF.list <- gof_all(dat.min10, param$estimate, distr)
plot_all(dat.min10,GOF.list, param, distr, method)

GOF_list <- gof_all(sdat$flom_DOGN, param$estimate, distr)
plot_all(sdat$flom_DOGN, GOF.list, param, distr, method)

dat.max10 <- load_1sample(sdat$flom_DOGN, k, "max")
param <- fit_distr(dat.max10, distr = distr, method = method)        
GOF.list <- gof_all(dat.max10, param$estimate, distr)
plot_all(dat.max10,GOF.list, param, distr, method)

GOF.list <- gof_all(sdat$flom_DOGN, param$estimate, distr)
plot_all(sdat$flom_DOGN,GOF.list, param, distr, method)

##############################################################################################
# RELIABILITY CALCULATIONS WITH BS (BRIER SCORE)
lim <- length(sdat$flom_DOGN) - min_years_2param - 1
L <- length(sdat$flom_DOGN)
x.axis <- c()

distr <- "gamma"
method <- "Lmom"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)
div <- list(max = c(), mean = c())
BS <- c() 
for (k in 1:lim)  {

  x.axis[k] <- L - k
  BS[k] <- gof_bs(sdat$flom_DOGN, param, k, distr, method)  
}
ymin <- min(BS)*1.5
ymax <- max(BS)*1.5
windows()
plot(x.axis, BS, col = "blue", ylim = c(ymin, ymax))  
par(new = TRUE)

distr <- "gl"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  
  x.axis[k] <- L - k
  BS[k] <- gof_bs(sdat$flom_DOGN, param, k, distr, method)  
}
plot(x.axis, BS, col = "black", ylim = c(ymin, ymax))  
legend("topright", inset = .05, c("gamma","gl"), col = c("blue","black"),
       lty = c(1,1),lwd = c(2,3), merge = TRUE, bg = "gray90")

##############################################################################################
# RELIABILITY CALCULATIONS WITH QS (QUANTILE SCORE)
lim <- length(sdat$flom_DOGN) - min_years_2param - 1
L <- length(sdat$flom_DOGN)
x.axis <- c()

distr <- "gamma"
method <- "mle"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)
div <- list(max = c(), mean = c())
QS <- c() 
for (k in 1:lim)  {
  
  x.axis[k] <- L - k
  QS[k] <- gof_qs(sdat$flom_DOGN, param, k, distr, method)  
}
ymin <- min(QS)*1.5
ymax <- max(QS)*1.5
windows()
plot(x.axis, QS, col = "blue", ylim = c(ymin, ymax))  
par(new = TRUE)

distr <- "gev"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  
  x.axis[k] <- L - k
  QS[k] <- gof_qs(sdat$flom_DOGN, param, k, distr, method)  
}
plot(x.axis, QS, col = "black", ylim = c(ymin, ymax))  
legend("topright", inset = .05, c("gamma", "gl"), col = c("blue", "black"),
       lty = c(1,1),lwd = c(2,3), merge = TRUE, bg = "gray90")

##############################################################################################
# STABILITY CALCULATIONS
## PLOT WITH DIFFERENT ESTIMATION METHODS -> mle better for gamma, gumbel and pearson, Lmom better for gev, mom better for gl (mle extremely slow)
lim <- length(sdat$flom_DOGN) - min_years_2param - 1
L <- length(sdat$flom_DOGN)
x.axis <- c()

distr <- "gev"
method <- "mle"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)
div <- list(max = c(), mean = c())
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean  
  x.axis[k] <- L - k
  print(k)
}

windows()
ymax <- max(div$mean) * 1.5
plot(x.axis, div$mean, col = "blue", ylim = c(0, ymax),xlab = "", ylab = "")
par(new = TRUE)

method <- "Lmom"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean  
  x.axis[k] <- L - k
}
plot(x.axis,div$mean,col = "black", ylim = c(0,ymax),xlab = "", ylab = "")
par(new = TRUE)

method <- "mom"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean   
  x.axis[k] <- L - k
}
plot(x.axis, div$mean, col = "red", ylim = c(0, ymax),
     main = "Percentage variation in return levels (50, 100, 200, 500 years) for the gev distr", 
     xlab = "Number of data points", ylab = "Variation in %")
legend("topright", inset = .05, c("mle", "Lmom", "mom" ),  col = c("blue", "black", "red"),
        lty = c(1, 1),lwd = c(2, 3), merge = TRUE, bg = "gray90")

##############################################################################################
# PLOT WITH DIFFERENT DISTRIBUTIONS
# gamma always better than gumbel, gl maybe better than gev, gamma is close
method <- "mom"
distr <- "gl"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean  
  x.axis[k] <- L - k
}
windows()
ymax <- max(div$mean) * 1.1
plot(x.axis,div$mean,col = "blue", ylim = c(0,ymax))
par(new = TRUE)

distr <- "gev"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean  
  x.axis[k] <- L - k
}
plot(x.axis, div$mean, col = "black", ylim = c(0,ymax))
par(new = TRUE)

distr <- "gamma"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
for (k in 1:lim)  {
  temp <- gof_stability(sdat$flom_DOGN, param$estimate, k, distr, method)
  div$max[k] <- temp$max
  div$mean[k] <- temp$mean  
  x.axis[k] <- L - k
}
plot(x.axis,div$mean, col = "red", ylim = c(0,ymax))
legend("topright", inset = .05, c("gl", "gev", "gamma" ), col = c("blue", "black", "red"),
       lty = c(1, 1),lwd=c(2, 3), merge = TRUE, bg = "gray90")

##############################################################################################
# Bayesian estimation: to double check, improve and integrate as new method into "fit_distr" function
## gev
sdat <- sdat_load(dat, station.nb.vect[20]) # 20 is good example
source('PLOTS.r')

distr <- "gamma"
method <- "bayes"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method)      
GOF.list <- gof_all(sdat$flom_DOGN, param$estimate, distr)
plot_all(sdat$flom_DOGN, GOF.list, param, distr = distr, method = method)

method <- "mle"
param <- fit_distr(sdat$flom_DOGN, distr = distr, method = method) 
GOF.list <- gof_all(sdat$flom_DOGN, param$estimate, distr)
plot_all(sdat$flom_DOGN, GOF.list, param, distr = distr, method = method)


gev.bayes <- gev_bayes_OLD(as.vector(sdat$flom_DOGN))
plot(gev.bayes)
mupars <- as.vector(gev.bayes$parameters[, 1, 1:3])
spars <- as.vector(gev.bayes$parameters[, 2, 1:3])
kpars <- as.vector(gev.bayes$parameters[, 3, 1:3])
mmrp <<- c(100)
test <- get_posterior_gev(sdat$flom_DOGN, mupars, spars, kpars)


## gl
gl.bayes <- gl_bayes(as.vector(sdat$flom_DOGN))
plot(gl.bayes)

mupars <- as.vector(gl.bayes$parameters[, 1, 1:3])
spars <- as.vector(gl.bayes$parameters[, 2, 1:3])
kpars <- as.vector(gl.bayes$parameters[, 3, 1:3])
mmrp <<- c(100)
test <- get_posterior_gl(sdat$flom_DOGN, mupars, spars, kpars)

## pearson OK, but not with prior
pearson.bayes <- pearson_bayes(as.vector(sdat$flom_DOGN))
plot(pearson.bayes)

mupars <- as.vector(pearson.bayes$parameters[, 1, 1:3])
spars <- as.vector(pearson.bayes$parameters[, 2, 1:3])
kpars <- as.vector(pearson.bayes$parameters[, 3, 1:3])
mmrp <<- c(100)
test <- get_posterior_pearson(sdat$flom_DOGN, mupars, spars, kpars)

## gumbel
gumbel.bayes <- gumbel_bayes(as.vector(sdat$flom_DOGN))
plot(gumbel.bayes)

mupars <- as.vector(gumbel.bayes$parameters[, 1, 1:3])
spars <- as.vector(gumbel.bayes$parameters[, 2, 1:3])
mmrp <<- c(100)
test <- get_posterior_gumbel(sdat$flom_DOGN, mupars, spars) 
plot(test)


##############################################################################################
# Direct copy of Kolbkorn for gev

myprior < - dnorm(dat,0,0.2)

  fit_gev_bayes_s1 <- BayesianMCMC(na.omit(set1[,i]), nbpas=5000, nbchaines=3, confint=c(0.05, 0.95), dist="gev",apriori=myprior)
  fit_gev_bayes_s2 <- BayesianMCMC(na.omit(set2[,i]), nbpas=5000, nbchaines=3, confint=c(0.05, 0.95), dist="gev",apriori=myprior)
  
  mupars1 <- as.vector(fit_gev_bayes_s1$parameters[, 1, 1:3])
  mupars2 <- as.vector(fit_gev_bayes_s2$parameters[, 1, 1:3])
  
  spars1 <- as.vector(fit_gev_bayes_s1$parameters[, 2, 1:3])
  spars2 <- as.vector(fit_gev_bayes_s2$parameters[, 2, 1:3])
  
  kpars1 <- as.vector(fit_gev_bayes_s1$parameters[, 3, 1:3])
  kpars2 <- as.vector(fit_gev_bayes_s2$parameters[, 3, 1:3])
  
  #  Calculate quantiles
  ss1_post_rl <- get_posterior_gev(myrp, mupars1, spars1, kpars1)
  ss1_ml_rl <- invF.gev(F = (1 - 1 / myrp), fit_gev_bayes_s1$parametersML[1], fit_gev_bayes_s1$parametersML[2],
                        fit_gev_bayes_s1$parametersML[3])
  
  ss2_post_rl<-get_posterior_gev(myrp,mupars2,spars2,kpars2)
  ss2_ml_rl<-invF.gev(F=(1-1/myrp),fit_gev_bayes_s2$parametersML[1],fit_gev_bayes_s2$parametersML[2],fit_gev_bayes_s2$parametersML[3])
  ### Calculate span
  
  #NB. These have to be put into a matrix!!!!!
  span_quantiles_gev_post[i,] <- 2 * abs(ss1_post_rl - ss2_post_rl) / (ss1_post_rl + ss2_post_rl) 
  span_quantiles_gev_ml[i,] <- 2 * abs(ss1_ml_rl - ss2_ml_rl) / (ss1_ml_rl + ss2_ml_rl) 
  
  # Get the predictive distribution approximated by spline
  emprp<-c(1.01, 1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 4, 8, 10, 15, 20, 50, 70, 100, 150, 200, 300,
           500, 700, 1000, 5000, 10000, 1000000)
  ss1_post_emp <- get_posterior_gev(emprp, mupars1, spars1, kpars1)
  ss2_post_emp <- get_posterior_gev(emprp, mupars2, spars2, kpars2)
  
  postemp1 <- splinefun(x = ss1_post_emp, y = 1-1/emprp,method = "fmm",ties = mean)
  postemp2 <- splinefun(x = ss2_post_emp, y = 1-1/emprp,method = "fmm",ties = mean)
  
  FF12_gev_post[i] <- (min(postemp1(max_set2[i]), 1))^n2[i]
  FF21_gev_post[i] <- (min(postemp2(max_set1[i]), 1))^n1[i]
  FF12_gev_ML[i] <- (F.gev(max_set2[i], fit_gev_bayes_s1$parametersML[1], fit_gev_bayes_s1$parametersML[2],
                           fit_gev_bayes_s1$parametersML[3]))^n2[i]
  FF21_gev_ML[i] <- (F.gev(max_set1[i], fit_gev_bayes_s2$parametersML[1], fit_gev_bayes_s2$parametersML[2],
                           fit_gev_bayes_s2$parametersML[3]) )^n1[i]
  if(is.na(FF12_gev_ML[i])) FF12_gev_ML[i] = 1.0
  if(is.na(FF21_gev_ML[i])) FF21_gev_ML[i] = 1.0
  
  