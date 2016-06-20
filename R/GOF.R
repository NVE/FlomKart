# GOF.R
# MOST OF THOSE FUNCTIONS ARE SUPERCEDED BY FUNCTIONS IN F4NC.R

# Gathers some the goodness of fit functions 

gof_rlevel <- function(dat.full, dat.sample, param.full, param.sample, distr = "distr") {
# Divergence in return levels

  r.periods <- c(50, 100, 200, 500)

  # Computing return levels for the specified return periods
  if(distr == 'gumbel') { 
    r.levels.full <- invF.gumb(F = (1 - 1 / r.periods), param.full[1], param.full[2])
    r.levels.sample <- invF.gumb(F = (1 - 1 / r.periods), param.sample[1], param.sample[2])
  }
  if(distr == 'gamma') { 
    r.levels.full <- qgamma(p = (1 - 1 / r.periods), shape=param.full[1], scale =1/param.full[2])
    r.levels.sample <- qgamma(p = (1 - 1 / r.periods), shape=param.sample[1], scale =1/param.sample[2])
  }
  if(distr == 'gev') {
    r.levels.full <- invF.GEV(F = (1 - 1 / r.periods), param.full[1], param.full[2], param.full[3])
    r.levels.sample <- invF.GEV(F = (1 - 1 / r.periods), param.sample[1], param.sample[2], param.sample[3])
  }
  if(distr == 'gl') {
    r.levels.full <- invF.genlogis(F = (1 - 1 / r.periods), param.full[1], param.full[2], param.full[3])
    r.levels.sample <- invF.genlogis(F = (1 - 1 /r.periods), param.sample[1], param.sample[2], param.sample[3])
  }
  if(distr == 'pearson') {
    r.levels.full <- invF.gamma(F = (1 - 1 / r.periods), param.full[1], param.full[2], param.full[3])
    r.levels.sample <- invF.gamma(F = (1 - 1 / r.periods), param.sample[1], param.sample[2], param.sample[3])
  }
  # Relative difference in %, compared to the density function of the full dataset 
  div.vect <- na.omit(abs(r.levels.full - r.levels.sample) / r.levels.full) * 100
  
  div <- list(max = c(), mean = c())
  
  # We take the max and the mean over the whole distribution
  div$max <- max(div.vect)
  div$mean <- mean(div.vect)
  invisible(div)
}


gof_div <- function(dat.full, dat.sample, param.full, param.sample, distr = "distr") {
# Divergence in density

  xmax <- max(dat.full)
  # Density fuction for the full and reduced datasets
  if(distr == 'gumbel') { 
    density.full <-   dgumbel(0:xmax, param.full[1], param.full[2])
    density.sample <- dgumbel(0:xmax, param.sample[1], param.sample[2])
  }
  if(distr == 'gamma')  { 
    density.full <-   dgamma(0:xmax, param.full[1], param.full[2])
    density.sample <- dgamma(0:xmax, param.sample[1], param.sample[2])
  }
  if(distr == 'gev')    {
    density.full <-   dgev(0:xmax, param.full[3], param.full[1], param.full[2])
    density.sample <- dgev(0:xmax, param.sample[3], param.sample[1], param.sample[2])
  }
  if(distr == 'gl')    {
    density.full <-   f.genlogis(0:xmax, param.full[1], param.full[2], param.full[3])
    density.sample <- f.genlogis(0:xmax, param.sample[1], param.sample[2], param.sample[3])
  }
  if(distr == 'pearson')    {
    density.full <-   f.gamma(0:xmax, param.full[1], param.full[2], param.full[3])
    density.sample <- f.gamma(0:xmax, param.sample[1], param.sample[2], param.sample[3])
  }
  
  # Relative difference in %, compared to the density function of the full dataset 
  div.vect <- na.omit(abs(density.full - density.sample) / density_full) * 100
  div <- list(max = c(),mean = c())
  # We take the max and the mean over the whole distribution
  div$max <- max(div.vect)
  div$mean <- mean(div.vect)
  invisible(div)
}

gof_stability <- function(dat.full, param.full, k, distr = "distr", method = "method") {
# sampling function calling the stability measures
  
  L <- length(dat.full) - 1         # Length of initial dataset
  stab <- list(max = c(), mean = c())

  for (j in k:L) {  # i < j so i:j is the NA window we create
    i <- j - k + 1
    dat.sample <- dat.full
    dat.sample[i:j] <- NA  # Creation of the NA window and na.omit to create a reduced dataset
    dat.sample <- as.vector(na.omit(dat.sample))

    # We paste the parameters "distr" and "method" to assimilate the right function to the temp FUN function
    FUN <- match.fun(paste(as.character(distr), "_", as.character(method), collapse = "", sep = ""))
    param.sample <- FUN(dat.sample)$estimate
    if (is.na(param.sample[1]) == TRUE) {
      stab$max[j] <- NA
      stab$mean[j] <- NA
      next(j)  # We could add a warning too
    }
    # Fill vectors with max and min divergence (density function or return level)
    temp <- gof_rlevel(dat.full, dat.sample, param.full, param.sample, 
                        as.character(distr))
    stab$max[j] <-  temp$max
    stab$mean[j] <-  temp$mean
  }
  stab$max <- max(na.omit(stab$max)) 
  stab$mean <- mean(na.omit(stab$mean))
  invisible(stab)
}

gof_area <- function(dat, station.nb.vect) {
# Calculation of the area between the uniform distribution ECDF and the FF distribution ECDDF
# According to Blanchet.

  distributions = c("gumbel","gamma","gev","gl","pearson")
  fitting_methods = c("mle","Lmom","mom")
  method = "mom"
  i <- 1
  for ( distr in distributions ) {
    # for (method in fitting_methods) {
      print( paste( "fitting parameters for ", distr))
      
      ff <- run_all_ff(dat, station.nb.vect, distr, method)
      sorted.ff.1 <- sort(ff$ff1)
      sorted.ff.2 <- sort(ff$ff2)
      uniform.1 <- seq(0, 1, length.out = length(sorted.ff.1))
      uniform.2 <- seq(0, 1, length.out = length(sorted.ff.2))
      
#       area.1$distr <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
#       area.2$distr <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
#       area.1.abs$distr  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
#       area.2.abs$distr <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
      area.1[i] <- sum(sorted.ff.1 - uniform.1) / 2 / length(sorted.ff.1)
      area.2[i] <- sum(sorted.ff.2 - uniform.2) / 2 / length(sorted.ff.2)
      area.1.abs[i]  <- sum(abs(sorted.ff.1 - uniform.1)) / 2 / length(sorted.ff.1)
      area.2.abs[i] <- sum(abs(sorted.ff.2 - uniform.2)) / 2 / length(sorted.ff.2)
      i <- i + 1
    # }
  }
results <- list(area1 = area.1, area1abs = area.1.abs,  area2 = area.2, area2abs = area.2.abs)
return(results)
}


gof_aic <- function(dat, param, distr = "distr") {
# Aikake Information Criterion. TO IMPLEMENT
    
}


gof_bic <- function(dat, param, distr = "distr") {
# Bayesian Information Criterion. TO IMPLEMENT 
}