# Seasonal plot for report and other plotting experiments with ggplot

# TEST
station.nb.vect <- unique(dat$snumber)
# dfx <- data.frame(family = seq(1,1, length.out = 4* length(station.nb.vect)), 
#                   item = c(station.nb.vect, station.nb.vect, station.nb.vect, station.nb.vect), 
#                   score = c( rep("winter", length(station.nb.vect)), rep("spring", length(station.nb.vect)), 
#                              rep("summer", length(station.nb.vect)), rep("fall", length(station.nb.vect))),
#                   value = rep(0, 4 * length(station.nb.vect)) )

dfx <- data.frame(family = c(seq(1, length(station.nb.vect), length.out = length(station.nb.vect)),
                             seq(1, length(station.nb.vect), length.out = length(station.nb.vect)),
                             seq(1, length(station.nb.vect), length.out = length(station.nb.vect)),
                             seq(1, length(station.nb.vect), length.out = length(station.nb.vect))
), 
item = c(station.nb.vect, station.nb.vect, station.nb.vect, station.nb.vect), 
score = c( rep("winter", length(station.nb.vect)), rep("spring", length(station.nb.vect)), 
           rep("summer", length(station.nb.vect)), rep("fall", length(station.nb.vect))),
value = rep(0, 4 * length(station.nb.vect)) )

for (i in seq(along = station.nb.vect)) {
  sdat <- sdat_load(dat, station.nb.vect[i]) 
  mydate = as.POSIXlt(as.character(sdat$date_DOGN))
  sdat$month = mydate$mon + 1
  
  dfx$value[i] <- length(which(sdat$month < 4)) / length(sdat$month)
  dfx$value[i + length(station.nb.vect)] <- length(which(sdat$month > 3 & sdat$month < 7)) / length(sdat$month)
  dfx$value[i + 2 * length(station.nb.vect)] <- length(which(sdat$month > 6 & sdat$month < 10)) / length(sdat$month)
  dfx$value[i + 3 * length(station.nb.vect)] <- length(which(sdat$month > 9)) / length(sdat$month)
  
}

dfx <- na.omit(dfx)

indices.spring <- which(dfx$score == "spring")
indices.summer <- which(dfx$score == "summer")
indices.fall <- which(dfx$score == "fall")
indices.winter <- which(dfx$score == "winter")

spring_sorted <- sort(dfx$value[indices.spring],  index.return = TRUE)
summer_sorted <- dfx$value[spring_sorted$ix + 2*440]
fall_sorted <- dfx$value[spring_sorted$ix + 3*440]
winter_sorted <- dfx$value[spring_sorted$ix]

h.spring <- cut(spring_sorted$x, breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))  
h.summer <- cut(summer_sorted, breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) 
h.fall <- cut(fall_sorted, breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))  
h.winter <- cut(winter_sorted, breaks = c(-0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) 

h.spring <- cut(spring_sorted$x, breaks = c(-0.001, 0.25, 0.5, 0.75, 1))  
h.summer <- cut(summer_sorted, breaks = c(-0.001, 0.25, 0.5, 0.75, 1)) 
h.fall <- cut(fall_sorted, breaks = c(-0.001, 0.25, 0.5, 0.75, 1))  
h.winter <- cut(winter_sorted, breaks = c(-0.001, 0.25, 0.5, 0.75, 1)) 

h.spring <- cut(spring_sorted$x, breaks = c(-0.001, 0.05, 0.5, 0.95, 1))  
h.summer <- cut(summer_sorted, breaks = c(-0.001, 0.05, 0.5, 0.95, 1)) 
h.fall <- cut(fall_sorted, breaks = c(-0.001, 0.05, 0.5, 0.95, 1))  
h.winter <- cut(winter_sorted, breaks = c(-0.001, 0.05, 0.5, 0.95, 1)) 



test <- data.frame(station = seq(1,1,440), spring = h.spring, summer = h.summer, fall = h.fall, winter = h.winter)
ggplot(data = test) + 
  geom_bar(mapping = aes(x = test$spring, fill = test$fall), position = "fill", binwidth = 1) +
  coord_polar(theta = "y") +
  ggtitle('Position = "fill"')

ggplot(data = test) + 
  geom_bar(mapping = aes(x = 1, fill = test$winter), position = "fill", binwidth = 0.8) +
  geom_bar(mapping = aes(x = 2.5, fill = test$spring), position = "fill", binwidth = 0.8) +
  geom_bar(mapping = aes(x = 4, fill = test$summer), position = "fill", binwidth = 0.8) +
  geom_bar(mapping = aes(x = 6, fill = test$fall), position = "fill", binwidth = 0.8) +
  # coord_polar(theta = "y") +
  ggtitle('Seasonal flood ratios')

####################################### THIS IS A GOOD PLOT
# test as box plot 
ggplot(data = test, mapping = aes(x = Season, y = Ratio)) + 
  geom_boxplot(mapping = aes(x = "Spring", y = spring_sorted$x),  binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "Summer", y = summer_sorted), binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "Fall", y = fall_sorted), binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "Winter", y = winter_sorted), binwidth = 0.8) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) + 
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  theme(axis.text.y = element_text(size = rel(1.8), angle = 90)) + 
  theme(axis.text.x = element_text(size = rel(1.8), angle = 00)) 
# coord_polar(theta = "y") +
# labs(x = "Season", y = "Ratio")
ggtitle('Seasonal flood ratios')

########################################## TEST OF BOX PLOT FOR ALL STATIONS THAT HAVE MORE THAN N data
# Works well without GGPLOT!
min_data <- 50
max_data <- 80

# sdat1 <- sdat_load(dat, 200604)              # This is Elverum, the longest dataset
# sdat2 <- sdat_load(dat, 6200005) 
flood_list <- list()
station_list <- c()
for (i in 1:length(station.nb.vect)) {
  # assign(paste("sdat",i,sep = ""), sdat_load(dat, station.nb.vect[i]) )
  if (!is.null(sdat_load(dat, station.nb.vect[i]))) {
    
    flood_list[[i]] <- sdat_load(dat, station.nb.vect[i])$flom_DOGN / mean(na.omit(sdat_load(dat, station.nb.vect[i])$flom_DOGN))
    station_list[i] <- as.character(sdat_load(dat, station.nb.vect[i])$name[1])
  } else {
    flood_list[[i]] <- NA
    station_list[i] <- NA
  }
}
to_plot <- c()
k <- 0
for (i in 1:length(station.nb.vect)) {
  
  if (length(flood_list[[i]]) > min_data && length(flood_list[[i]]) < max_data) {
    k <- k+1
    to_plot[k] <- i
  }
}

par(las = 1) # all axis labels horizontal
boxplot(flood_list[to_plot], names = station_list[to_plot], horizontal = TRUE, outwex = TRUE)

max_vect <- colMaxs(flood_list[to_plot])

boxplot(sort(flood_list[to_plot], decreasing = TRUE), horizontal = TRUE)

     
#####################################################################################################




ggplot(mapping = aes(x = "Station", y = "Flood_data")) +
  geom_boxplot(mapping = aes(x = paste(station_list[to_plot[1]]), y = flood_list[to_plot[1]]),  binwidth = 0.8)

+
  
  geom_boxplot(mapping = aes(x = "prout2", y = sdat2$flom_DOGN), binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "prout3", y = sdat25$flom_DOGN), binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "prout4", y = sdat26$flom_DOGN), binwidth = 0.8) +
  geom_boxplot(mapping = aes(x = "prout5", y = sdat50$flom_DOGN), binwidth = 0.8) 


#########################

g <- ggplot(mapping = aes(x = "Station", y = "Flood_data")) 
g <- g+  geom_boxplot(mapping = aes(x = "prout", y = sdat1$flom_DOGN),  binwidth = 0.8) 

g <- g+  geom_boxplot(mapping = aes(x = "prout2", y = sdat2$flom_DOGN), binwidth = 0.8) 
g
# plots <- list()
k <- 0
temp_dat <- list()
# seq(along = station.nb.vect)
g <- ggplot(temp_dat, mapping = aes(x = "Station", y = "Flood_data"))

for (i in 1:20) {
  sdat <- sdat_load(dat, station.nb.vect[i])
  if (length(sdat$flom_DOGN) > min_data) {
    #     plots[[k]]  <- ggplot(mapping = aes(x = Station, y = Flood_data)) +
    #       geom_boxplot(mapping = aes(x = k, y = sdat$flom_DOGN),  binwidth = 0.8) + 
    #       xlim(0,10)
    k <- k + 1
    temp_dat[[k]] <- data.frame(y = sdat$flom_DOGN)  
    g <- g +
      geom_boxplot(mapping = aes(x = paste("prout", k, sep = ""),
                                 y = temp_dat[[k]]),  binwidth = 0.8)  
    # geom_jitter(mapping = aes(x = k, y = sdat$flom_DOGN)) +
    
    
    # geom_bar(aes(fill = k))
    # facet_grid(.~ k)
    
    
    print(sdat$name[1])
    print(length(sdat$flom_DOGN))
  }
}
# g <- g+ facet_wrap(~mix, scales="free_x") 
g

g <- ggplot(mapping = aes(x = "Station", y = "Flood_data")) +
  
  geom_boxplot(mapping = aes(x = paste("prout", k, sep = ""),
                             y = temp_dat),  binwidth = 0.8) 

# ggtitle('Seasonal flood ratios')
g


# Modify lst into data frames of varying dimension
temp_dat2 <- lapply(series, function(x) {
  data.frame(series = factor(x, levels = series),
             y = temp_dat[[x]])
})

# Stack the data frames
lst <- do.call(rbind, lst)

# Make the plot
ggplot(lst, aes(x = series, y = y)) +
  geom_boxplot()








print(k)
print(plots)
plots
g

ggtitle('Boxplot for stations with more than 80 years of data')

##################################################################
# TRY WITHOUT GGPLOT
min_data <- 50

sdat1 <- sdat_load(dat, 200604)              # This is Elverum, the longest dataset
sdat2 <- sdat_load(dat, 6200005) 
# test <- sdat_load(dat, station.nb.vect[445]) 

ggplot(mapping = aes(x = "Station", y = "Flood_data")) +
  geom_boxplot(mapping = aes(x = "prout", y = sdat1$flom_DOGN),  binwidth = 0.8) +
  
  geom_boxplot(mapping = aes(x = "prout2", y = sdat2$flom_DOGN), binwidth = 0.8) 

########################################################

plot(h.spring)
plot(h.summer)
plot(h.fall)
plot(h.winter)


# Doesnt work
plot(h.spring + h.summer + h.fall + h.winter,
     ylim = c(0, 1), type = "h", lwd = 2, col = "black", xlab = "", ylab = "")

plot(winter_sorted + fall_sorted + summer_sorted + spring_sorted$x,
     ylim = c(0, 1), type = "h", lwd = 2, col = "black", xlab = "", ylab = "")
par(new = TRUE)
plot(fall_sorted + summer_sorted + spring_sorted$x,
     ylim = c(0, 1), type = "h", lwd = 2, col = "red", xlab = "", ylab = "")
par(new = TRUE)
plot(summer_sorted + spring_sorted$x,
     ylim = c(0, 1), type = "h", lwd = 2, col = "blue", xlab = "", ylab = "")
par(new = TRUE)

plot(spring_sorted$x, 
     ylim = c(0, 1), type = "h", lwd = 2, col = "grey", xlab = "stations", ylab = "seasonal ratio")



# indices.summer <- which(dfx$score == "summer")
# temp <- dfx$value[indices.summer]
# indices.summer <- spring_sorted$ix + 440


# indices.winter <- which(dfx$score == "winter")
# plot(dfx$value[indices.spring] + dfx$value[indices.summer] + dfx$value[indices.fall] + dfx$value[indices.winter])
# hist(dfx$value[indices.winter])