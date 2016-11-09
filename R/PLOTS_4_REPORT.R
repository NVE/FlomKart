# PLOTS_4_REPORT.r
# This file gathers functions created to produce specific plots to be used in the report
# 1. functions with leaflet
# 2. Distribution of the length of data at stations
# 3. Attempt at a Lmom plot.


# norway_map <- leaflet() %>%
#   addTiles() %>%
#   setView(10.45, 59.57, zoom = 10) %>%
#   addMarkers(10.45, 59.57, popup = 'Oslo') %>%
# addMarkers(12, 61.1, popup = 'Oslo')


#' norway_map
#' @description Figure: MAP OF NORWAY WITH ALL THE STATIONS. WORKS
#' @param dat
#' @param station.nb.vect
#'
#' @return
#' @import leaflet
#' @import magrittr
#' @export
#'
#' @examples
norway_map <- function(dat,station.nb.vect) {

  temp <- list(long = c(), lat = c(), name = c(), width = c())
   for (i in seq(along = station.nb.vect)) {
    j <- which(dat$snumber == station.nb.vect[i])
    w <- length(j) / 20
    j <- j[1]
    if (!is.null(dat$longitude[j]) && !is.null(dat$lat[j])) {
    temp$long <- c(temp$long, dat$longitude[j])
    temp$lat <- c(temp$lat, dat$latitude[j])
    temp$name <- c(temp$name, as.character(dat$name[j]))
    temp$width <- c(temp$width, w)
    } else {next(i)}

  }
  map <- leaflet() %>% addTiles()
  setView(map, 13, 64, zoom = 5)
  # addMarkers(map, data = temp, lng = ~ long, lat = ~ lat, popup = paste(as.character(temp$name)))
  addCircleMarkers(map, data = temp, lng = ~ long, lat = ~ lat,
                   popup = paste(as.character(temp$name)), radius = ~ width,
                   color = "red")
}



#' norway_map2
#' @description MAP OF NORWA WITH COLOR CHANGING ACCORDING TO THE NUMBER OF DATA.
#' @param dat
#' @param station.nb.vect
#'
#' @return
#' @export
#'
#' @examples
norway_map2 <- function(dat,station.nb.vect) {

  library(leaflet)
  library(magrittr)

  temp <- list(long = c(), lat = c(), name = c(), width = c(), color = c())
  map <- leaflet() %>% addTiles()
  setView(map, 13, 64, zoom = 5)


  for (i in seq(along = station.nb.vect)) {

    j <- which(dat$snumber == station.nb.vect[i])
    w <- length(j) / 20
    j <- j[1]
    if (!is.null(dat$longitude[j]) && !is.null(dat$lat[j])) {
      temp$long <- c(temp$long, dat$longitude[j])
      temp$lat <- c(temp$lat, dat$latitude[j])
      temp$name <- c(temp$name, as.character(dat$name[j]))
      temp$width <- c(temp$width, w)

    } else {next(i)}

  }

  temp$colorindex <- temp$width * 20
  m <- max(temp$colorindex) * 1.5
  # Color_Assets <- colorFactor(rainbowPalette(m), levels = seq(1,m))
  Color_Assets <- colorFactor(heat.colors(m), levels = seq(1,m))



  pal <- colorNumeric(
    palette = heat.colors(m),
    domain = seq(1,m))

  pal <- colorNumeric(
    palette = heat.colors(10),
    domain = seq(1,m))


  qpal <- colorQuantile("RdYlBu", temp$colorindex, n = 5)

#   addCircleMarkers(map, data = temp, lng = ~ long, lat = ~ lat,
#                    popup = paste(as.character(temp$name)), radius = ~ width,
#                    color = ~Color_Assets(temp$colorindex), stroke = FALSE, fillOpacity = 0.5)  %>%

  addCircleMarkers(map, data = temp, lng = ~ long, lat = ~ lat,
                   popup = paste(as.character(temp$name), temp$colorindex), radius = 5,
                   color = ~qpal(temp$colorindex), stroke = FALSE, fillOpacity = 0.5)  %>%
  addLegend(map, position = "bottomright", pal = qpal, values = temp$colorindex,
            title = "Length of flood record",
            opacity = 1)


}





##################################################

# OTHER OPTIONS WITH LEAFLET
#
# # addCircleMarkers(-122.42, 37.78, popup = 'Bay Area', radius = 5, color = 'red')
#
# # For plotting many markers according to a dataset
# norway_map %>%
#   addTiles() %>%
#   setView(-122.42, 37.78, zoom = 13) %>%
#   addMarkers(data = data, lng = ~ X, lat = ~ Y, popup = data$Category)
# norway_map
#
# # For clustering
# norway_map %>%
#   addTiles() %>%
#   setView(-122.42, 37.78, zoom = 13) %>%
#   addCircleMarkers(data = data, lng = ~ X, lat = ~ Y, radius = 5,
#                    color = ~ ifelse(Category == 'BRIBERY', 'red', 'blue'),
#                    clusterOptions = markerClusterOptions())


