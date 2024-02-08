library(R.matlab)
library(rgdal)
library(sf)
######################## DATA CLEANING FILE ####################################
## Use raster to download USA and Canada raw shapefiles and combine and reduce size
## NOTE: This block isextremely slow because of large initial shapefiles
canadian_provinces <- raster::getData(country = 'Canada', level = 1)
canadian_provinces <- readRDS('gadm36_CAN_1_sp.rds')
rgdal::writeOGR(canadian_provinces, '.', 'canada', driver = 'ESRI Shapefile')
##Use mapshaper.org to reduce shapefile size, then read in
canadian_provinces <- st_read('canada.shp')


## Assumes that the current working directory is community-detect-raw-data
greatcircle <- function(p1, p2){
  ## Convert to radians
  p1 <- p1 * pi / 180
  p2 <- p2 * pi / 180
  return(6378 * acos(
    sin(p1[, 1]) * sin(p2[, 1]) + cos(p1[, 1]) * cos(p2[, 1]) * cos(p1[, 2] - p2[, 2])
  ))
}
## AIRPORT REACHABILITY NETWORK
conns <- read.delim('reachability.txt')[-(1:5), ]
conns <- as.data.frame(matrix(as.numeric(unlist(strsplit(conns, ' '))), 
                              ncol = 3, byrow = TRUE))
n <- max(conns) + 1
links <- list(A = matrix(0, n, n), X = 0, dis = matrix(0, n, n))
links$A[(conns$V1 + 1) + conns$V2 * n] <- 1
sum(links$A)
## Make network undirected
links$A[upper.tri(links$A)] <- pmax(links$A[upper.tri(links$A)], t(links$A)[upper.tri(links$A)])
links$A[lower.tri(links$A)] <- t(links$A)[lower.tri(links$A)]
sum(links$A)
adj <- read.csv('reachability-meta.csv')

links$X <- (cbind(log(adj$metro_pop), adj$longitude, adj$latitude))
## Calculate great-circle distance between airports, then scale
for (i in 1:n){
  links$dis[i, ] <- greatcircle(adj[i, 4:5], adj[, 4:5])
}
diag(links$dis) <- 0
links$dis <- links$dis / sd(links$dis)
## Include city size in distance calculation
size_diff <- abs(outer(log(adj$metro_pop), log(adj$metro_pop), '-')) 
size_diff <- size_diff / sd(size_diff)
## Combine and scale
links$dis <- sqrt(links$dis^2 + size_diff^2)
links$dis <- links$dis / sd(links$dis)
saveRDS(links, 'airports.RData')

## FACEBOOK 
meta <- readMat('facebook-ego-0.10.mat')
links <- list(A = as.matrix(meta$A), X = 0, dis = 0)

## Use only top ten most used covariates for more simplicity
top_inds <- order(colMeans(meta$features), decreasing = TRUE)[1:10]
links$X <- scale(meta$features[, top_inds])
links$dis <- as.matrix(dist(links$X, diag = TRUE, upper = TRUE))
saveRDS(links, 'facebook.RData')



