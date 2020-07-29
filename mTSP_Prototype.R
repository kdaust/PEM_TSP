##Kiri Daust, July 2020

      library(clhs)
      library(sf)
      library(raster)
      library(sp)
      library(gdistance)
      library(foreach)
      library(data.table)
      library(fasterize)
      library(reticulate)

load("BoundaryLayers.Rdata")
s <- sampleRegular(ancDat, size = 1000000, sp = TRUE) # sample raster
s <- s[!is.na(s$clhc_covars.1) & !is.na(s$layer),]

t1 <- clhs(s, size = 25, iter = 5000, simple = F, progress = T, cost = 'layer')
pts <- t1$sampled_data

### add starting point
temp <- st_as_sf(pts) %>% st_transform(st_crs(cost_start_points)) ##cost_start_points from 03a_Stage1_SampleDesign
colnames(temp)[6] <- "geom"
st_geometry(temp) <- "geom"
pts2 <- rbind(temp["geom"],cost_start_points[1,"geom"])
pts2$ID <- seq_along(pts2$geom)
pts <- as(pts2,"Spatial")

load("CorrectedCostTransition.Rdata")
# tr <- transition(cost, transitionFunction = function(x) 1/mean(x), directions = 8)
# tr <- geoCorrection(tr, type = "c")
test <- costDistance(tr,pts,pts)
dMat2 <- as.matrix(test)
source_python("./_functions/mTSP.py")
##don't forget python starts indexing at 0 not 1
vrp <- py_mTSP(dat = dMat2,num_days = 4L, start_node = 25L, max_cost = 200L)
vrp

acost <- accCost(tr, pts[26,])
plot(acost)
plot(pts, add = T)
# p1 <- vrp$`3`+1
# p2 <- pts[p1,]
# plot(p2, add = T, col = "purple")
# tpLines <- shortestPath(tr,p2,p2,output = "SpatialLines")
# plot(tpLines, add = T)

###loops through results on VRP and create lines for shortest paths
paths <- foreach(j = 0:3, .combine = rbind) %do% {
  p1 <- vrp[[as.character(j)]]+1
  p2 <- pts[p1,]
  
  out <- foreach(i = 1:length(p1)-1, .combine = rbind) %do% {
    temp1 <- pts[p1[i],]
    temp2 <- pts[p1[i+1],]
    temp3 <- shortestPath(tr,temp1,temp2,output = "SpatialLines") %>% st_as_sf()
    #temp3$Segment = i
    temp3
  }
  out$Day = j
  out
}

st_write(paths, dsn = "TSP_Working",layer = "TSPRoutes", driver = "ESRI Shapefile")
st_write(pts2, dsn = "TSP_Points",layer = "SampleLocations", driver = "ESRI Shapefile")
writeRaster(acost,"AccCost.tif", driver = "GTiff")

