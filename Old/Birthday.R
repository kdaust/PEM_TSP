library(sf)
paths <- st_read(dsn = "C:/Users/kirid/Downloads/x_malkow/x_malkow/gisData/vec/Trails", layer = "Trail_paths")
paths <- st_transform(paths, 4326)
points <- st_read(dsn = "C:/Users/kirid/Downloads/x_malkow/x_malkow/gisData/vec/Trails", layer = "Trail_Junctions")
points <- st_transform(points, 4326)
c <- st_coordinates(points)
library(stplanr)
library(dodgr)
paths$type <- 1
graph <- weight_streetnet(paths, type_col = "type", id_col = "id", wt_profile = 1)


paths <- st_buffer(paths,dist = 2)
library(fasterize)
r <- raster("C:/Users/kirid/Downloads/x_malkow/x_malkow/gisData/tif/RasterTrails.tif")
r <- fasterize(paths,r, background = 10000)
plot(r)

r[r == 1] <- 0.004
plot(r)
library(gdistance)
tr <- transition(r,transitionFunction = function(x) 1/mean(x), directions = 8)

points2 <- points[-c(1,21,24),]
rownames(points2) <- NULL
c <- st_coordinates(points2)
d <- dodgr_dists (graph, from = c, to = c)

library(reticulate)
source_python("./_functions/mTSP.py")
##don't forget python starts indexing at 0 not 1
#vrp <- py_mTSP2(dat = d,num_days = 10L, start_node = 20L, max_cost = 8L, plot_time = 1L)
vrp <- py_mTSP(dat = d,num_days = 6L, start_node = 20L, max_cost = 2000L, plot_time = 250)
totDistTSP <- sum(unlist(vrp[[2]]))
vrp <- vrp[[1]]
library(foreach)


points2$Day <- NA
points2$Order <- NA
for(i in 0:5){
  p1 <- vrp[[as.character(i)]]+1
  points2$Day[p1] <- i
  points2$Order[p1[-1]] <- 1:(length(p1)-1)
}
st_write(points2, dsn = "VRP",layer = "Plots", append = T,overwrite = T, driver = "ESRI Shapefile")

daddy <- list(d1 = c("C","B","J"),d2= c("S","M"),d3 = c("D","E","F","G"),
              d4 = c("K","L","N","P"),d5 = c("H","I","Shed"),d6 = c("Q","R","nearQ","O"))

paths <- foreach(j = 0:4, .combine = rbind) %do% {
  p1 <- vrp[[as.character(j)]]+1
  p2 <- points2[p1,]
  
  out <- foreach(i = 1:(length(p1)-1), .combine = rbind) %do% {
    temp1 <- c[p1[i],]
    temp2 <- c[p1[i+1],]
    temp3 <- dodgr_paths(graph, from = temp1, to = temp2)
    verts <- dodgr_vertices(graph)
    v1 <- verts[match(temp3[[1]][[1]], verts$id),]
    v2 <- st_multipoint(x = as.matrix(v1[,c("x","y")]), dim = "XY") %>% 
      st_cast("LINESTRING") %>% 
      st_sfc()
    temp <- st_sf(v2, ID = i)
    temp
  }
  out$Day = j
  out
}

st_crs(paths) <- st_crs(points)
st_write(paths, dsn = "VRP",layer = "Test1",driver = "ESRI Shapefile")


names <- c("rubus","rosa","lathyrus","vicia","trifolium","castillega","viola","ranunculus","arnica",
           "aquillegia","hieracium","achillea","bellis","galium","taraxacum")
n1 <- strsplit(names,split = "")
out <- data.table()
for(i in 1:length(n1)){
  t1 <- n1[[i]] %>% as.data.table() %>% t()
  out <- rbind(out,t1, fill = T)
}

