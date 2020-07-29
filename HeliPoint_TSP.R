###Heli point sample plan
library(clhs)
library(sf)
library(raster)
library(sp)
library(gdistance)
library(foreach)
library(data.table)
library(fasterize)
library(reticulate)
library(here)

datLoc <- here("InputData")
covars <- paste(datLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
                          "25m_MRVBF_Classified_IS64Low6Up2.tif","dem.tif"), sep = "/")
ancDat <- raster::stack(covars)

rd1 <- 0.0003125
rd2 <- 0.000625
track <- 0.00125
walkFast <- 0.0125
walkSlow <- 0.01667
slopeAdjust <- function(slope){1+((slope-25)*0.02)}

# read in slope data
slope_raster <-  grep("^slope", list.files(datLoc))
slope <- raster(list.files(datLoc, full.name = TRUE)[slope_raster])

boundary <- st_read(paste0(datLoc,"/bec_edited.gpkg"))
boundary <- boundary[,"MAP_LABEL"]
boundary <- boundary[grep("ESSF",boundary$MAP_LABEL),]
b2 <- st_union(boundary)
b2 <- fasterize(boundary, slope)
slope <- mask(slope, b2)
ancDat <- mask(ancDat,b2)

heliDrop <- st_read(paste0(datLoc,"/Deception_Heli_Samples.gpkg"))
heliDrop <- heliDrop[,"name"]
heliDrop <- st_transform(heliDrop, 3005)
start <- as(heliDrop, "Spatial")

#slope <- raster(paste0(shapes_path,"/slope.tif"))
cost <- tan(slope)*100
cost[cost < 25] <- walkFast 
cost[cost >= 25] <- walkFast*slopeAdjust(cost[cost >= 25])
names(cost) <- "cost"

tr <- transition(cost, transitionFunction = function(x) 1/mean(x), directions = 8) 

acost <- accCost(tr, start)
plot(acost)

lays <- stack(ancDat,acost)
names(lays) <- c("DAH","LFC","MRVBF","DEM","cost")

s <- sampleRegular(lays , size = 500000, sp = TRUE) # sample raster
s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
#s@data <- cbind(s@data,s@coords)

spoints <- clhs(s,size = 50, cost = "cost", iter = 5000, simple = F,progress = T)
pnts <- spoints$sampled_data
plot(pnts, add = T)
p2 <- st_as_sf(pnts)
heliDrop <- st_zm(heliDrop)
heliDrop <- st_transform(heliDrop, st_crs(pnts))
dropPnts <- as(heliDrop,"Spatial")
pnts <- pnts[,"DAH"]
colnames(pnts@data) <- "name"
pnts <- rbind(pnts, dropPnts)
pnts2 <- st_as_sf(pnts)

test <- costDistance(tr,pnts,pnts)
dMat2 <- as.matrix(test)
dMat2 <- dMat2*60
dMat2[is.infinite(dMat2)] <- 1000

source_python("./mTSP.py")
vrp <- py_mTSP(dat = dMat2,num_days = 20L, start = c(50:59,50:59), end = c(50:59,50:59), max_cost = 8L*60L, plot_time = 45L)
vrp <- py_mTSP(dat = dMat2,num_days = 10L, start = 50:59, end = 50:59, max_cost = 9L*60L, plot_time = 45L, penalty = 500L)

result <- vrp[[1]]

paths <- foreach(j = 0:(length(result)-1), .combine = rbind) %do% {
  if(length(result[[as.character(j)]]) > 2){
    cat("Drop site",j,"...\n")
    p1 <- result[[as.character(j)]]+1
    out <- foreach(i = 1:(length(p1)-1), .combine = rbind) %do% {
      temp1 <- pnts2[p1[i],]
      temp2 <- pnts2[p1[i+1],]
      temp3 <- shortestPath(tr,st_coordinates(temp1),
                            st_coordinates(temp2),output = "SpatialLines") %>% st_as_sf()
      temp3$Segment = i
      temp3
    }
    out$DropSite = j
    out
  }
  
}

st_write(paths, dsn = "Heli_DropAllowed.gpkg", layer = "Paths", append = T, driver = "GPKG")  

p2$PID <- seq_along(p2$DAH)
p2 <- p2[,"PID"]
p2$DropLoc <- NA
p2$Order <- NA
for(i in 0:(length(result)-1)){
  p1 <- result[[as.character(i)]]+1
  p1 <- p1[-c(1,length(p1))]
  p2$DropLoc[p1] <- i
  p2$Order[p1] <- 1:(length(p1))
}
st_write(p2, dsn = "Heli_DropAllowed.gpkg",layer = "Points", append = T,overwrite = T, driver = "GPKG")
writeRaster(acost, "CostSurface.tif",format = "GTiff")
