###Heli point sample plan
## Kiri Daust, July 2020  
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

Rcpp::sourceCpp("CppCLHS.cpp")
source("FastCLHS_R.R")

datLoc <- here("InputData") 
### landscape levels covariates
covars <- paste(datLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
                          "25m_MRVBF_Classified_IS64Low6Up2.tif","dem.tif"), sep = "/")
ancDat <- raster::stack(covars)

##in this case we're only using walkFast
rd1 <- 0.0003125
rd2 <- 0.000625
track <- 0.00125
walkFast <- 0.0125
walkSlow <- 0.01667
slopeAdjust <- function(slope){1+((slope-25)*0.02)}

# read in slope data
slope_raster <-  grep("^slope", list.files(datLoc))
slope <- raster(list.files(datLoc, full.name = TRUE)[slope_raster])

# read in already sampled locations
included <- st_read(paste0(datLoc,"/ESSF_samples.gpkg"))

##read in buffer
buff <- st_read("InputData/ESSF_Buffer.gpkg")

## clip to just ESSF
# boundary <- st_read(paste0(datLoc,"/bec_edited.gpkg"))
# boundary <- boundary[,"MAP_LABEL"]
# boundary <- boundary[grep("ESSF",boundary$MAP_LABEL),]
# b2 <- st_union(boundary)
# b2 <- fasterize(boundary, slope)
buff <- fasterize(buff, slope)
ancDat <- mask(ancDat,buff)

##read in drop points
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

### unsliced (cost layer includes all drop points)
acost <- accCost(tr, start)
lays <- stack(ancDat,acost)
names(lays) <- c("DAH","LFC","MRVBF","DEM","cost")

s <- sampleRegular(lays , size = 500000, sp = TRUE) # sample raster
s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
s2 <- st_as_sf(s)
## have to add already sampled points to data
included <- st_transform(included, st_crs(s))
incPnts <- raster::extract(lays, included, sp = T)
incPnts <- incPnts[,-(1:3)]
incPnts <- st_as_sf(incPnts)
s <- rbind(incPnts, s2)

### get sample locations
spoints <- clhs_dist(s,size = 50+nrow(incPnts), minDist = 260, maxCost = 2, include = 1:nrow(incPnts), 
                  cost = "cost", iter = 5000, simple = F,progress = T)
############################################################################

pnts <- spoints$sampled_data

plot(acost)
plot(pnts, add = T)


ac2 <- acost
ac2[ac2 > 2.5] <- 20
plot(ac2)

p2 <- st_as_sf(pnts)
heliDrop <- st_zm(heliDrop)
heliDrop <- st_transform(heliDrop, st_crs(pnts))
pnts <- pnts[,"DAH"]
colnames(pnts) <- c("name","geom")
st_geometry(pnts) <- "geom"
pnts <- rbind(pnts, heliDrop)
pnts2 <- as(pnts, "Spatial")

## create distance matrix between sample points
test <- costDistance(tr,pnts2,pnts2)
dMat2 <- as.matrix(test)
dMat2 <- dMat2*60
dMat2[is.infinite(dMat2)] <- 1000

source_python("./mTSP.py")
maxTime <- 8L ##hours
plotTime <- 45L ##mins
## note that indexing in python starts at 0, not 1
## to not allow dropped sites, set penalty > 10000

##penalty
objVals <- spoints[["final_obj"]][1:50]
objVals <- objVals - max(objVals)
minPen <- (maxTime*60L)/2L
maxPen <- (maxTime*60L)*2L
objVals <- scales::rescale(objVals, to = c(minPen,maxPen))
objVals <- as.integer(objVals)

## this one for two routes at each drop
vrp <- py_mTSP(dat = dMat2,num_days = 20L, start = c(50:59,50:59), 
               end = c(50:59,50:59), max_cost = maxTime*60L, plot_time = plotTime,penalty =  maxTime*60L+5L)

## one route at each drop
# dMat2[51:60,] <- 0
# dMat2[,51:60] <- 0
vrp <- py_mTSP(dat = dMat2,num_days = 10L, start = 50:59, end = 50:59, 
               max_cost = maxTime*60L, plot_time = plotTime, penalty =  objVals, arbDepot = F)

result <- vrp[[1]]

## create spatial paths
paths <- foreach(j = 0:(length(result)-1), .combine = rbind) %do% {
  if(length(result[[as.character(j)]]) > 2){
    cat("Drop site",j,"...\n")
    p1 <- result[[as.character(j)]]+1
    out <- foreach(i = 1:(length(p1)-1), .combine = rbind) %do% {
      temp1 <- pnts[p1[i],]
      temp2 <- pnts[p1[i+1],]
      temp3 <- shortestPath(tr,st_coordinates(temp1),
                            st_coordinates(temp2),output = "SpatialLines") %>% st_as_sf()
      temp3$Segment = i
      temp3
    }
    out$DropSite = j
    out
  }
  
}
filename <- "Heli_ObjPenalty.gpkg"

st_write(paths, dsn = filename, layer = "Paths", append = T, driver = "GPKG")  

## label points
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
st_write(p2, dsn = filename,layer = "Points", append = T,overwrite = T, driver = "GPKG")
temp <- mask(acost, buff)
writeRaster(temp, "CostSurface.tif",format = "GTiff")
#################################################################

## sliced clhs
getSample <- function(index){
  acost <- accCost(tr, start[index,])
  
  tempBuff <- st_buffer(heliDrop[index,],dist = 4000)
  tempBuffR <- fasterize(tempBuff, acost)
  acost <- mask(acost, tempBuffR, updatevalue = 10000)
  lays <- stack(ancDat,acost)
  names(lays) <- c("DAH","LFC","MRVBF","DEM","cost")
  
  s <- sampleRegular(lays , size = 500000, sp = TRUE) # sample raster
  s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
  return(s)
}

s <- getSample(1)
spoints <- clhs(s,size = 5, cost = "cost", iter = 5000, simple = F,progress = T)
for(site in 2:10) {
  prevSampled <- spoints$sampled_data
  s <- getSample(site)
  s <- rbind(prevSampled,s)
  spoints <- clhs(s,size = length(prevSampled)+5,include = 1:length(prevSampled), 
                  cost = "cost", iter = 5000, simple = F,progress = T)
}
##now go to line 74