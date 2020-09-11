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

source("FastCLHS_R.R")
source_python("./mTSP.py")

datLocGit <- here("InputData/Boundary") ## Data
covLoc <- here("BoundaryCovs") ## Too big for git data
### landscape levels covariates
# covars <- paste(covLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
#                           "25m_MRVBF_Classified_IS64Low6Up2.tif","dem.tif"), sep = "/")# ,"DEM_25m.tif"
covars <- paste(covLoc, c("DAH_25m.tif","LandForm_Seive10_25m.tif", 
                          "Boundary_25m_MRVBF64Up2Low6.tif","slope_25m.tif"), sep = "/")
layerNames <- c("DAH","LFC","MRVBF","SLP","cost") ##need to change this if you change the layers
ancDat <- raster::stack(covars)
proj4string(ancDat) <- "+init=epsg:3005"

###read in roads
rdsAll <- st_read(dsn = paste0(datLocGit,"/Boundary_PEM_prep.gdb"), layer = "A2_BCGW_ROAD_ATLAS_MPAR")
rdsAll <- rdsAll[,c("ROAD_SURFACE","ROAD_CLASS")]
rdsAll <- rdsAll[rdsAll$ROAD_SURFACE != "overgrown",]
colnames(rdsAll)[2] <- "road_surface"
rdsAll <- as.data.table(rdsAll) %>% st_as_sf()

##road speed
rSpd <- fread("road_spd_boundary.csv")
#rSpd[,Values := (1/speed)/40]
rSpd[,Values := speed] ##convert to m/h
rdsAll <- merge(rdsAll, rSpd, by = "road_surface", all = F)
rdsAll <- rdsAll[,"Values"]
allRast <- raster(rdsAll, resolution = 25)

rdsAll <- st_buffer(rdsAll,dist = 35, endCapStyle = "SQUARE", joinStyle = "MITRE")
rdsAll <- st_cast(rdsAll, "MULTIPOLYGON")
rdsRast <- fasterize(rdsAll, allRast, field = "Values")

## dem for transtion layer
alt <- raster(paste0(covLoc, "/dem.tif"))
proj4string(alt) <- "+init=epsg:3005"

alt2 <- projectRaster(alt,rdsRast,method = 'ngb')
altAll <- merge(rdsRast, alt2)


trFn <- function(x){
  if(x[1] > 500 & x[2] > 500){
    x[1]-x[2]
  }else{
    min(x[1],x[2])
  }
}

rdIdx <- which(values(altAll) < 100)
slpIdx <- which(values(altAll) > 500)
rdAdj <- adjacent(altAll, cells = rdIdx, pairs = T, directions = 8)
adj <- adjacent(altAll, cells = slpIdx, pairs = T, directions = 8)
adj <- adj[!adj[,1] %in% rdIdx ,]
adj <- adj[!adj[,2] %in% rdIdx ,]

tr <- transition(altAll,trFn,directions = 8, symm = F) ##altDiff and speed (km/h)

tr1 <- geoCorrection(tr) ##divided by 25 - slope and conductance (km/h/m)
tr1[adj] <- 0.25*(6*exp(-3.5*abs(tr1[adj] + 0.08))) ##tobler's hiking function * 3/5 - gives km/h
tr1 <- tr1*1000 ##now roads are correct conductance (h/m), and walking in m/h
tr2 <- geoCorrection(tr1) ##have to geocorrect this part again
tr1[adj] <- tr2[adj] ##tr1 values are now all conductance in h/metre

# tr1[rdAdj] <- 1/0.02
# tr1[cbind(rdAdj[,2],rdAdj[,1])] <- 1/0.02
# tr1 <- geoCorrection(tr1)

cities <- st_read(paste0(datLocGit,"/Boundary_PEM_prep.gdb"), layer = "A2_BCGW_MAJOR_CITIES_500m")
start <- cities[cities$NAME == "Osoyoos","NAME"] %>% as("Spatial")
acost <- accCost(tr1,start)
plot(acost)
writeRaster(acost, "TestAcost.tif", "GTiff",overwrite = T)

acost2 <- projectRaster(acost, ancDat)
lays <- stack(ancDat,acost2)
names(lays) <- layerNames
bgc <- st_read(dsn = here("InputData/Boundary"), layer = "Boundary_BGC_dissolved")
bgc <- bgc[bgc$MAP_LABEL == "ESSFdc2",]
bgc <- st_cast(bgc,"POLYGON")
bgc$ID <- seq_along(bgc$MAP_LABEL)
bgc <- bgc[bgc$ID %in% c(1,2),]
bgc <- st_buffer(bgc, dist = -144)
maskRast <- fasterize(bgc,ancDat[[1]])
lays <- mask(lays,maskRast)
lays <- crop(lays,bgc)
cutblks <- st_read(here("InputData/Boundary/cutblock1.gpkg"))
cutblks <- st_crop(cutblks,bgc)
cutblks <- st_buffer(cutblks,144)
cutblkRast <- fasterize(cutblks, lays$DAH)
lays <- mask(lays, cutblkRast, inverse = T)
rdsAll <- st_read(dsn = paste0(datLocGit,"/Boundary_PEM_prep.gdb"), layer = "A2_BCGW_ROAD_ATLAS_MPAR")
rdsAll <- rdsAll[,c("ROAD_SURFACE","ROAD_CLASS")]
rdsAll <- rdsAll[rdsAll$ROAD_SURFACE != "overgrown",]
colnames(rdsAll)[2] <- "road_surface"
rds <- st_crop(rdsAll, bgc)
rds <- st_buffer(rds, dist = 100)
lays <- mask(lays, rds, inverse = T)

s <- sampleRegular(lays , size = 5000000, sp = TRUE) # sample raster
s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
s <- st_as_sf(s)
temp <- st_drop_geometry(s) %>% as.matrix()
templhs <- c_clhs(temp,size = 10, include = NULL, 
                     i_cost = 5, iter = 20000)

idx <- templhs$indeces
pnts <- s[idx,]
plot(acost2)
plot(pnts, add = T)
spPnts <- as(pnts,"Spatial")
origPnts <- pnts

create_acost <- function(altAll){
  rdIdx <- which(values(altAll) < 100)
  slpIdx <- which(values(altAll) > 500)
  rdAdj <- adjacent(altAll, cells = rdIdx, pairs = T, directions = 8)
  adj <- adjacent(altAll, cells = slpIdx, pairs = T, directions = 8)
  adj <- adj[!adj[,1] %in% rdIdx ,]
  adj <- adj[!adj[,2] %in% rdIdx ,]
  
  tr <- transition(altAll,trFn,directions = 8, symm = F) ##altDiff and speed (km/h)
  
  tr1 <- geoCorrection(tr) ##divided by 25 - slope and conductance (km/h/m)
  tr1[adj] <- 0.25*(6*exp(-3.5*abs(tr1[adj] + 0.08))) ##tobler's hiking function * 3/5 - gives km/h
  tr1 <- tr1*1000 ##now roads are correct conductance (h/m), and walking in m/h
  tr2 <- geoCorrection(tr1) ##have to geocorrect this part again
  tr1[adj] <- tr2[adj] ##tr1 values are now all conductance in h/metre
  return(tr1)
}

pntsTSP <- foreach(lhsPnt = 1:nrow(pnts), .combine = rbind) %do% {
  cat("Processing point",lhsPnt,"\n")
  pnt <- pnts[lhsPnt,]
  pntBuff <- st_buffer(pnt, dist = 1000)
  altSmall <- mask(altAll,pntBuff) %>% crop(pntBuff)
  trTemp <- create_acost(altAll = altSmall)
  acostTemp <- accCost(trTemp,as(pnt,"Spatial"))
  rdsTemp <- rdsAll[pntBuff,]
  rdsTemp$RdNum <- NA
  ints <- st_intersects(rdsTemp,sparse = F)
  
  ### figure out connected roads
  flag1 <- T
  count <- 0
  idxAll <- 1:ncol(ints)
  idxOpts <- idxAll
  while(flag1){
    idx <- which(ints[,idxOpts[1]])
    flag2 <- T
    if(length(idx) > 1){
      while(flag2){
        temp <- which(apply(ints[,idx], 1, FUN = function(x){any(x)}))
        if(setequal(idx,temp)){
          flag2 <- F
        }
        idx <- union(idx,temp)
      }
    }
    rdsTemp$RdNum[idx] <- count
    count <- count+1
    idxOpts <- setdiff(idxOpts,idx)
    if(length(idxOpts) == 0){
      flag1 <- F
    }
  }
  
  plot(acostTemp)
  plot(rdsTemp, add = T)
  rdPnts <- foreach(rdNum = unique(rdsTemp$RdNum), .combine = rbind) %do% {
    rds <- rdsTemp[rdsTemp$RdNum == rdNum,]
    rds <- st_buffer(rds, dist = 35)
    acost2 <- mask(acostTemp, rds)
    acost2 <- trim(acost2)
    minCost <- min(values(acost2), na.rm = T)
    cPnts <- rasterToPoints(acost2, fun = function(x){x == minCost}, spatial = T)
    cPnts
  }
  rdPnts <- st_as_sf(rdPnts)
  rdPnts$LHSpnt = lhsPnt
  rdPnts
}



p2 <- pntsTSP
pnts <- pntsTSP[,"LHSpnt"]
colnames(pnts) <- c("name","geometry")
st_geometry(pnts) <- "geometry"
start <- st_as_sfc(start) %>% st_transform(st_crs(pnts))
startPnts <- st_as_sf(data.frame(name = "Start",geometry = start))
pnts <- rbind(pnts, startPnts)
pnts2 <- as(pnts, "Spatial")

## create distance matrix between sample points
test <- costDistance(tr1,pnts2,pnts2)
dMat2 <- as.matrix(test)
dMat2 <- dMat2*60
dMat2[is.infinite(dMat2)] <- 1000

addCost <- st_drop_geometry(pntsTSP)
addCost$layer <- addCost$layer*60 ##covert to minutes

nPoints <- nrow(dMat2)
for(i in 1:nPoints){
  for(j in 1:nPoints){
    if(i != j){
      if(i == nPoints | j == nPoints){
        dMat2[i,j] <- dMat2[i,j] + addCost$layer[min(i,j)]
      }else{
        dMat2[i,j] <- dMat2[i,j] + addCost$layer[i] + addCost$layer[j]
      }
    }
  }
}

dupPoints <- list()
for(ix in unique(pntsTSP$LHSpnt)){
  temp <- rownames(pntsTSP[pntsTSP$LHSpnt == ix,]) %>% as.numeric()
  temp <- as.integer(temp - 1)
  dupPoints[[ix]] <- temp
}

##penalty based on quality of points
objVals <- templhs[["final_obj"]]
objVals <- max(objVals) - objVals

maxTime <- 11L ##hours
## time per transect
plotTime <- 40L ##mins
temp <- dMat2[1:10,1:10]
maxDist <- sum(temp[upper.tri(temp)])
minPen <- maxDist * 2
maxPen <- maxDist * 5
objVals <- scales::rescale(objVals, to = c(minPen,maxPen))
objVals <- as.integer(objVals)

n = nrow(dMat2)-1
ndays <- 2L
indStart <- as.integer(rep(n,ndays))
##run vehicle routing problem from python script
## GCS is global span cost coefficient
vrp <- py_mTSP(dat = dMat2,num_days = ndays, start = indStart, end = indStart, 
               max_cost = maxTime*60L, plot_time = plotTime, duplicates = dupPoints,
               penalty =  objVals, arbDepot = F, GSC = 8L)
result <- vrp[[1]]

## create spatial paths
paths <- foreach(j = 0:(length(result)-1), .combine = rbind) %do% {
  if(length(result[[as.character(j)]]) > 2){
    cat("Drop site",j,"...\n")
    p1 <- result[[as.character(j)]]+1
    out <- foreach(i = 1:(length(p1)-1), .combine = rbind) %do% {
      temp1 <- pnts[p1[i],]
      temp2 <- pnts[p1[i+1],]
      temp3 <- shortestPath(tr1,st_coordinates(temp1),
                            st_coordinates(temp2),output = "SpatialLines") %>% st_as_sf()
      if(p1[i] != max(p1)){
        temp4 <- origPnts[pntsTSP$LHSpnt[p1[i]],]
        temp5 <- shortestPath(tr1,st_coordinates(temp1),
                              st_coordinates(temp4),output = "SpatialLines") %>% st_as_sf()
        temp3 <- rbind(temp3,temp5)
      }
      temp3$Segment = i
      temp3
    }
    out$DropSite = j
    out
  }
  
}

paths <- st_transform(paths, 3005)
st_write(paths, dsn = "BoundaryTSP_Fix.gpkg", layer = "Paths", append = T, driver = "GPKG")  

## label points
idx <- templhs$indeces
pnts <- s[idx,]
p2 <- pnts
p2$PID <- seq_along(p2$DAH)
p2 <- p2[,"PID"]
p2$DropLoc <- NA
p2$Order <- NA
for(i in 0:(length(result)-1)){
  p1 <- result[[as.character(i)]]+1
  p1 <- p1[-1]
  p2$DropLoc[p1] <- i
  p2$Order[p1] <- 1:length(p1)
}
p2 <- st_transform(p2, 3005)
st_write(p2, dsn = "BoundaryTSP_Fix.gpkg",layer = "Points", append = T,overwrite = T, driver = "GPKG")
