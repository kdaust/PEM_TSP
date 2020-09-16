source_python("./mTSP.py")

datLocGit <- here("InputData/Boundary") ## Data
covLoc <- here("BoundaryCovs") ## Too big for git data

###read in roads
rdsAll <- st_read(dsn = paste0(datLocGit,"/DRA_edited_stage2_boundarytsa.gpkg"))
rdsAll <- rdsAll[,c("ROAD_SURFACE","ROAD_CLASS","trail")]
rdsAll <- rdsAll[rdsAll$ROAD_SURFACE != "overgrown",]
colnames(rdsAll)[2] <- "road_surface"
rdsAll$road_surface[!is.na(rdsAll$trail)] <- "trail"
rdsAll <- as.data.table(rdsAll) %>% st_as_sf()

##road speed
rSpd <- fread("road_spd_boundary.csv")
#rSpd[,Values := (1/speed)/40]
rSpd[,Values := speed] ##convert to m/h
rdsAll <- merge(rdsAll, rSpd, by = "road_surface", all = F)
rdsAll <- rdsAll[,"Values"]
allRast <- raster(rdsAll, resolution = 25)

rdsAll <- st_buffer(rdsAll,dist = 25, endCapStyle = "SQUARE", joinStyle = "MITRE")
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
tr1[adj] <- (3/5)*(6*exp(-3.5*abs(tr1[adj] + 0.08))) ##tobler's hiking function * 3/5 - gives km/h
tr1 <- tr1*1000 ##now roads are correct conductance (h/m), and walking in m/h
tr2 <- geoCorrection(tr1) ##have to geocorrect this part again
tr1[adj] <- tr2[adj] ##tr1 values are now all conductance in h/metre

# tr1[rdAdj] <- 1/0.02
# tr1[cbind(rdAdj[,2],rdAdj[,1])] <- 1/0.02
# tr1 <- geoCorrection(tr1)

cities <- st_read(paste0(datLocGit,"/Boundary_PEM_prep.gdb"), layer = "A2_BCGW_MAJOR_CITIES_500m")
start <- cities[cities$NAME == "Osoyoos","NAME"] %>% as("Spatial")

p2 <- st_read(dsn = "./Stage2_TSP/s2_MSdm1_IDFdm1_pts.gpkg")
p2$id <- seq_along(p2$id)
pnts <- p2
colnames(pnts)[1] <- c("name")
pnts <- pnts["name"]
colnames(pnts)[2] <- "geometry"
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

maxTime <- 10L ##hours
## time per transect
plotTime <- 60L ##mins
temp <- dMat2[-nrow(dMat2),-ncol(dMat2)]
maxDist <- sum(temp[upper.tri(temp)])
pen <- as.integer(rep(maxDist * 5, nrow(p2))) ##penalty

n = nrow(dMat2)-1
ndays <- as.integer(n/4)
indStart <- as.integer(rep(n,ndays))
##run vehicle routing problem from python script
## GCS is global span cost coefficient
vrp <- py_mTSP(dat = dMat2,num_days = ndays, start = indStart, end = indStart, 
               max_cost = maxTime*60L, plot_time = plotTime,
               penalty =  pen, arbDepot = F, GSC = 1L)
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
      temp3$Segment = i
      temp3
    }
    out$DropSite = j
    out
  }
  
}

paths <- st_transform(paths, 3005)
st_write(paths, dsn = "BoundaryTSP_Stage2.gpkg", layer = "Paths", append = T, driver = "GPKG")  
st_write(pnts, dsn = "BoundaryTSP_Stage2.gpkg", layer = "Points", append = T)
