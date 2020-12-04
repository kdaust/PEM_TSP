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
covars <- paste(covLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
                          "25m_MRVBF_Classified_IS64Low6Up2.tif","dem.tif"), sep = "/")# ,"DEM_25m.tif"
layerNames <- c("DAH","LFC","MRVBF","DEM","cost") ##need to change this if you change the layers
ancDat <- raster::stack(covars)
proj4string(ancDat) <- "+init=epsg:3005"


## dem for transtion layer
alt <- raster(paste0(covLoc, "/dem.tif"))
proj4string(alt) <- "+init=epsg:3005"

allRast <- raster(paste0(covLoc,"/Road_Rast_Small.tif"))    
allRast[allRast == 255] <- NA
allRast <- trim(allRast)

###read in roads
rdsAll <- st_read(paste0(datLocGit,"/road_access_for_cost.gpkg"))
rdsAll <- rdsAll[,"DESCRIPTIO"]
colnames(rdsAll)[1] <- "road_surface"
rdsAll <- as.data.table(rdsAll) %>% st_as_sf()

##road speed
rSpd <- fread("road_speed.csv")
#rSpd[,Values := (1/speed)/40]
rSpd[,Values := speed] ##convert to m/h
rdsAll <- merge(rdsAll, rSpd, by = "road_surface", all = F)
rdsAll <- rdsAll[,"Values"]
rdsAll <- st_buffer(rdsAll,dist = 35, endCapStyle = "SQUARE", joinStyle = "MITRE")
st_write(rdsAll, "TestBufferRds.gpkg")
rdsAll <- st_cast(rdsAll, "MULTIPOLYGON")
rdsRast <- fasterize(rdsAll, allRast, field = "Values")

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

start <- st_read("SmithersStart.gpkg") %>% as("Spatial")
acost <- accCost(tr1,start)
plot(acost)
writeRaster(acost, "TestAcost.tif", "GTiff")

acost2 <- projectRaster(acost, ancDat)
acost2 <- mask(acost2, alt)
lays <- stack(ancDat,acost2)
names(lays) <- layerNames
s <- sampleRegular(lays , size = 5000000, sp = TRUE) # sample raster
s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
s <- st_as_sf(s)

templhs <- clhs_dist(s,size = 50, minDist = 260, 
                     maxCost = 1.7, include = NULL, 
                     cost = "cost", iter = 5000, simple = F,progress = T)

pnts <- templhs$sampled_data

plot(acost2)
plot(pnts, add = T)

p2 <- st_as_sf(pnts)
pnts <- pnts[,"DAH"]
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

##penalty based on quality of points
objVals <- templhs[["final_obj"]]
objVals <- max(objVals) - objVals

maxTime <- 8L ##hours
## time per transect
plotTime <- 45L ##mins
minPen <- (maxTime*60L)/2L
maxPen <- (maxTime*60L)*2L
objVals <- scales::rescale(objVals, to = c(minPen,maxPen))
objVals <- as.integer(objVals)

##fixed start, arbitrary end

indStart <- as.integer(rep(50,10))
##run vehicle routing problem from python script
## GCS is global span cost coefficient
vrp <- py_mTSP(dat = dMat2,num_days = 10L, start = indStart, end = indStart, 
               max_cost = maxTime*60L, plot_time = plotTime, penalty =  objVals, arbDepot = F, GSC = 8L)
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
st_write(paths, dsn = "RoadTSP.gpkg", layer = "Paths", append = T, driver = "GPKG")  

## label points
p2 <- pnts
p2$PID <- seq_along(p2$name)
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
st_write(p2, dsn = "RoadTSP.gpkg",layer = "Points", append = T,overwrite = T, driver = "GPKG")

### number of points

SLcov <- c("twi.tif","valley_depth_2.tif","tca2.tif","swi_area_mod.tif","cov.tif","tpi.tif")
covars <- paste(covLoc, SLcov, sep = "/")# ,"DEM_25m.tif"
layerNames <- c("twi","valley","tca2","swi","cov","tpi") ##need to change this if you change the layers
ancDat <- raster::stack(covars)
proj4string(ancDat) <- "+init=epsg:3005"

###theoretical # to fill cov space
sampledDat <- st_read(paste0(datLocGit,"/transect1_30m_pts_att.gpkg"))
sampledDat <- sampledDat["mapunit1"]
temp <- raster::extract(ancDat, sampledDat)
dat <- data.table(cbind(as.character(sampledDat$mapunit1),temp))
setnames(dat, old = "V1", new = "Unit")
dat <- dat[!is.na(Unit),]
dat2 <- dat[,lapply(.SD, function(x){sd(as.numeric(x), na.rm = T)}),by = .(Unit)]
dat <- dat[,lapply(.SD, function(x){mean(as.numeric(x), na.rm = T)}),by = .(Unit)]
dat <- dat[grep("SBSmc2",Unit),]
dat2 <- dat2[grep("SBSmc2",Unit),]
##############################33

nums <- round(seq(10, 800, length.out = 15))
nums <- rep(nums, each = 15)

sumKSTest <- function(full,small){
  out <- 0
  for(i in 1:ncol(full)){
    ks <- ks.test(full[,i],small[,i])
    out <- out+ks$statistic
  } 
  return(out)
}
  
fullSet <- sampleRegular(ancDat, size = 1000000,useGDAL = T)
X <- fullSet
X <- X[complete.cases(X),]
actSize <- nrow(X)

worker.init <- function(){
  Rcpp::sourceCpp("CppCLHS.cpp")
}

require(doParallel)
cl <- makePSOCKcluster(detectCores()-2)
clusterCall(cl, worker.init)
registerDoParallel(cl)

test_res <- foreach(n = nums, .combine = rbind, .noexport = c("c_cor","obj_fn"), 
                    .packages = c("Rcpp","LaplacesDemon","foreach")) %dopar% {
                      smallSet <- clhs_fast(X,size = n, iter = 1000, simple = F,progress = T)
                      smallSet <- as.matrix(smallSet$sampled_data)
                      res1 <- sumKSTest(X,smallSet)
                      res2  <- meanKLD(X,smallSet, nb = 20)
                      data.frame(Num = n, KLyx = res2[1], KLxy = res2[2],KS = res1)
}

boxplot(KS ~ Num, data = test_res)
boxplot(KLxy ~ Num, data = test_res)
boxplot(KLyx ~ Num, data = test_res)

temp <- as.data.table(test_res)
temp <- temp[,.(y = quantile(KS,0.5)), by = .(Num)]
temp[,y := (y-min(y))/(max(y)-min(y))]
plot(y ~ Num, data = temp)
spFun <- splinefun(y = temp$Num, x = temp$y)

x = temp$Num
y = temp$y
plot(x, y, xlab="sample number",ylab = "KL Stat")          # Initial plot of the data
start <- list(k = 1,b1 = 0.05,b0 = 0)
fit1 <- nls(y ~ k*exp(-b1*x) + b0, start = start)
lines(x, fitted(fit1), col="red")

xx<- seq(1, 800,1)
jj <- predict(fit1,list(x=xx))

normalized = (jj-min(jj))/(max(jj)-min(jj))###standardise
plot(xx,normalized, type = "l")
approx(x = normalized, y = xx, xout = 0.10)$y###get 90th quantile




# mat <- matrix(data = c(NA,3,6,NA,5,4,80,80,80,NA,5,5,NA,6,4),nrow = 5,byrow = T)
# r <- crop(alt, extent(alt, 1,5,1,3))
# values(r) <- mat
# 
# rdIdx <- which(values(r) < 5)
# slpIdx <- which(values(r) > 500)
# rdAdj <- adjacent(r, cells = rdIdx, pairs = T, directions = 8)
# adj <- adjacent(r, cells = slpIdx, pairs = T, directions = 8)
# adj <- adj[!adj[,1] %in% rdIdx ,]
# adj <- adj[!adj[,2] %in% rdIdx ,]
# 
# tr <- transition(r,function(x){min(x[1],x[2])},directions = 8, symm = F)
# tr1 <- geoCorrection(tr)
# tr1[adj] <- (1.5*exp(-3.5*abs(tr1[adj] + 0.1))) ##tobler's hiking function
# tr1[adj] <- 1/0.36
# tr1[rdAdj] <- 1/0.02
# tr1[cbind(rdAdj[,2],rdAdj[,1])] <- 1/0.02
# tr1 <- geoCorrection(tr1)
# acost <- accCost(tr1,fromCoords = c(0.2,0.5))
