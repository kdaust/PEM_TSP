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
covLoc <- here("Covariates")
### landscape levels covariates
covars <- paste(covLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
                          "25m_MRVBF_Classified_IS64Low6Up2.tif","DEM_25m.tif"), sep = "/")# 
ancDat <- raster::stack(covars)
proj4string(ancDat) <- "+init=epsg:3005"
##in this case we're only using walkFast
rd1 <- 0.0003125
rd2 <- 0.000625
track <- 0.00125
walkFast <- 0.01667# 0.0125
walkSlow <- 0.01667
slopeAdjust <- function(slope){1+((slope-25)*0.02)}

# read in slope data

slope_raster <-  grep("^slope", list.files(covLoc))
slope <- raster(list.files(covLoc, full.name = TRUE)[slope_raster])

proj4string(slope) <- "+init=epsg:3005"

alt <- raster(paste0(datLoc, "/dem.tif"))
proj4string(alt) <- "+init=epsg:3005"
# read in already sampled locations
included <- st_read(paste0(datLocGit,"/ESSF_sampledpairs.gpkg"))

##read in buffer
buff <- st_read(paste0(datLoc,"/ESSF_Buffer.gpkg"))

## clip to just ESSF
boundary <- st_read(paste0(datLoc,"/bec_edited.gpkg"))
boundary <- boundary[,"MAP_LABEL"]
boundary <- boundary[boundary$MAP_LABEL == "ESSFmc",] ## set as mask for individual BGC
b2 <- st_union(boundary)
b2 <- fasterize(boundary, slope)
buff <- fasterize(buff, slope)
buff <- mask(buff,b2)
ancDat <- mask(ancDat,buff)

##read in drop points
heliDrop <- st_read(paste0(datLoc,"/DropCombined.gpkg"))
heliDrop <- heliDrop[,"name"]
heliDrop <- st_transform(heliDrop, 3005)
#heliDrop <- st_zm(heliDrop)
heliDrop <- heliDrop[boundary,]
#heliDrop <- st_transform(heliDrop, st_crs(pnts))
start <- as(heliDrop, "Spatial")

#slope <- raster(paste0(shapes_path,"/slope.tif"))

# cost <- tan(slope)*100
# cost[cost < 25] <- walkFast 
# cost[cost >= 25] <- walkFast*slopeAdjust(cost[cost >= 25])
# names(cost) <- "cost"

altDiff <- function(x){x[1]-x[2]}
tr <- transition(alt, transitionFunction = altDiff, directions = 8, symm = F) 
tr <- geoCorrection(tr)
adj <- adjacent(alt, cells = 1:ncell(alt), pairs = T, directions = 8)
tr[adj] <- 1.5*exp(-3.5*abs(tr[adj] + 0.1)) ##tobler's hiking function
tr <- geoCorrection(tr)

# read in already sampled locations
included <- st_read(paste0(datLocGit,"/ESSF_sampledpairs.gpkg"))
included <- st_transform(included, st_crs(ancDat))

source_python("./mTSP.py")
maxTime <- 8L ##hours
plotTime <- 45L ##mins

listComb <- function(a,b){
  t1 <- rbind(a$sampled_data,b$sampled_data)
  t2 <- b$obj
  t3 <- c(a$final_obj,b$final_obj)
  list(sampled_data = t1, obj = t2, final_obj = t3)
}

createLayout <- function(startPnts, toInclude, nPoints, iter = 1000){
  spoints <- foreach(x = 1:2, .combine = listComb) %do% {
    acost <- accCost(tr, startPnts[x,])
    acost <- acost/3600
    lays <- stack(ancDat,acost)
    names(lays) <- c("DAH","LFC","MRVBF","DEM","cost")
    
    incPnts <- raster::extract(lays, toInclude, sp = T)
    incPnts <- incPnts[,-(1:3)]
    incPnts <- st_as_sf(incPnts)
    if(x == 2){
      inc2 <- templhs$sampled_data
      temp <- raster::extract(acost, inc2["geometry"], sp = T) %>% st_as_sf()
      inc2$cost <- temp$layer
      incPnts <- rbind(incPnts, inc2)
    }
    
    s <- sampleRegular(lays , size = 500000, sp = TRUE) # sample raster
    s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
    s2 <- st_as_sf(s)
    ## have to add already sampled points to data
    s <- rbind(incPnts, s2)
    
    ### get sample locations
    templhs <- clhs_dist(s,size = nPoints/2+nrow(incPnts), minDist = 260, maxCost = 2, include = 1:nrow(incPnts), 
                         cost = "cost", iter = iter, simple = F,progress = T)
    
    list(sampled_data = templhs$sampled_data,obj = templhs$obj, final_obj = templhs$final_obj)
  }

  pnts <- spoints$sampled_data
  
  plot(acost)
  plot(pnts, add = T)
  
  p2 <- st_as_sf(pnts)
  pnts <- pnts[,"DAH"]
  colnames(pnts) <- c("name","geometry")
  st_geometry(pnts) <- "geometry"
  startPnts <- st_as_sf(startPnts) %>% st_transform(st_crs(pnts))
  pnts <- rbind(pnts, startPnts)
  pnts2 <- as(pnts, "Spatial")
  
  ## create distance matrix between sample points
  test <- costDistance(tr,pnts2,pnts2)
  dMat2 <- as.matrix(test)
  dMat2 <- dMat2/60
  dMat2[is.infinite(dMat2)] <- 1000
  
  
  ## note that indexing in python starts at 0, not 1
  ## to not allow dropped sites, set penalty > 10000
  
  ##penalty
  objVals <- spoints[["final_obj"]][1:nPoints]
  objVals <- max(objVals) - objVals
  minPen <- (maxTime*60L)/2L
  maxPen <- (maxTime*60L)*2L
  objVals <- scales::rescale(objVals, to = c(minPen,maxPen))
  objVals <- as.integer(objVals)

  # dMat2[51:60,] <- 0
  # dMat2[,51:60] <- 0
  ##fixed start, arbitrary end
  dMat2 <- rbind(dMat2, rep(0, ncol(dMat2)))
  dMat2 <- cbind(dMat2, rep(0, nrow(dMat2)))
  indStart <- as.integer(nPoints:(nPoints+nrow(startPnts)-1))
  indEnd <- as.integer(rep(nrow(dMat2)-1,nrow(startPnts)))
  vrp <- py_mTSP(dat = dMat2,num_days = 2L, start = indStart, end = indEnd, 
                 max_cost = maxTime*60L, plot_time = plotTime, penalty =  objVals, arbDepot = T)
  
  return(list(route = vrp, objective = spoints$obj[iter],pnts = pnts))
}

worker.init <- function(){
  Rcpp::sourceCpp("CppCLHS.cpp")
}

require(doParallel)
cl <- makePSOCKcluster(detectCores()-2)
clusterCall(cl, worker.init)
registerDoParallel(cl)

dropInd <- 1:nrow(heliDrop)
combs <- combn(dropInd, m = 2)

outStats <- foreach(i = 1:ncol(combs), .combine = c, 
                    .packages = c("Rcpp","reticulate","sf","raster","gdistance","foreach"),
                    .noexport = "c_cor") %dopar% {
  source_python("./mTSP.py")
  res <- createLayout(startPnts = start[c(combs[1,i],combs[2,i]),],toInclude = included, nPoints = 10)
  route <- res$route
  temp <- data.frame(start = paste(combs[1,i],combs[2,i],sep = "_"),cost = sum(unlist(route[[2]])), 
             num = paste(length(route[[1]][["0"]]),length(route[[1]][["1"]])), objFun = res$objective,
             totNum = sum(length(route[[1]][["0"]]),length(route[[1]][["1"]])))
  out <- list(list(stats = temp, solution = route,points = res$pnts))
  names(out) <- i
  out
}

stats <- foreach(i = 1:ncol(combs), .combine = rbind) %do% {
  outStats[[i]][["stats"]]
}
stats <- stats[stats$num == "6 6",]
ids <- rownames(stats[order(stats$cost),])

for(x in 1:length(ids)){
  writeLayout(id = ids[x], filename = paste0("ESSFtest2_",x,".gpkg"))
}

writeLayout <- function(id,filename){
  vrp <- outStats[[id]][["solution"]]
  pnts <- outStats[[id]][["points"]]
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

  paths <- st_transform(paths, 3005)
  st_write(paths, dsn = filename, layer = "Paths", append = T, driver = "GPKG")  
  
  ## label points
  p2 <- pnts
  p2$PID <- seq_along(p2$name)
  p2 <- p2[,"PID"]
  p2$DropLoc <- NA
  p2$Order <- NA
  for(i in 0:(length(result)-1)){
    p1 <- result[[as.character(i)]]+1
    p1 <- p1[-c(1,length(p1))]
    p2$DropLoc[p1] <- i
    p2$Order[p1] <- 1:(length(p1))
  }
  p2 <- st_transform(p2, 3005)
  st_write(p2, dsn = filename,layer = "Points", append = T,overwrite = T, driver = "GPKG")
  return(TRUE)
}

# ## this one for two routes at each drop
# vrp <- py_mTSP(dat = dMat2,num_days = 20L, start = c(50:59,50:59), 
#                end = c(50:59,50:59), max_cost = maxTime*60L, plot_time = plotTime,penalty =  maxTime*60L+5L)
# ### unsliced (cost layer includes all drop points)


# temp <- mask(acost, buff)
# writeRaster(temp, "CostSurface.tif",format = "GTiff")

#################################################################

# ## sliced clhs
# getSample <- function(index){
#   acost <- accCost(tr, start[index,])
#   
#   tempBuff <- st_buffer(heliDrop[index,],dist = 4000)
#   tempBuffR <- fasterize(tempBuff, acost)
#   acost <- mask(acost, tempBuffR, updatevalue = 10000)
#   lays <- stack(ancDat,acost)
#   names(lays) <- c("DAH","LFC","MRVBF","DEM","cost")
#   
#   s <- sampleRegular(lays , size = 500000, sp = TRUE) # sample raster
#   s <- s[!is.na(s$DAH) & !is.infinite(s$cost),]
#   return(s)
# }
# 
# s <- getSample(1)
# spoints <- clhs(s,size = 5, cost = "cost", iter = 5000, simple = F,progress = T)
# for(site in 2:10) {
#   prevSampled <- spoints$sampled_data
#   s <- getSample(site)
#   s <- rbind(prevSampled,s)
#   spoints <- clhs(s,size = length(prevSampled)+5,include = 1:length(prevSampled), 
#                   cost = "cost", iter = 5000, simple = F,progress = T)
# }
# ##now go to line 74