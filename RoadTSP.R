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
source_python("./mTSP.py")
  
datLocGit <- here("InputData") ## Data
covLoc <- here("Covariates") ## Too big for git data
### landscape levels covariates
covars <- paste(covLoc, c("25m_DAH_3Class.tif","25m_LandformClass_Default_Seive4.tif",
                          "25m_MRVBF_Classified_IS64Low6Up2.tif","dem.tif"), sep = "/")# ,"DEM_25m.tif"
layerNames <- c("DAH","LFC","MRVBF","DEM","cost") ##need to change this if you change the layers
ancDat <- raster::stack(covars)
proj4string(ancDat) <- "+init=epsg:3005"

# read in slope data
slope_raster <-  grep("^slope", list.files(covLoc))
slope <- raster(list.files(covLoc, full.name = TRUE)[slope_raster])
proj4string(slope) <- "+init=epsg:3005"

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
rSpd[,Values := (1/(speed/60))/40]
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
  if(x[1]|x[2] < 10){
    return(min(x[1],x[2]))
  }else{
    return(x[1]-x[2])
  }
}

rdIdx <- which(values(altAll) < 5)
slpIdx <- which(values(altAll) > 500)
adj <- adjacent(altAll, cells = slpIdx, pairs = T, directions = 8)
# adj <- adj[!adj[,1] %in% rdIdx ,]
# adj <- adj[!adj[,2] %in% rdIdx ,]

tr <- transition(altAll,trFn,directions = 8, symm = F)
tr1 <- geoCorrection(tr)
tr1[adj] <- (1.5*exp(-3.5*abs(tr1[adj] + 0.1)))/60 ##tobler's hiking function
tr1 <- geoCorrection(tr1)

start <- st_read("SmithersStart.gpkg") %>% as("Spatial")
acost <- accCost(tr1,start)
plot(acost)
writeRaster(acost, "TestAcost.tif", "GTiff")
