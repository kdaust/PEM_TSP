### check variable coverage of sample points
library(sf)
library(raster)
library(sp)
library(gdistance)
library(foreach)
library(data.table)
library(fasterize)
library(reticulate)
library(here)
library(LearnGeom)
library(tidyverse)
library(goftest)
library(stars)
source("FastCLHS_R.R")


##Kolmogorov-Smirnov test
sumKSTest <- function(full,small){
  out <- 0
  for(i in 1:ncol(full)){
    ks <- suppressWarnings(ks.test(full[,i],small[,i]))
    out <- out+ks$statistic
  } 
  return(out)
}

adTest <- function(full,small){
  out <- 0
  for(i in 1:ncol(full)){
    fn <- ecdf(full[,i])
    temp <- ad.test(small[,i], null = fn)
    out <- out + temp$statistic
  }
  return(out)
}

fillTest <- function(full, small){
  strata <- apply(
    full, 
    2, 
    function(x) {
      quantile(x, probs = seq(0, 1, 0.01), na.rm = TRUE)
    }
  )
  fillQual <- foreach(var = 1:ncol(full),.combine = cbind) %do% {
    .Call(graphics:::C_BinCount, small[,var], strata[,var], TRUE,TRUE)
  }
  
  return(length(fillQual[fillQual == 0])) 
}

datLocGit <- here("InputData") ## Data
covLoc <- "/media/kiridaust/MrBig/boundary/BoundaryTSA_AOI/1_map_inputs/covariates/5m/"
covLoc <- "E:/boundary/BoundaryTSA_AOI/1_map_inputs/covariates/5m/"

SLcov <- c("rid_level.tif","twi.tif","mrvbf.tif","convergence.tif","flow_accum_ft.tif",
           "aspect.tif","dem.tif")
covars <- paste(covLoc, SLcov, sep = "/")
ancDatSL <- raster::stack(covars)
proj4string(ancDatSL) <- "+init=epsg:3005"

ancDat <- read_stars(covars, proxy = T)

### BGC1
bgc <- st_read(dsn = datLocGit, layer = "Boundary_BGC_dissolved")
bgcSub <- bgc[bgc$MAP_LABEL == "IDFdm1",]
ancDat <- ancDat[bgcSub]
test <- st_as_stars(ancDat)
rast <- raster(covars[1])
rast <- fasterize(bgcSub, rast)

ancDatSL <- mask(ancDatSL,rast)

fullSet <- sampleRegular(ancDatSL, size = 1000000,useGDAL = T)
X <- fullSet
X <- X[complete.cases(X),]
X <- X[,-8]

createLHS <- function(X,n){
  temp <- c_clhs(X,size = n, i_cost = NULL, iter = 5000)
  small <- as.matrix(temp$sampled_data)
  return(small)
}

createLHS_30m <- function(X,n,ancDat){
  xMat <- as.matrix(st_drop_geometry(X))
  temp <- c_clhs(xMat,size = n,i_cost = ncol(xMat),iter = 5000)
  dat <- X[temp$indeces,]
  small <- foreach(j = 1:nrow(dat),.combine = rbind) %do% {
    tri <- Tri_build(id = j, x = st_coordinates(dat[j,])[1],y = st_coordinates(dat[j,])[2])
    tPs <- st_sample(tri, size = 25, type = "regular")
    temp <- suppressWarnings(st_sf(data.frame(id = j, geometry = tPs)) %>%
                               st_cast("POINT"))
    temp
    
  }
  temp <- suppressWarnings(raster::extract(ancDat, small))
  temp <- temp[!is.na(temp[1,]),]
  return(temp)
}

worker.init <- function(){
  Rcpp::sourceCpp("CppLHS.cpp")
}

require(doParallel)
cl <- makePSOCKcluster(detectCores()-2)
clusterCall(cl, worker.init)
registerDoParallel(cl)

collectStats <- function(X,sampleFun,nums,testFun1 = sumKSTest,testFun2 = fillTest, testFun3 = adTest){
  test_res <- foreach(n = nums, .combine = rbind, .noexport = c("CppLHS"), 
                      .packages = c("Rcpp","LaplacesDemon","foreach","goftest"),
                      .export = c("c_clhs")) %dopar% {
                        smallSet <- sampleFun(X,n)
                        res1 <- testFun1(X,smallSet)
                        res2 <- testFun2(X, smallSet)
                        #res3 <- testFun3(X[,1:6], smallSet)
                        data.frame(Num = n,KS = res1, fill = res2)
                      }
  return(test_res)
}

nums <- round(seq(10, 1500, length.out = 20))
nums <- rep(nums, each = 5)
statsLHS <- collectStats(X,createLHS, nums)
boxplot(fill ~ Num, data = statsLHS)
boxplot(KS ~ Num, data = statsLHS)

sampleDat <- st_read(dsn = paste0(datLocGit,"/transect1_30m_pts_att.gpkg"))
sampleDat <- sampleDat["id"]
sampleDat <- sampleDat[grep("MSdm1",sampleDat$id),]
temp <- raster::extract(ancDatSL,sampleDat)

orig <- as.data.table(X)
orig <- orig[complete.cases(orig),]
smallSet <- as.data.table(temp)
smallSet <- smallSet[complete.cases(smallSet),]
out <- foreach(pos = 1:ncol(orig), .combine = rbind) %do% {
  densOrig <- density(orig[[pos]])
  densSample <- density(smallSet[[pos]])
  densDat <- data.table(x = densOrig$x, y = densOrig$y, Type = "Orig")
  tempDat <- data.table(x = densSample$x, y = densSample$y, Type = "Sample")
  densDat <- rbind(densDat, tempDat)
  densDat[,Var := names(smallSet)[pos]]
  densDat
}

ggplot(out, aes(x = x, y = y, group = Type, colour = Type))+
  geom_line()+
  facet_wrap(.~Var, scales = "free")

temp <- temp[,-8]
sumKSTest(X,temp)
fillTest(X,temp)
