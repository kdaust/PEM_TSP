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

SLcov <- c("max_fp_l.tif","rid_level.tif","twi.tif","mrvbf.tif","convergence.tif","flow_accum_ft.tif",
           "aspect.tif","dem.tif")
covars <- paste(covLoc, SLcov, sep = "/")
ancDatSL <- raster::stack(covars)
proj4string(ancDatSL) <- "+init=epsg:3005"

### BGC1
bgc <- st_read(dsn = datLocGit, layer = "Boundary_BGC_dissolved")
bgcSub <- bgc[bgc$MAP_LABEL == "IDFdm1",]
rast <- raster(covars[1])
rast <- fasterize(bgcSub, rast)

ancDatSL <- mask(ancDatSL,rast)
fullSet <- sampleRegular(ancDatSL, size = 1000000,useGDAL = T)
X <- fullSet
X <- X[complete.cases(X),]

createLHS <- function(X,n){
  temp <- c_clhs(X,size = n, i_cost = NULL, iter = 5000)
  small <- as.matrix(temp$sampled_data)
  return(small)
}

collectStats <- function(X,sampleFun,nums,testFun1 = sumKSTest,testFun2 = fillTest, testFun3 = adTest){
  test_res <- foreach(n = nums, .combine = rbind, .noexport = c("CppLHS"), 
                      .packages = c("Rcpp","LaplacesDemon","foreach","goftest"),
                      .export = c("c_clhs")) %dopar% {
                        smallSet <- sampleFun(X,n)
                        res1 <- testFun1(X[,1:6],smallSet)
                        res2 <- testFun2(X[,1:6], smallSet)
                        #res3 <- testFun3(X[,1:6], smallSet)
                        data.frame(Num = n,KS = res1, fill = res2)
                      }
  return(test_res)
}

nums <- round(seq(10, 800, length.out = 20))
nums <- rep(nums, each = 5)