##PEM TSP functions
fillTest <- function(full, small){
  strata <- full[,lapply(.SD,function(x) {quantile(x, probs = seq(0, 1, 0.005), na.rm = TRUE)})]
  
  fillQual <- foreach(var = 1:ncol(full),.combine = cbind) %do% {
    .Call(graphics:::C_BinCount, small[[var]], strata[[var]], TRUE,TRUE)
  }
  
  return(length(fillQual[fillQual == 0])) 
}

fillTest_Transect <- function(full, small){
  strata <- lapply(
    full$rasterbands, 
    function(x) {
      quantile(x, probs = seq(0, 1, 0.01), na.rm = TRUE)
    }
  )
  fillQual <- foreach(var = 1:ncol(small),.combine = cbind) %do% {
    .Call(graphics:::C_BinCount, small[,var], strata[[var]], TRUE,TRUE)
  }
  
  return(length(fillQual[fillQual == 0])) 
}

Tri_build <- function(id, x, y){
  tris <- CreateRegularPolygon(3, c(x,y), 145) # number of sides, center pt and length of sides
  tris <- tris[c(1:3,1),]
  lines <- st_linestring(tris)
  
  g = st_sfc(lines)
  out <- st_sf(g,ID = id)
  #st_set_crs(newproj)
  return(out)
} 

KL_stat_clhs <- function(sample,ancSmall){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 9)})]
  ##split into same bins as full data set
  for(i in 1:ncol(sampDat)){
    sampDat[[i]] <- cut(sampDat[[i]],breaks = fullBreaks[[i]], labels = F, include.lowest = T)
  }
  sampBin <- copy(unite(sampDat,sampLHS,sep = "_"))
  sampBin[ancSmall, ID := i.ID, on = c(sampLHS = "LHSVar")]
  sampBin <- na.omit(sampBin)
  sampBin <- sampBin[,.(SampNum = .N), by = .(sampLHS)]##count
  #merge to compare
  fullDat <- copy(ancSmall)
  fullDat[sampBin, SampNum := i.SampNum, on = c(LHSVar = "sampLHS")]
  fullDat[is.na(SampNum), SampNum := 0]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  fullDat[,Num := Num/sum(Num)]
  tempMat <- rbind(fullDat$Num,fullDat$SampNum)
  
  kl2 <- suppressMessages(KL(tempMat)) ## return KL divergence
  return(kl2) #, ks2p, js2, js2, ks2$statistic
}

KL_stat_noextra <- function(sample,ancSmall){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 9)})]
  ##split into same bins as full data set
  for(i in 1:ncol(sampDat)){
    sampDat[[i]] <- cut(sampDat[[i]],breaks = fullBreaks[[i]], labels = F, include.lowest = T)
  }
  sampBin <- copy(unite(sampDat,sampLHS,sep = "_"))
  sampBin[ancSmall, ID := i.ID, on = c(sampLHS = "LHSVar")]
  sampBin <- na.omit(sampBin)
  sampBin <- sampBin[,.(SampNum = .N), by = .(sampLHS)]##count
  #merge to compare
  fullDat <- copy(ancSmall)
  fullDat[sampBin, SampNum := i.SampNum, on = c(LHSVar = "sampLHS")]
  fullDat[is.na(SampNum), SampNum := 0]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  fullDat[,Num := Num/sum(Num)]
  fullDat[Num < SampNum, SampNum := Num]
  fullDat[Num > SampNum, SampNum := SampNum/sum(SampNum)]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  tempMat <- rbind(fullDat$Num,fullDat$SampNum)
  
  kl2 <- suppressMessages(KL(tempMat)) ## return KL divergence
  return(kl2) #, ks2p, js2, js2, ks2$statistic
}

binfill_clhs <- function(sample,ancSmall){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 9)})]
  ##split into same bins as full data set
  for(i in 1:ncol(sampDat)){
    sampDat[[i]] <- cut(sampDat[[i]],breaks = fullBreaks[[i]], labels = F, include.lowest = T)
  }
  sampBin <- copy(unite(sampDat,sampLHS,sep = "_"))
  sampBin[ancSmall, ID := i.ID, on = c(sampLHS = "LHSVar")]
  sampBin <- na.omit(sampBin)
  sampBin <- sampBin[,.(SampNum = .N), by = .(sampLHS)]##count
  #merge to compare
  fullDat <- copy(ancSmall)
  fullDat[sampBin, SampNum := i.SampNum, on = c(LHSVar = "sampLHS")]
  fullDat[is.na(SampNum), SampNum := 0]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  fullDat[,Num := Num/sum(Num)]
  fullDat[,CumAmount := cumsum(Num)]
  top80 <- fullDat[CumAmount < 0.5,]
  numFull <- nrow(top80[SampNum > 0,])
  return(numFull/nrow(top80))
}

plotUniDists <- function(sample){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 9)})]
  ##split into same bins as full data set
  for(i in 1:ncol(sampDat)){
    sampDat[[i]] <- cut(sampDat[[i]],breaks = fullBreaks[[i]], labels = F, include.lowest = T)
  }
  sampBin <- copy(unite(sampDat,sampLHS,sep = "_"))
  sampBin[ancSmall, ID := i.ID, on = c(sampLHS = "LHSVar")]
  sampBin <- na.omit(sampBin)
  sampBin <- sampBin[,.(SampNum = .N), by = .(sampLHS)]##count
  #merge to compare
  fullDat <- copy(ancSmall)
  fullDat[sampBin, SampNum := i.SampNum, on = c(LHSVar = "sampLHS")]
  fullDat[is.na(SampNum), SampNum := 0]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  fullDat[,Num := Num/sum(Num)]
  tempMat <- rbind(fullDat$Num,fullDat$SampNum)
  kl2 <- KL(tempMat)
  return(list(data = fullDat, KL = kl2))
  
}

Tri_build <- function(id, x, y){
  tris <- CreateRegularPolygon(3, c(x,y), 145) # number of sides, center pt and length of sides
  tris <- tris[c(1:3,1),]
  lines <- st_linestring(tris)
  
  g = st_sfc(lines)
  out <- st_sf(g,ID = id)
  #st_set_crs(newproj)
  return(out)
} 

rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

## Feature rotation
rotFeature <- function(Feature, PivotPt, Bearing) {
  # where Feature is the Shape to be rotated, eg:  #Feature <- tri
  # Bearing is the compass bearing to rotate to    #PivotPt <- pt.sf
  # PivotPt is the point to rotate around          #Bearing <- 15 
  
  ## extract the geometry
  Feature_geo <- st_geometry(Feature)
  PivotPoint  <- st_geometry(PivotPt)
  
  ## Convert bearing from degrees to radians
  d <- ifelse(Bearing > 180, pi * ((Bearing-360)/ 180) ,  pi * (Bearing / 180))
  
  rFeature <- (Feature_geo - PivotPoint) * rot(d)   + PivotPoint
  rFeature <- st_set_crs(rFeature, st_crs(Feature))
  
  Feature$geometry <- st_geometry(rFeature) ## replace the original geometry
  return(Feature)
}



#sample = s1
plotUniDists <- function(sample){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  ancBin <- ancVals[,lapply(.SD, function(x){cut(x,breaks = seq(min(x),max(x),length.out = 9),labels = F,include.lowest = T)})]
  ancBin <- copy(unite(ancBin, LHSVar, sep = "_"))
  ancSmall <- ancBin[,.(Num = .N), by = .(LHSVar)] ##count occurences
  ancSmall <- ancSmall[Num > 150,] ##remove uncommon
  setorder(ancSmall,-Num)
  ancSmall[,ID := seq_along(Num)]
  ancBin[ancSmall, ID := i.ID, on = "LHSVar"]
  ancBin <- na.omit(ancBin)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 9)})]
  ##split into same bins as full data set
  for(i in 1:ncol(sampDat)){
    sampDat[[i]] <- cut(sampDat[[i]],breaks = fullBreaks[[i]], labels = F, include.lowest = T)
  }
  sampBin <- copy(unite(sampDat,sampLHS,sep = "_"))
  sampBin[ancSmall, ID := i.ID, on = c(sampLHS = "LHSVar")]
  sampBin <- na.omit(sampBin)
  sampBin <- sampBin[,.(SampNum = .N), by = .(sampLHS)]##count
  #merge to compare
  fullDat <- copy(ancSmall)
  fullDat[sampBin, SampNum := i.SampNum, on = c(LHSVar = "sampLHS")]
  fullDat[is.na(SampNum), SampNum := 0]
  fullDat[,SampNum := SampNum/sum(SampNum)]
  fullDat[,Num := Num/sum(Num)]
  tempMat <- rbind(fullDat$Num,fullDat$SampNum)
  kl2 <- KL(tempMat)
  return(list(data = fullDat, KL = kl2))
  
}

##create TSP to calculate time
calcCost_clhs <- function(pnts,objVals,plotTime = 50L, minPerDay = 3L){
  n = nrow(pnts)
  p2 <- st_as_sf(pnts)
  pnts <- pnts[,1]
  colnames(pnts) <- c("name","geometry")
  st_geometry(pnts) <- "geometry"
  startPnts <- st_as_sf(data.frame(name = "Start",geometry = start_sf))
  pnts <- rbind(pnts, startPnts)
  pnts2 <- as(pnts, "Spatial")
  
  ## create distance matrix between sample points
  test <- costDistance(trSmall,pnts2,pnts2)
  dMat2 <- as.matrix(test)
  dMat2 <- dMat2*60
  dMat2[is.infinite(dMat2)] <- 1000
  
  ##penalty based on quality of points
  objVals <- max(objVals) - objVals
  
  maxTime <- 8L ##hours
  ## time per transect
  temp <- dMat2[1:n,1:n]
  maxDist <- sum(temp[upper.tri(temp)])
  minPen <- maxDist
  maxPen <- maxDist * 4
  objVals <- scales::rescale(objVals, to = c(minPen,maxPen))
  objVals <- as.integer(objVals)
  
  ndays <- as.integer(ceiling(n/minPerDay)+1)
  
  pen = objVals
  
  indStart <- as.integer(rep(n,ndays))
  ##run vehicle routing problem from python script
  ## GCS is global span cost coefficient
  vrp <- py_mTSP(dat = dMat2,num_days = ndays, start = indStart, end = indStart, 
                 max_cost = maxTime*60L, plot_time = plotTime, penalty =  pen, arbDepot = F, GSC = 5L)
  time <- vrp[[2]]
  totTime <- sum(as.numeric(unlist(time)))
  return(totTime)
}