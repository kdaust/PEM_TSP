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

KL_stat_clhs <- function(sample){
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
  
  kl2 <- KL(tempMat) ## return KL divergence
  js2 <- JSD(tempMat) ## returns JS divergenxe
  temp <- t(tempMat)# %>% dplyr::select()
  ks2 <- ks.test(temp[,1], temp[,2], alternative = "g")
  ks2p <- ks2$p.value## returns KS statistic and p value
  return(kl2, ks2p, js2, js2, ks2$statistic)
}

KS_stat_clhs <- function(sample){
  sampDat <- as.data.table(sample)
  sampDat <- na.omit(sampDat)
  nsamp <- nrow(sampDat)
  
  ancBin <- ancVals[,lapply(.SD, function(x){cut(x,breaks = seq(min(x),max(x),length.out = 15),labels = F,include.lowest = T)})]
  ancBin <- copy(unite(ancBin, LHSVar, sep = "_"))
  ancSmall <- ancBin[,.(Num = .N), by = .(LHSVar)] ##count occurences
  ancSmall <- ancSmall[Num > 150,] ##remove uncommon
  setorder(ancSmall,-Num)
  ancSmall[,ID := seq_along(Num)]
  ancBin[ancSmall, ID := i.ID, on = "LHSVar"]
  ancBin <- na.omit(ancBin)
  
  fullBreaks <- ancVals[,lapply(.SD,function(x){seq(min(x),max(x),length.out = 15)})]
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
  #fullDat[,Diff := abs(Num - SampNum)]
  t1 <- cumsum(fullDat$SampNum)
  t2 <- cumsum(fullDat$Num)
  return(max(abs(t2-t1)))
  #return(mean(fullDat$Diff))
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