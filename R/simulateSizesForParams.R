#'Get params based on simData
#'
#'@param d a model run dataFrame
#'@return a list
#'@export
#'
updateSimData <- function(simData){
  simData$isYOY <- c(rep(0,simData$nOccYOY),rep(1,simData$nOccAdult))
  simData$seasons <- c(3,4,1,2,3,4,1,2,3,4) # will automate depending on nOcc_
  simData$itersForSim <- sort(sample(1:simData$nSamplesAvail, simData$nReps))
  simData$nOcc <- simData$nOccYOY + simData$nOccAdult

  simData$lengthMean <- dGIn[[simData$species]][[1]]$constants$lengthMean
  simData$lengthSD <- dGIn[[simData$species]][[1]]$constants$lengthSD

  simData$sampInt <- c(70,122,70,101,70,122,70,101,70,122)
  # d2 %>%
  # group_by(season) %>%
  # summarise(sampInt = mean(sampleInterval))

  return(simData)
}

#'Get params based on simData
#'
#'@param d a model run dataFrame
#'@return a list
#'@export
#'
getNBetas <- function(){

  nBetas <- dim(mcmc[[simData$species]]$sims.list$grBeta)[2]
  nBetasBNT <- dim(mcmc[[simData$species]]$sims.list$grBetaBNT)[2]
  nBetasATS <- dim(mcmc[[simData$species]]$sims.list$grBetaATS)[2]
  nSigmaBetas <- dim(mcmc[[simData$species]]$sims.list$sigmaBeta)[2]
#sizeMat <- intMat <- betaMat <- isYOYMat <- seasonMat <- matrix(NA, nrow = simData$numInd, ncol = simData$nOccYOY + simData$nOccAdult)
  return(list(
    nBetas=nBetas,nBetasBNT=nBetasBNT,nBetasATS=nBetasATS,nSigmaBetas=nSigmaBetas
  ))
}


#'Create grid of input variables
#'
#'@param d a model run dataFrame
#'@return a list of matrices
#'@export
#'
createMatrices <- function(){

  sizeMat <- intMat <- sigmaIntMat <- isYOYMat <-  array(NA, c(simData$nReps, simData$numInd, simData$nOcc))
  seasonMat <- sampIntMat <- indREMat <-  grMat <- array(NA, c(simData$nReps, simData$numInd, simData$nOcc))
  betaMat <- sigmaBetaMat <- array(NA, c(simData$nReps, nB$nBetas, simData$numInd, simData$nOcc))
  betaMatBNT <- array(NA, c(simData$nReps, nB$nBetas, simData$numInd, simData$nOcc))
  betaMatATS <- array(NA, c(simData$nReps, nB$nBetas, simData$numInd, simData$nOcc))

  indRE <- rnorm(simData$numInd * simData$nReps,0,mcmc[[simData$species]]$mean$sigmaIndRE)

  # isYOY
  isYOYMat[,,1:simData$nOccYOY] <- 1
  isYOYMat[,,(simData$nOccYOY + 1):(simData$nOcc)] <- 2

  #season and sampleInterval
  for(i in 1:simData$nReps){
    for(r in 1:simData$numInd){
      seasonMat[i,r,] <- simData$seasons
      sampIntMat[i,r,] <- simData$sampInt
    }
  }

  #size
  sizeMat[,,1] <- rnorm(simData$numInd * simData$nReps, simData$meanInit, simData$sdInit)

  #indRE, fills by i, then rows
  indREMat[] <- indRE

  # ints and betas
  for(i in 1:simData$nReps){
    for(rows in 1:simData$numInd){
      for(cols in 1:simData$nOcc){

        spp <- simData$species
        iter <- simData$itersForSim[i]
        y <- isYOYMat[i,rows,cols]
        s <- seasonMat[i,rows,cols]
        r <- simData$riverIndex

        intMat[i,rows,cols] <- mcmc[[spp]]$sims.list$grInt[iter,y,s,r]
        sigmaIntMat[i,rows,cols] <- mcmc[[spp]]$sims.list$grInt[iter,y,s,r]

        for(b in 1:nB$nBetas) {    betaMat[i,b,rows,cols] <-    mcmc[[spp]]$sims.list$grBeta[iter,b,y,s,r] }
        for(b in 1:nB$nBetasBNT) { betaMatBNT[i,b,rows,cols] <- mcmc[[spp]]$sims.list$grBetaBNT[iter,b,y,s,r] }
        for(b in 1:nB$nBetasATS) { betaMatATS[i,b,rows,cols] <- mcmc[[spp]]$sims.list$grBetaATS[iter,b,y,s,r] }
        for(b in 1:nB$nSigmaBetas) { sigmaBetaMat[i,b,rows,cols] <-    mcmc[[spp]]$sims.list$sigmaBeta[iter,b,y,s,r] }

      }
    }
  }

  return(list(isYOYMat=isYOYMat,seasonMat=seasonMat,sampIntMat=sampIntMat,sizeMat=sizeMat,grmat=grMat,indREMat=indREMat,
              intMat=intMat,sigmaIntMat=sigmaIntMat,betaMat=betaMat,
              betaMatBNT=betaMatBNT,betaMatATS=betaMatATS,sigmaBetaMat=sigmaBetaMat))
}

#'Create grid of input variables
#'
#'@param d a model run dataFrame
#'@return a data frame
#'@export

getIndVars <- function(){

  indVars <- expand.grid( temp = simData$indVarSeq,
                          flow = simData$indVarSeq,
                          bkt =  simData$indVarSeq,
                          bnt =  simData$indVarSeq,
                          ats =  simData$indVarSeq
  )
  indVars$rowNum <- 1:nrow(indVars)
  return(indVars)
}



#'Generate predictions of body sizes
#'
#'@param d a model run dataFrame
#'@return a matrix
#'@export
#'
projectSizes <- function(m,simData,testRow){

  d <- indVars %>% filter(rowNum == testRow)

  for(i in 1:simData$nReps){
    print(paste0(i," out of ",simData$nReps,"  ",speciesForProj))
    print(d)
    for(rows in 1:simData$numInd){
      for(cols in 1:(simData$nOcc-1)){

        gr <- m$intMat[i,rows,cols] +
          m$betaMat[i,1,rows,cols] * d$temp +
          m$betaMat[i,2,rows,cols] * d$flow +
          m$betaMat[i,3,rows,cols] * d$temp * d$flow +
          m$betaMat[i,4,rows,cols] * d$temp^2 +
          m$betaMat[i,5,rows,cols] * d$flow^2 +
          m$betaMat[i,6,rows,cols] * d$bkt +
          m$betaMat[i,7,rows,cols] * ((m$sizeMat[i,rows,cols] - simData$lengthMean) / simData$lengthSD) +
          m$betaMatBNT[i,1,rows,cols] * d$bnt +
          m$betaMatBNT[i,2,rows,cols] * d$bnt * d$bkt +
          m$betaMatATS[i,1,rows,cols] * d$ats +
          m$indREMat[i,rows,cols]
        #temp + flow + bkt + bnt + ats + temp * flow + temp^2 + flow^2 + bnt*bkt
        sigma <- exp(
          m$sigmaIntMat[i,rows,cols] +
            m$sigmaBetaMat[i,1,rows,cols] * d$temp +
            m$sigmaBetaMat[i,2,rows,cols] * d$flow +
            m$sigmaBetaMat[i,3,rows,cols] * d$temp * d$flow
        )
        # print(c(sigma,exp(sigma)))
        m$sizeMat[i,rows,cols + 1] <-  rnorm(1, m$sizeMat[i,rows,cols] + gr * m$sampIntMat[i,rows,cols], sigma)
        m$grMat[i,rows,cols] <- gr
      }
    }
  }
  return( list(sizeMat = m$sizeMat, grMat = m$grMat) )

}


#'Summarize predictions of body sizes
#'
#'@param d a model run dataFrame
#'@return a matrix
#'@export
#'
summarizeLens <- function(sizeMat,testRow){
  predLen <- array2df(sizeMat, levels = list(iter=1:simData$nReps,ind=1:simData$numInd,occ=1:simData$nOcc), label.x="len")

  p <- predLen %>%
    group_by(iter,occ) %>%
    summarize( lenMean = mean(len),
               lenSD = sd(len),
               n=n()) %>%
    ungroup()

  d <- indVars %>% filter(rowNum == testRow)
  p <- cbind(p,d)

  p2 <- p %>%
    group_by(occ) %>%
    summarize( lenMeanMean = mean(lenMean),
               lenMeanSD = sd(lenMean),
               lenSDMean = mean(lenSD),
               lenSDSD = sd(lenSD),
               n=n()) %>%
    ungroup()

  p2 <- cbind(p2,d)

  return(list(p,p2))
}

