
#'Generate predictions from a set of iterations of a jags model run
#'
#'@param d a model run dataFrame
#'@param limits symetrical lower and upper limit for the predictions
#'@param nPoints the number of points between the limits. Use an odd # to ensure a value of 0
#'@param itersForPred, the iterations from the model for prediction
#'@param varsToEstimate which variables to use for estimates (len,count,phi,temp,flow)
#'@return a data frame
#'@export

getPrediction <- function(d, limits = 2, nPoints = 5, itersForPred, constants, sampleIntervalMeanIn, varsToEstimate, ii=1){

  # get grInt in df format
  #grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
  grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="int")

  # get grBeta in df format and merge in grInt
  grBeta3 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="est")

  grBeta2 <- spread( grBeta3, key = beta, value = est, sep = "" ) %>%
    left_join( .,grInt )

  # grBetaBNT, 2 rivers only
  grBetaBNT <- array2df(d$sims.list$grBetaBNT, levels = list(iter=NA,betaBNT=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:(4-(4-constants$nRivers))]), label.x="est")

  grBeta1 <- spread( grBetaBNT, key = betaBNT, value = est, sep = "" ) %>%
    left_join( .,grBeta2 )

  # grBetaATS, no river
  grBetaATS <- array2df(d$sims.list$grBetaATS, levels = list(iter=NA,betaATS=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:(4-(4-constants$nRivers))]), label.x="est")

  grBeta <- spread( grBetaATS, key = betaATS, value = est, sep = "" ) %>%
    left_join( .,grBeta1 )


  sampleIntervalMean <- array2df(sampleIntervalMeanIn, levels = list( river=riverOrderedIn[1:constants$nRivers],season=1:constants$nSeasons,species=1:constants$nSpecies), label.x = "interval")
  grBeta <- left_join(grBeta,sampleIntervalMean)

  # prediction template
  x <- seq( -limits,limits,length.out = nPoints )

  # create inital template of standardized values to predict over
  expIndex <- 0

  if ("temp" %in% varsToEstimate) {
    tempData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    tempData <- 0

  if ("flow" %in% varsToEstimate) {
    flowData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    flowData <- 0

  if ("phi" %in% varsToEstimate) {
    phiData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
  phiData <- 0

  if ("count" %in% varsToEstimate) {
    countData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
  countData <- 0

  if ("len" %in% varsToEstimate) {
    lenData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
   }  else
   lenData <- 0

  if ("cBKT" %in% varsToEstimate) {
    cBKTData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cBKTData <- 0

  if ("cBNT" %in% varsToEstimate) {
    cBNTData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cBNTData <- 0

  if ("cATS" %in% varsToEstimate) {
    cATSData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cATSData <- 0

  predTemplate <- data.frame(
                            #  len = lenData,
                             # count = countData,
                              flow =  flowData,
                              temp =  tempData,
                              cBKT =  cBKTData,
                              cBNT =  cBNTData,
                              cATS =  cATSData
                            )


    # expand predTemplate across grBeta rows for a given set of iterations
    grBetaIter <- grBeta %>% filter( iter %in% itersForPred )

    # repeat predTemplate nrow( grBetaIter ) times
    predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
    # repeat each row of grBetaIter nrow(predTemplate) times
    grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
    # Put them together
    preds <- cbind( predTemplateLong,grBetaLong )

  # This model structure needs to match that in grModelx.jags
  if(speciesGr == 'ats'){
    preds <- preds %>%
    mutate( predGr =
            ( int +

              beta1 * temp +
              beta2 * flow +
              beta3 * temp * flow +

              beta4 * temp^2 +
              beta5 * flow^2 +

              beta6 * cBKT +

              betaBNT1 * cBNT +
           #   betaBNT2 * cBNT^2 +
          #    betaBNT3 * cBNT * cBKT +

              betaATS1 * cATS
           #   betaATS2 * cATS^2

            ) #* interval
    ) %>%
    dplyr::select( flow,temp,cBKT,cBNT,cATS,iter,isYOY,season,river,predGr )

  } else {
  preds <- preds %>%
    mutate( predGr =
              ( int +

                  beta1 * temp +
                  beta2 * flow +
                  beta3 * temp * flow +

                  beta4 * temp^2 +
                  beta5 * flow^2 +

                  beta6 * cBKT +

                  betaBNT1 * cBNT +
                  betaBNT2 * cBNT^2 +
                  betaBNT3 * cBNT * cBKT +

                  betaATS1 * cATS

              ) #* interval
    ) %>%
    dplyr::select( flow,temp,cBKT,cBNT,cATS,iter,isYOY,season,river,predGr )
  }
  return(preds)
}



#'Generate predictions of sigma the sd in growth from a set of iterations of a jags model run
#'
#'@param d a model run dataFrame
#'@param limits symetrical lower and upper limit for the predictions
#'@param nPoints the number of points between the limits. Use an odd # to ensure a value of 0
#'@param itersForPred, the iterations from the model for prediction
#'@return a data frame
#'@export

getPredictionSigma <- function(d, limits = 2, nPoints = 5, itersForPred, constants, sampleIntervalMeanIn, varsToEstimate, ii=1){

  # get grInt in df format
  #grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
  grInt <- array2df(d$sims.list$sigmaInt, levels = list(iter=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="int")

  # get grBeta in df format and merge in grInt
  grBeta3 <- array2df(d$sims.list$sigmaBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="est")

  grBeta2 <- spread( grBeta3, key = beta, value = est, sep = "" ) %>%
    left_join( .,grInt )

  # grBetaBNT, 2 rivers only
  grBetaBNT <- array2df(d$sims.list$grBetaBNT, levels = list(iter=NA,betaBNT=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:(4-(4-constants$nRivers))]), label.x="est")

  grBeta1 <- spread( grBetaBNT, key = betaBNT, value = est, sep = "" ) %>%
    left_join( .,grBeta2 )

  # grBetaATS, no river
  grBetaATS <- array2df(d$sims.list$grBetaATS, levels = list(iter=NA,betaATS=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:(4-(4-constants$nRivers))]), label.x="est")

  grBeta <- spread( grBetaATS, key = betaATS, value = est, sep = "" ) %>%
    left_join( .,grBeta1 )


  sampleIntervalMean <- array2df(sampleIntervalMeanIn, levels = list( river=riverOrderedIn[1:constants$nRivers],season=1:constants$nSeasons,species=1:constants$nSpecies), label.x = "interval")
  grBeta <- left_join(grBeta,sampleIntervalMean)

  # prediction template
  x <- seq( -limits,limits,length.out = nPoints )

  # create inital template of standardized values to predict over
  expIndex <- 0

  if ("temp" %in% varsToEstimate) {
    tempData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    tempData <- 0

  if ("flow" %in% varsToEstimate) {
    flowData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    flowData <- 0

  if ("phi" %in% varsToEstimate) {
    phiData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    phiData <- 0

  if ("count" %in% varsToEstimate) {
    countData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    countData <- 0

  if ("len" %in% varsToEstimate) {
    lenData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    lenData <- 0

  if ("cBKT" %in% varsToEstimate) {
    cBKTData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cBKTData <- 0

  if ("cBNT" %in% varsToEstimate) {
    cBNTData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cBNTData <- 0

  if ("cATS" %in% varsToEstimate) {
    cATSData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    cATSData <- 0

  predTemplate <- data.frame(
  #  len = lenData,
    # count = countData,
    flow =  flowData,
    temp =  tempData
  #  cBKT =  cBKTData,
  #  cBNT =  cBNTData,
  #  cATS =  cATSData
  )


  # expand predTemplate across grBeta rows for a given set of iterations
  grBetaIter <- grBeta %>% filter( iter %in% itersForPred )

  # repeat predTemplate nrow( grBetaIter ) times
  predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
  # repeat each row of grBetaIter nrow(predTemplate) times
  grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
  # Put them together
  preds <- cbind( predTemplateLong,grBetaLong )

  # This model structure needs to match that in grModelx.jags
  if(speciesGr == 'ats'){
  preds <- preds %>%
    mutate( predGrSigma =
              exp( int

                #  beta1 * temp
                 ## beta2 * flow +
                #  beta3 * temp * flow

              ) #* interval
    ) %>%
    dplyr::select( flow,temp,iter,isYOY,season,river,predGrSigma )

  } else {
  preds <- preds %>%
    mutate( predGrSigma =
              exp( int +

                     beta1 * temp +
                     beta2 * flow +
                     beta3 * temp * flow

              ) #* interval
    ) %>%
    dplyr::select( flow,temp,iter,isYOY,season,river,predGrSigma )
  }

  return(preds)

}




#'Generate predictions of sigma the sd in growth from a set of iterations of a jags model run
#'
#'@param d a model run dataFrame
#'@param limits symetrical lower and upper limit for the predictions
#'@param nPoints the number of points between the limits. Use an odd # to ensure a value of 0
#'@param itersForPred, the iterations from the model for prediction
#'@return a data frame
#'@export

getPredictionSigma_Counts <- function(d, limits = 2, nPoints = 5, itersForPred, constants, sampleIntervalMeanIn, varsToEstimate, ii=1){

  # get grInt in df format
  #grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
  grInt <- array2df(d$sims.list$sigmaInt, levels = list(iter=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="int")

  # get grBeta in df format and merge in grInt
  # grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta1 <- array2df(d$sims.list$sigmaBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),season=1:constants$nSeasons,river=riverOrderedIn[1:constants$nRivers]), label.x="est")

  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" ) %>%
    left_join( .,grInt )

  sampleIntervalMean <- array2df(sampleIntervalMeanIn, levels = list( river=riverOrderedIn[1:constants$nRivers],season=1:constants$nSeasons,species=1:constants$nSpecies), label.x = "interval")
  grBeta <- left_join(grBeta,sampleIntervalMean)
  # prediction template
  x <- seq( -limits,limits,length.out = nPoints )

  # create inital template of standardized values to predict over
  expIndex <- 0

  if ("temp" %in% varsToEstimate) {
    tempData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    tempData <- 0

  if ("flow" %in% varsToEstimate) {
    flowData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    flowData <- 0

  if ("phi" %in% varsToEstimate) {
    phiData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    phiData <- 0

  if ("count" %in% varsToEstimate) {
    countData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    countData <- 0

  if ("len" %in% varsToEstimate) {
    lenData <- rep(x, each = nPoints ^ expIndex)
    expIndex <- expIndex + 1
  }  else
    lenData <- 0

  predTemplate <- data.frame(
    len = lenData,
    count = countData,
    flow =  flowData,
    temp =  tempData
  )


  # expand predTemplate across grBeta rows for a given set of iterations
  grBetaIter <- grBeta %>% filter( iter %in% itersForPred )

  # repeat predTemplate nrow( grBetaIter ) times
  predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
  # repeat each row of grbetaIter nrow(predTemplate) times
  grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
  # Put them together
  preds <- cbind( predTemplateLong,grBetaLong )

  # This model structure needs to match that in grModelx.jags
  preds <- preds %>%
    mutate( predGrSigma =
              exp(
                int +
                beta1 * (count) +
                beta2 * (temp) +
                beta3 * (flow) +
                beta4 * (count) * (temp) +
                beta5 * (temp) * (flow) +
                beta6 * (count) * (flow)
              )
    )

  return(preds)
}




# getPredictionSigma <- function(d, limits = 2, nPoints = 5, itersForPred){
#
#   # get grInt in df format
#   grInt <- array2df(d$sims.list$sigmaInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
#
#   # get grBeta in df format and merge in grInt
#   grBeta1 <- array2df(d$sims.list$sigmaBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
#   grBeta <- spread( grBeta1, key = beta, value = est, sep = "" )  %>%
#     left_join( .,grInt )
#
#   # prediction template
#   x <- seq( -limits,limits,length.out = nPoints )
#
#   # create inital template of standardized values to predict over
#   predTemplate <- data.frame( #len =   rep(x, each = nPoints ^ 3),
#                               count = rep(x, each = nPoints ^ 2),
#                               flow =  rep(x, each = nPoints ^ 1),
#                               temp =      x
#   )
#
#
#   # expand predTemplate across grBeta rows for a given set of iterations
#   grBetaIter <- grBeta %>% filter( iter %in% itersForPred )
#
#   # repeat predTemplate nrow( grBetaIter ) times
#   predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
#   # repeat each row of grbetaIter nrow(predTemplate) times
#   grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
#   # Put them together
#   preds <- cbind( predTemplateLong,grBetaLong )
#
#   # This model structure needs to match that in grModel.jags
#   preds <- preds %>%
#     mutate( predGr = int +
#
#               beta1 * count +
#
#               beta2 * count^2 +
#
#               beta3 * temp +
#               beta4 * flow +
#               beta5 * temp^2 +
#               beta6 * flow^2 +
#               beta7 * temp * flow
#     )
#
#   return(preds)
# }

#'Plot predictions from a set of iterations of a jags model run
#'
#'@param p a model run set of predictions from getPredictions()
#'@param varsToPlot variables for ploting. Eith singles or pairs of (flow,temp,count)
#'@param isYOYGG, isYOY, 0 or 1
#'@param speciesGG, bkt, bnt or ats
#'@return a ggplot object
#'@export
#'

plotPred <- function(p, depVar, varsToPlot, isYOYGG, speciesGG) {
  all = c('temp','flow','cBKT','cBNT','cATS')

  #print(c(depVar,(as.name(depVar))))

  if(length(varsToPlot) == 1) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]
    notPlot[3] <- all[!(all %in% varsToPlot)][3]
    notPlot[4] <- all[!(all %in% varsToPlot)][4]
 #   notPlot[5] <- all[!(all %in% varsToPlot)][5]

 #   pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0 ) %>%
#                 distinct(eval(as.name(varsToPlot[1])), iter, isYOY, river, species, season, .keep_all = TRUE)
    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0,
                                          eval(as.name(notPlot[4])) == 0
                                        #, eval(as.name(notPlot[5])) == 0
                        ) %>%
      distinct(eval(as.name(varsToPlot[1])), iter, isYOY, river, season, .keep_all = TRUE)

    ggOut <- ggplot( pGG, aes( eval(as.name(varsToPlot[1])), eval(as.name(depVar)), group = iter,color=(iter) ) ) +
      geom_line( alpha=0.25 ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG,', species = ',speciesGG)) +
      geom_hline(yintercept = 0) +
    #  ylim(-0.25,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      facet_grid(river~season)
  }

  if(length(varsToPlot) == 2) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]
    notPlot[3] <- all[!(all %in% varsToPlot)][3]
  #  notPlot[4] <- all[!(all %in% varsToPlot)][4]

#    pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0 ) %>%
#      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), iter, isYOY, river, species, season, .keep_all = TRUE)

    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0,
                                          eval(as.name(notPlot[3])) == 0
                        #, eval(as.name(notPlot[4])) == 0
                        ) %>%
      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), iter, isYOY, river, season, .keep_all = TRUE)

    pGG$iterGroup <- paste0(pGG$iter,pGG[[varsToPlot[2]]])

    ggOut <- ggplot(pGG, aes(eval(as.name(varsToPlot[1])), eval(as.name(depVar)), group = iterGroup, name = varsToPlot[2])) +
      geom_line(aes( color = (eval(as.name(varsToPlot[2]))) ), alpha = 0.25 ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG,', species = ',speciesGG)) +
      geom_hline(yintercept = 0) +
      ylim(-0.33,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      #    scale_fill_continuous(name=varsToPlot[2]) +
      facet_grid(river~season)
  }
  return(ggOut)
}


#'Plot predictions from means of a jags model run
#'
#'@param p a model run set of predictions from getPredictions()
#'@param varsToPlot variables for ploting. Eith singles or pairs of (flow,temp,count)
#'@param isYOYGG, isYOY, 0 or 1
#'@return a ggplot object
#'@export
#'
plotPredMeans <- function(p, depVar, varsToPlot, isYOYGG) {
  all = c('temp','flow','cBKT','cBNT','cATS')

  #print(c(depVar,(as.name(depVar))))

  if(length(varsToPlot) == 1) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]
    notPlot[3] <- all[!(all %in% varsToPlot)][3]
    notPlot[4] <- all[!(all %in% varsToPlot)][4]


    #   pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0 ) %>%
    #                 distinct(eval(as.name(varsToPlot[1])), iter, isYOY, river, species, season, .keep_all = TRUE)
    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0,
                        eval(as.name(notPlot[4])) == 0) %>%
      distinct(eval(as.name(varsToPlot[1])), speciesGG, isYOY, riverGG, seasonGG, .keep_all = TRUE)

    ggOut <- ggplot( pGG, aes( eval(as.name(varsToPlot[1])), eval(as.name(depVar)), color=(speciesGG) ) ) +
      geom_line( size = 1.2 ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG)) +
      geom_hline(yintercept = 0) +
      #  ylim(-0.25,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      scale_y_continuous('Growth rate (mm/d)') +
   #   theme_publication() +
      facet_grid(riverGG~seasonGG)
  }

  if(length(varsToPlot) == 2) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]
    notPlot[3] <- all[!(all %in% varsToPlot)][3]


    #    pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0 ) %>%
    #      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), iter, isYOY, riverGG, species, season, .keep_all = TRUE)

    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0,
                        eval(as.name(notPlot[3])) == 0) %>%
      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), speciesGG, isYOY, riverGG, seasonGG, .keep_all = TRUE)

    pGG$speciesGroup <- paste0(pGG$speciesGG,pGG[[varsToPlot[2]]])

    ggOut <- ggplot(pGG, aes(eval(as.name(varsToPlot[1])), eval(as.name(depVar)), group = speciesGroup, name = varsToPlot[2])) +
      geom_line(aes( color = factor(speciesGG), alpha = (as.numeric(eval(as.name(varsToPlot[2])))+2)/3.5 ), size = 1.2 )  +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG)) +
      geom_hline(yintercept = 0) +
      #ylim(-0.33,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      scale_y_continuous('Growth rate (mm/d)') +
 #     theme_publication() +
      #    scale_fill_continuous(name=varsToPlot[2]) +
      facet_grid(riverGG~seasonGG)
  }
  print(ggOut)
  return(ggOut)
}


#'Plot observed/predicted and return RMSE
#'
#'@param residLimit THis is the lower limit for flagging a residual as an outlier
#'@return rmse and outliers, print text and plot graph
#'@export
#'
#'
getRMSE <- function(residLimit = 0.6){
  estLen <- array2df(dG[[1]]$q50$lengthExp, label.x = "estLen") %>% rename(rowNumber = d1)

  ddGIn <- ddG[[ii]][[2]]
  ddGIn$isEvalRow <- ifelse( ddGIn$lOcc == 0,TRUE, FALSE )

  ddGIn <- left_join( ddGIn, estLen, by = 'rowNumber') %>%
    mutate( resid = abs(estLen - lengthDATAOriginalLnStd),
            isOutlier = resid > residLimit )

  outlierFish1 <- ddGIn %>% filter(isOutlier) %>% select(tag)

  outlierFish <- ddGIn %>% filter(tag %in% outlierFish1$tag) %>%
    select(tag,season,observedLength,lengthDATAOriginalLnStd,estLen,resid,rowNumber,isOutlier,species)

  gg <- ggplot( ddGIn, aes( lengthDATAOriginalLnStd, estLen, color = isOutlier ) ) +
    geom_point( alpha = 0.2 ) +
    geom_abline(intercept = 0, slope = 1) +
    ylim(-2,4) +
    facet_wrap(~leftOut)
  print(gg)

  rmse <- ddGIn %>%
    mutate( resid = estLen - lengthDATAOriginalLnStd ) %>%
    group_by(leftOut) %>%
    summarise( rmse = sqrt( sum(resid^2,na.rm=T)/length(resid) ) )

  return(list(rmse = rmse, outliers = outlierFish))
}


#'Plot observed/predicted and return RMSE
#'
#'@param residLimit THis is the lower limit for flagging a residual as an outlier
#'@return rmse and outliers, print text and plot graph
#'@export
#'
#'
getRMSE_Nimble <- function(d,residLimit = 0.6, ii = 1){
  estLen <- array2df(d$q50$lengthExp, label.x = "estLen") %>%
    rename(rowNumber = d1) %>%
  #  filter( estLen != 0 ) #
    mutate( estLen = na_if(estLen, 0) )

  ddGIn <- ddG[[ii]][[2]]
  ddGIn$isEvalRow <- ifelse( ddGIn$lOcc == 0,TRUE, FALSE )

  ddGIn <- left_join( ddGIn, estLen, by = 'rowNumber') %>%
    #mutate( resid = abs(estLen - lengthDATAOriginalLnStd),
    mutate( resid = abs(estLen - observedLength),
            isOutlier = resid > residLimit )

  outlierFish1 <- ddGIn %>% filter(isOutlier) %>% dplyr::select(tag)

  outlierFish <- ddGIn %>% filter(tag %in% outlierFish1$tag) %>%
    dplyr::select(tag,season,sampleName,observedLength,lengthDATALnStd,observedLength,estLen,resid,rowNumber,isOutlier,species)

  gg <- ggplot( ddGIn, aes( observedLength, estLen, color = factor(year) ) ) +
        geom_point( alpha = 0.8 ) +
        geom_abline(intercept = 0, slope = 1) +
 #   ylim(-2,4) +
        facet_grid(isYOY+riverN~season+leftOut)
  print(gg)

  rmse <- ddGIn %>%
    mutate( resid = estLen - observedLength ) %>%
    group_by(leftOut) %>%
    summarise( rmse = sqrt( sum(resid^2,na.rm=T)/length(resid) ) )

  return(list(rmse = rmse, outliers = outlierFish))
}
