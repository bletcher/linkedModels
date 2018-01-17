
#'Generate predictions from a set of iterations of a jags model run
#'
#'@param d a model run dataFrame
#'@param limits symetrical lower and upper limit for the predictions
#'@param nPoints the number of points between the limits. Use an odd # to ensure a value of 0
#'@param itersForPred, the iterations from the model for prediction
#'@param varsToEstimate which variables to use for estimates (len,count,phi,temp,flow)
#'@return a data frame
#'@export

getPrediction <- function(d, limits = 2, nPoints = 5, itersForPred, varsToEstimate){

  # get grInt in df format
  #grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
  grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),season=1:nSeasons,river=riverOrderedIn), label.x="int")

  # get grBeta in df format and merge in grInt
 # grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),season=1:nSeasons,river=riverOrderedIn), label.x="est")

  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" ) %>%
    left_join( .,grInt )

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
    mutate( predGr =
              int +
              beta1 * len +
              beta2 * count +
              beta3 * temp +
              beta4 * flow +

          #    beta4 * len^2 +
              beta5 * count^2 +
              beta6 * temp^2 +
              beta7 * flow^2 +

              beta8 * temp * flow +
              beta9 * temp * count +
              beta10 * count * flow +

      #        beta12 * len * count +
      #        beta13 * len * flow +
      #        beta14 * len * temp +

              beta11 * count * temp * flow
       #       beta16 * len * temp * flow +
      #        beta17 * count * len * flow +
       #       beta18 * count * temp * len

              # int +
              # beta1 * len +
              # beta2 * count +
              # beta3 * temp +
              # beta4 * flow +
              #
              # beta5 * len^2 +
              # beta6 * count^2 +
              # beta7 * temp^2 +
              # beta8 * flow^2 +
              #
              # beta9 * temp * flow +
              # beta10 * temp * count +
              # beta11 * count * flow +
              #
              # beta12 * len * count +
              # beta13 * len * flow +
              # beta14 * len * temp +
              #
              # beta15 * count * temp * flow +
              # beta16 * len * temp * flow +
              # beta17 * count * len * flow +
              # beta18 * count * temp * len

    )

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

getPredictionSigma <- function(d, limits = 2, nPoints = 5, itersForPred, varsToEstimate){

  # get grInt in df format
  #grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")
  grInt <- array2df(d$sims.list$sigmaInt, levels = list(iter=NA,isYOY=c(0,1),season=1:nSeasons,river=riverOrderedIn), label.x="int")

  # get grBeta in df format and merge in grInt
  # grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta1 <- array2df(d$sims.list$sigmaBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),season=1:nSeasons,river=riverOrderedIn), label.x="est")

  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" ) %>%
    left_join( .,grInt )

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
    mutate( predGr = int +

              beta1 * count +
              beta2 * temp +
              beta3 * flow +
              beta4 * count * temp +
              beta5 * temp * flow +
              beta6 * count * flow

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

plotPred <- function(p, varsToPlot, isYOYGG, speciesGG) {
  all = c('len','temp','flow','count')


  if(length(varsToPlot) == 1) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]
    notPlot[3] <- all[!(all %in% varsToPlot)][3]

 #   pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0 ) %>%
#                 distinct(eval(as.name(varsToPlot[1])), iter, isYOY, river, species, season, .keep_all = TRUE)
    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0, eval(as.name(notPlot[3])) == 0  ) %>%
      distinct(eval(as.name(varsToPlot[1])), iter, isYOY, river, season, .keep_all = TRUE)

    ggOut <- ggplot(pGG, aes(eval(as.name(varsToPlot[1])),predGr, group = iter,color=(iter))) +
      geom_line( alpha=0.25 ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG,', species = ',speciesGG)) +
      geom_hline(yintercept = 0) +
  #    ylim(-0.25,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      facet_grid(river~season)
  }

  if(length(varsToPlot) == 2) {
    notPlot <- NA
    notPlot[1] <- all[!(all %in% varsToPlot)][1]
    notPlot[2] <- all[!(all %in% varsToPlot)][2]

#    pGG <- p %>% filter(isYOY == isYOYGG, species == speciesGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0 ) %>%
#      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), iter, isYOY, river, species, season, .keep_all = TRUE)

    pGG <- p %>% filter(isYOY == isYOYGG, eval(as.name(notPlot[1])) == 0, eval(as.name(notPlot[2])) == 0 ) %>%
      distinct(eval(as.name(varsToPlot[1])), eval(as.name(varsToPlot[2])), iter, isYOY, river, season, .keep_all = TRUE)

    pGG$iterGroup <- paste0(pGG$iter,pGG[[varsToPlot[2]]])

    ggOut <- ggplot(pGG, aes(eval(as.name(varsToPlot[1])), predGr, group = iterGroup, name = varsToPlot[2])) +
      geom_line(aes( color = (eval(as.name(varsToPlot[2]))) ), alpha = 0.25 ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      ggtitle(paste0('isYOY = ',isYOYGG,', species = ',speciesGG)) +
      geom_hline(yintercept = 0) +
  #    ylim(-0.33,1.33) +
      scale_x_continuous(varsToPlot[1]) +
      #    scale_fill_continuous(name=varsToPlot[2]) +
      facet_grid(river~season)
  }
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
