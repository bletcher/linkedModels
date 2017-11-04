
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
#  grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")

  # get grBeta in df format and merge in grInt
  grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" )  #%>%
 #   left_join( .,grInt )

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

  predTemplate <- data.frame( len =   lenData,
                              count = countData,
                         #     phi =   phiData,
                              flow =  flowData,
                              temp =  tempData
                            )

  # predTemplate <- data.frame( len =   rep(x, each = nPoints ^ 4),
  #                             count = rep(x, each = nPoints ^ 3),
  #                             phi =  rep(x, each = nPoints ^ 2),
  #                             flow =  rep(x, each = nPoints ^ 1),
  #                             temp =      x
  #
  #                             #len =  rep(x, each = nPoints ^ 1),
  #                             #biomass =      x
  #)


    # expand predTemplate across grBeta rows for a given set of iterations
    grBetaIter <- grBeta %>% filter( iter %in% itersForPred)

    # repeat predTemplate nrow( grBetaIter ) times
    predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
    # repeat each row of grbetaIter nrow(predTemplate) times
    grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
    # Put them together
    preds <- cbind( predTemplateLong,grBetaLong )

  # This model structure needs to match that in grModel.jags
  preds <- preds %>%
    mutate( predGr =
              #int +

              beta1 * len +
              beta2 * count +
           #   beta3 * phi +
              beta3 * temp +
              beta4 * flow +
              beta5 * len^2 +
              beta6 * count^2 +
          #    beta8 * phi^2 +
              beta7 * temp^2 +
              beta8 * flow^2 +
              beta9 * len * count +
           #   beta12 * len * phi +
              beta10 * temp * flow +
              beta11 * temp * flow * len

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

getPredictionSigma <- function(d, limits = 2, nPoints = 5, itersForPred){

  # get grInt in df format
  grInt <- array2df(d$sims.list$sigmaInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")

  # get grBeta in df format and merge in grInt
  grBeta1 <- array2df(d$sims.list$sigmaBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" )  %>%
    left_join( .,grInt )

  # prediction template
  x <- seq( -limits,limits,length.out = nPoints )

  # create inital template of standardized values to predict over
  predTemplate <- data.frame( #len =   rep(x, each = nPoints ^ 3),
                              count = rep(x, each = nPoints ^ 2),
                              flow =  rep(x, each = nPoints ^ 1),
                              temp =      x
  )


  # expand predTemplate across grBeta rows for a given set of iterations
  grBetaIter <- grBeta %>% filter( iter %in% itersForPred )

  # repeat predTemplate nrow( grBetaIter ) times
  predTemplateLong <- do.call("rbind", replicate( nrow( grBetaIter ), predTemplate, simplify = FALSE) )
  # repeat each row of grbetaIter nrow(predTemplate) times
  grBetaLong <- grBetaIter[rep(seq_len(nrow(grBetaIter)), each = nrow(predTemplate)),]
  # Put them together
  preds <- cbind( predTemplateLong,grBetaLong )

  # This model structure needs to match that in grModel.jags
  preds <- preds %>%
    mutate( predGr = int +

              beta1 * count +

              beta2 * count^2 +

              beta3 * temp +
              beta4 * flow +
              beta5 * temp^2 +
              beta6 * flow^2 +
              beta7 * temp * flow
    )

  return(preds)
}
