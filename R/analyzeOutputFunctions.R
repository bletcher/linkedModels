
#'Generate predictions from a set of iterations of a jags model run
#'
#'@param d a model run dataFrame
#'@param limits symetrical lower and upper limit for the predictions
#'@param nPoints the number of points between the limits. Use an odd # to ensure a value of 0
#'@param itersForPred, the iterations from the model for prediction
#'@return a data frame
#'@export

getPrediction <- function(d, limits = 2, nPoints = 5, itersForPred){

  # get grInt in df format
  grInt <- array2df(d$sims.list$grInt, levels = list(iter=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="int")

  # get grBeta in df format and merge in grInt
  grBeta1 <- array2df(d$sims.list$grBeta, levels = list(iter=NA,beta=NA,isYOY=c(0,1),species=species,season=1:nSeasons,river=riverOrderedIn), label.x="est")
  grBeta <- spread( grBeta1, key = beta, value = est, sep = "" )  %>%
    left_join( .,grInt )

  # prediction template
  x <- seq( -limits,limits,length.out = nPoints )

  # create inital template of standardized values to predict over
  predTemplate <- data.frame( len =   rep(x, each = nPoints ^ 3),
                              count = rep(x, each = nPoints ^ 2),
                              flow =  rep(x, each = nPoints ^ 1),
                              temp =      x
                            )


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
    mutate( predGr = int +
              beta1 * len +
              beta2 * count +
              beta3 * len^2 +
              beta4 * count^2 +
              beta5 * len * count +
              beta6 * temp +
              beta7 * flow +
              beta8 * temp^2 +
              beta9 * flow^2 +
              beta10 * temp * flow
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
