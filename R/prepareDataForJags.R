#'Get data ready to run the growth model
#'
#'@param d a data frame created using getCoreData()
#'@param modelType a character variable specifying the model type, ("detection', "growth"). This specifies the data that go into 'data'
#'@return A list for importing into the jags run. Run with run"modelType"Model()
#'@export

prepareDataForJags <- function(d,modelType){

  nInd <- n_distinct(d$tag, na.rm = T)
  nOcc <- n_distinct(d$sampleIndex, na.rm = T)
  minOcc <- min(d$sampleIndex)
  nRivers <- n_distinct(d$riverOrdered, na.rm = T)
  nSeasons <- n_distinct(d$season, na.rm = T)
  nYears <- n_distinct(d$year, na.rm = T)
  nSpecies <- n_distinct(d$species, na.rm = T)

  d <- d %>%
    group_by(tag) %>%
    mutate( minOcc = min(sampleNumber), fOcc = (sampleNumber == minOcc)*1,
            maxOcc = max(sampleNumber), lOcc = (sampleNumber == maxOcc)*1
    ) %>%
    ungroup()

  evalRows <- which( d$lOcc == 0 )
  firstObsRows <- which( d$fOcc == 1 )
  lastObsRows <- which( d$lOcc == 1 )

  nEvalRows <- length(evalRows)
  nFirstObsRows <- length(firstObsRows)
  nLastObsRows <- length(lastObsRows)

  nAllRows <- nEvalRows + nLastObsRows

  # add column for z initial values
  obsOcc <- d %>%
    dplyr::select( tag,sampleNumber,enc ) %>%
    filter(enc == 1) %>%
    group_by(tag) %>%
    mutate( minObsOcc = min(sampleNumber),
            maxObsOcc = max(sampleNumber)
    ) %>%
    dplyr::select( tag,minObsOcc,maxObsOcc ) %>%
    distinct()

  d <- left_join( d,obsOcc )
  d$zForInit <- ifelse( (d$sampleNumber > d$minObsOcc) & (d$sampleNumber <= d$maxObsOcc), 1, NA )

  load(file = paste0("./data/cutoffYOYInclSpring1DATA_",drainage,".RData"))
  cutoffYOYDATA <- cutoffYOYInclSpring1DATA # update as needed using getYOYCutoffs(cd,drainage)

  d$riverN <- as.numeric(d$riverOrdered)
  d$speciesN <- as.numeric(as.factor(d$species))

  means <- d %>%
    group_by( speciesN,season,riverN ) %>%
    summarize( meanLen = mean(observedLength, na.rm = T),
               sdLen = sd(observedLength, na.rm = T),
               sampleIntervalMean = mean(sampleInterval, na.rm = T)
             )

  propSampled <- getPropSampled(nSeasons,nRivers,nYears)

  ## standardize env variables
  meansEnv <- d %>%
    group_by( season,riverOrdered ) %>%
    summarize( meanTemperatureMean = mean(meanTemperature, na.rm = T),
               meanTemperatureSD = sd(meanTemperature, na.rm = T),
               meanFlowMean = mean(meanFlow, na.rm = T),
               meanFlowSD = sd(meanFlow, na.rm = T)

             )

  d <- left_join( d, meansEnv ) %>%
    mutate( tempStd = (meanTemperature - meanTemperatureMean)/meanTemperatureSD,
            flowStd = (meanFlow - meanFlowMean)/meanFlowSD
          )

  if ( modelType == "detection" ){
  data <- list( encDATA = d$enc,
                lengthDATA = d$observedLength,
                riverDATA = d$riverN,
                ind = d$tagIndex,
                nRivers = nRivers,
                nSpecies = nSpecies,
                #nInd = nInd,
                #nOcc = nOcc,
                #occ = d$sampleIndex - minOcc + 1,
                species = as.numeric(as.factor(d$species)),
                season = d$season,
                year = d$year - min(d$year) + 1,
                yearForCutoff = d$year - d$minYear + 1 + (d$minYear - 1997), # minYear is watershed-specific, -1997 because min year in cutoffYOYDATA is 1997
                nYears = nYears, #max(d$year) - min(d$year) + 1,
                nEvalRows = nEvalRows, evalRows = evalRows,
                nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
                nLastObsRows = nLastObsRows, lastObsRows = lastObsRows,
                nAllRows = nAllRows,
                nSeasons = nSeasons,
                lengthMean = array(means$meanLen, dim = c(nRivers,nSeasons,nSpecies)),
                lengthSD = array(means$sdLen, dim = c(nRivers,nSeasons,nSpecies)),
                sampleIntervalMean = array(means$sampleIntervalMean, dim = c(nRivers,nSeasons,nSpecies)),
                cutoffYOYDATA = cutoffYOYDATA,
                sampleInterval = d$sampleInterval,
                zForInit = d$zForInit, # z for firstObs gets set to zero in jags. Can't set values in inits for values assigned in jags
                propSampledDATA = propSampled$propSampledDATA,
                tempStd = d$tempStd,
                flowStd = d$flowStd
    )
  }
  if ( modelType == "growth" ){
    data <- list( encDATA = d$enc,
                  lengthDATA = d$observedLength,
                  riverDATA = d$riverN,
                  ind = d$tagIndex,
                  nRivers = nRivers,
                  nSpecies = nSpecies,
                  #nInd = nInd,
                  #nOcc = nOcc,
                  #occ = d$sampleIndex - minOcc + 1,
                  species = as.numeric(as.factor(d$species)),
                  season = d$season,
                  year = d$year - min(d$year) + 1,
                  yearForCutoff = d$year - d$minYear + 1 + (d$minYear - 1997), # minYear is watershed-specific, -1997 because min year in cutoffYOYDATA is 1997
                  nYears = nYears, #max(d$year) - min(d$year) + 1,
                  nEvalRows = nEvalRows, evalRows = evalRows,
                  nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
                  nLastObsRows = nLastObsRows, lastObsRows = lastObsRows,
                  nAllRows = nAllRows,
                  nSeasons = nSeasons,
                  lengthMean = array(means$meanLen, dim = c(nRivers,nSeasons,nSpecies)),
                  lengthSD = array(means$sdLen, dim = c(nRivers,nSeasons,nSpecies)),
                  sampleIntervalMean = array(means$sampleIntervalMean, dim = c(nRivers,nSeasons,nSpecies)),
                  cutoffYOYDATA = cutoffYOYDATA,
                  sampleInterval = d$sampleInterval,
                  zForInit = d$zForInit, # z for firstObs gets set to zero in jags. Can't set values in inits for values assigned in jags
                  propSampledDATA = propSampled$propSampledDATA,
                  countPStd = d$countPStd,
                  tempStd = d$tempStd,
                  flowStd = d$flowStd
    )
  }

  return(data)
}

