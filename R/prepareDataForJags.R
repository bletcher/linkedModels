#'Get data ready to run the detection or growth model
#'
#'@param d a data frame created using getCoreData()
#'@param modelType a character variable specifying the model type, ("detection', "growth"). This specifies the data that go into 'data'
#'@return A list for importing into the jags run. Run with run"modelType"Model()
#'@export

prepareDataForJags <- function(d,modelType){

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


  # Without error in lengthDATA (from Munch/sigourney), firstObs doesn't det estimated
  # There are a few fish with observations with fOcc=1 and NA for oberservedLength (e.g. 4 for west/bkt). Filter those fish out here
  noLenFOcc <- d %>% dplyr::select(tagIndex,observedLength,fOcc) %>% filter( (is.na(observedLength) & fOcc == 1) )
  d <- d %>% filter( !(tagIndex %in% noLenFOcc$tagIndex ) )
  print(paste("NA length for first observation", length(noLenFOcc$tagIndex)))
  print( noLenFOcc$tagIndex )

  d$leftOut <- F #placeholder
  d <- d %>% crossValidate( runCrossValidationTF )

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

  #load(file = paste0("./data/cutoffYOYInclSpring1DATA_",drainage,".RData"))
  #cutoffYOYDATA <- cutoffYOYInclSpring1DATA # update as needed using getYOYCutoffs(cd,drainage) in getAndPrepareDataWB.R

  d$riverN <- as.numeric(d$riverOrdered)
  d$speciesN <- as.numeric(as.factor(d$species))

  d$observedLengthLn <- log(d$observedLength)
  d$observedLengthOriginalLn <- log(d$observedLengthOriginal)

  means <- d %>%
    group_by( speciesN,season,riverN ) %>%
    summarize( meanLen = mean(observedLength, na.rm = T),
               sdLen = sd(observedLength, na.rm = T),
               meanLenOriginal = mean(observedLengthOriginal, na.rm = T),
               sdLenOriginal = sd(observedLengthOriginal, na.rm = T),
               meanLenLn = mean(observedLengthLn, na.rm = T),
               sdLenLn = sd(observedLengthLn, na.rm = T),
               meanLenOriginalLn = mean(observedLengthOriginalLn, na.rm = T),
               sdLenOriginalLn = sd(observedLengthOriginalLn, na.rm = T),
               sampleIntervalMean = mean(sampleInterval, na.rm = T)
             )

  d <- d %>%     # use original means so length and lengthOriginal have same std values
    mutate( lengthDATAStd = (observedLength -                 mean(observedLengthOriginal,na.rm = T)) / sd(observedLengthOriginal,na.rm = T),
            lengthDATAOriginalStd = (observedLengthOriginal - mean(observedLengthOriginal,na.rm = T)) / sd(observedLengthOriginal,na.rm = T),

            lengthDATALnStd = (observedLengthLn -                 mean(observedLengthOriginalLn,na.rm = T)) / sd(observedLengthOriginalLn,na.rm = T),
            lengthDATAOriginalLnStd = (observedLengthOriginalLn - mean(observedLengthOriginalLn,na.rm = T)) / sd(observedLengthOriginalLn,na.rm = T)
            )

  propSampled <- getPropSampled(nSeasons,nRivers,nYears)

  d <- addNPasses(d,drainage)
  d$nPasses <- ifelse( is.na(d$nPasses), 1, d$nPasses ) #nPasses gets NA when propSampled==0. Just set these to 1 so there are no NAs in the data. propSampled==0 takes care of thes in the jages code

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

  d$tagIndexJags <- as.numeric(as.factor(d$tagIndex)) #make sure have right tagIndices for subset of fish for Jags
  nInd <- n_distinct(d$tagIndexJags, na.rm = T)

  d$rowNumber <- 1:nrow(d)

  d$isYOYDATA <- ifelse( d$ageInSamples <= 3, 1, 2 )

  speciesByInd <- d %>% distinct(tagIndexJags,speciesN) %>% arrange(tagIndexJags)

  if ( modelType == "detection" ){
  data <- list( encDATA = d$enc,
                lengthDATA = d$lengthDATAStd,
                riverDATA = d$riverN,
                ind = d$tagIndexJags,
                nRivers = nRivers,
                nSpecies = nSpecies,
                nInd = nInd,
                #nOcc = nOcc,
                #occ = d$sampleIndex - minOcc + 1,
                species = as.numeric(factor(d$species, levels = c('bkt','bnt','ats'), ordered = T)), #this might screw up the indexing if a lower level is ignored in the species list   as.numeric(as.factor(d$species)),
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
              #  cutoffYOYDATA = cutoffYOYDATA,
                sampleInterval = d$sampleInterval,
                zForInit = d$zForInit, # z for firstObs gets set to zero in jags. Can't set values in inits for values assigned in jags
                propSampledDATA = propSampled$propSampledDATA,
                tempStd = d$tempStd,
                flowStd = d$flowStd,
                nPasses = d$nPasses
    )
  }


  #fill NAs for testing
 # d$lInterp =  na.approx(d$observedLength)
#  d$lenInit <- ifelse( is.na(d$observedLength), d$lInterp, NA )

#  d$initialIsYOY <- ifelse( is.na(d$lenInit), NA, ifelse( d$lInterp > 90, 2, 1 ) )
    # div <- 10
  # sep <- round(nrow(d)/div)
  # keep <- 1:sep
  # interp <- (sep + 1):nrow(d)
  # d$lInterp = c( d$lInterp[keep], na.approx(d$lInterp[interp]) )


  if ( modelType == "growth" ){
    data <- list( encDATA = d$enc,
                  lengthDATA = d$lengthDATALnStd, #d$lengthDATAStd,
                  riverDATA = d$riverN,
                  ind = d$tagIndexJags,
                  nRivers = nRivers,
                  nSpecies = nSpecies,
                  nInd = nInd,
                  #nOcc = nOcc,
                  #occ = d$sampleIndex - minOcc + 1,
                  species = as.numeric(factor(d$species, levels = c('bkt','bnt','ats'), ordered = T)), #this might screw up the indexing if a lower level is ignored in the species list   as.numeric(as.factor(d$species)),
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
               #   cutoffYOYDATA = cutoffYOYDATA,
                  sampleInterval = d$sampleInterval,
                  zForInit = d$zForInit, # z for firstObs gets set to zero in jags. Can't set values in inits for values assigned in jags
                  propSampledDATA = propSampled$propSampledDATA,
                  countPStd = d$countPStd,
                  tempStd = d$tempStd,
                  flowStd = d$flowStd,
   #               biomassDeltaAllSpp = d$meanBiomassAllSppStdDelta,
    #              biomassDelta = d$meanBiomassStdDelta,
           #       logitPhiStd = d$logitPhiStd,
  #                lForInit = d$lInterp,
                  isYOYDATA = d$isYOYDATA,
                  speciesByInd = speciesByInd$speciesN
    )
  }

  return(list(data,d))
}

