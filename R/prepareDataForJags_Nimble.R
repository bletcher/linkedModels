#'Get data ready to run the detection or growth model
#'
#'@param d a data frame created using getCoreData()
#'@param modelType a character variable specifying the model type, ("detection', "growth"). This specifies the data that go into 'data'
#'@return A list for importing into the jags run. Run with run"modelType"Model()
#'@export

prepareDataForJags_Nimble <- function(d,modelType){

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
  # There are a few fish with observations with fOcc=1 and NA for observedLength (e.g. 4 for west/bkt). Filter those fish out here
  noLenFOcc <- d %>% dplyr::select(tagIndex,observedLength,fOcc) %>% filter( (is.na(observedLength) & fOcc == 1) )
  d <- d %>% filter( !(tagIndex %in% noLenFOcc$tagIndex ) )
  print(paste("NA length for first observation, # =", length(noLenFOcc$tagIndex)))
  print( noLenFOcc$tagIndex )

  d$leftOut <- FALSE #placeholder
  d <- d %>% crossValidate( runCrossValidation_TF )

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
  d$speciesN <- as.numeric(factor(d$species, levels = c('bkt','bnt','ats'), ordered = T)) #as.numeric(as.factor(d$species))

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

  d$isYOYDATA <- ifelse( d$ageInSamples <= 3, 1, 2 )

##########################################################
# means within a sample - use for relative size at sample
##########################################################
  meansByYOY_Year <- d %>%
    group_by( isYOYDATA,speciesN,season,year,riverN ) %>%
    summarize( meanLen_ByYOY_Year = mean(observedLength, na.rm = T),
               sdLen_ByYOY_Year = sd(observedLength, na.rm = T),
               meanLenOriginal_ByYOY_Year = mean(observedLengthOriginal, na.rm = T),
               sdLenOriginal_ByYOY_Year = sd(observedLengthOriginal, na.rm = T),
               meanLenLn_ByYOY_Year = mean(observedLengthLn, na.rm = T),
               sdLenLn_ByYOY_Year = sd(observedLengthLn, na.rm = T),
               meanLenOriginalLn_ByYOY_Year = mean(observedLengthOriginalLn, na.rm = T),
               sdLenOriginalLn_ByYOY_Year = sd(observedLengthOriginalLn, na.rm = T),
               sampleIntervalMean_ByYOY_Year = mean(sampleInterval, na.rm = T)
    )
  d <- left_join( d, meansByYOY_Year )

  d <- d %>%     #  # use original means so length and lengthOriginal have same std values
    mutate( lengthDATAStd_ByYOY_Year = (observedLength -                 meanLenOriginal_ByYOY_Year) / sdLenOriginal_ByYOY_Year,
            lengthDATAOriginalStd_ByYOY_Year = (observedLengthOriginal - meanLenOriginal_ByYOY_Year) / sdLenOriginal_ByYOY_Year,

            lengthDATALnStd_ByYOY_Year = (observedLengthLn -                 meanLenOriginalLn_ByYOY_Year) / sdLenOriginalLn_ByYOY_Year,
            lengthDATAOriginalLnStd_ByYOY_Year = (observedLengthOriginalLn - meanLenOriginalLn_ByYOY_Year) / sdLenOriginalLn_ByYOY_Year
    )
#############################################
#############################################

  #done in addEnvironmental()
  #propSampled <- getPropSampled(nSeasons,nRivers,nYears,min(d$year))

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

  d$proportionSampled <- ifelse( is.na(d$proportionSampled), 1, d$proportionSampled )

  speciesByInd <- d %>% distinct(tagIndexJags,speciesN) %>% arrange(tagIndexJags)

  print(paste0("Number of input rows = ",nrow(d),", Number of encounters = ",sum(d$enc,na.rm=T)))

  # get ind that were only observed for 1 gr interval - for grIndRE
  d <- d %>%
    group_by(tag) %>%
    mutate( onlyOneGrInterval = (maxOcc == (minOcc + 1)) ) %>%
    ungroup()

  onlyOne <- d %>% group_by(tag) %>%
    mutate( onlyOneGrInterval = (maxOcc == (minOcc + 1)) ) %>%
    summarize(onlyOne_TF = unique(onlyOneGrInterval)) %>%
    ungroup()

  moreThanOneInterval <- which(!onlyOne$onlyOne_TF)
  exactlyOneInterval <- which(  onlyOne$onlyOne_TF)

  nMoreThanOneInterval <- length(moreThanOneInterval)
  nExactlyOneInterval <- length(exactlyOneInterval)

  ##############
  # from linearGR_ModelsMap.Rmd

  #d$grLength <- ifelse(d$species == 2 & d$river == 3, NA, d$grLength)
  #d$grLength <- ifelse(d$species %in% c('2','3') & d$isYOY == 1 & d$season == "2", NA, d$grLength)
  #d$grLength <- ifelse(d$grLength < -0.025, NA, d$grLength)


  ##############



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
                  propSampledDATA = d$proportionSampled, # propSampled$propSampledDATA,
                  tempStd = d$tempStd,
                  flowStd = d$flowStd,
                  nPasses = d$nPasses,
                  isYOYDATA = d$isYOYDATA
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

    # d <- d %>%
    #   mutate(
    #     cBKTStd = as.numeric(scale(nAllFishBySpeciesPBKT)),
    #     cBNTStd = as.numeric(scale(nAllFishBySpeciesPBNT)),
    #     cATSStd = as.numeric(scale(nAllFishBySpeciesPATS))
    #     #   cBKTStd2 = cBKTStd^2,
    #     #    cBNTStd2 = cBNTStd^2,
    #     #    cATSStd2 = cATSStd^2
    #   )

    ## standardize abundances by river
    nMeans <- d %>%
      group_by( riverOrdered ) %>%
      summarize( meanBKT = mean(nAllFishBySpeciesPBKT, na.rm = T),
                 meanBNT = mean(nAllFishBySpeciesPBNT, na.rm = T),
                 meanATS = mean(nAllFishBySpeciesPATS, na.rm = T),
                 sdBKT = sd(nAllFishBySpeciesPBKT, na.rm = T),
                 sdBNT = sd(nAllFishBySpeciesPBNT, na.rm = T),
                 sdATS = sd(nAllFishBySpeciesPATS, na.rm = T)

      )
#print(nMeans)
    d <- left_join( d, nMeans ) %>%
      mutate( cBKTStd = (nAllFishBySpeciesPBKT - meanBKT)/sdBKT,
              cBNTStd = (nAllFishBySpeciesPBNT - meanBNT)/sdBNT,
              cATSStd = (nAllFishBySpeciesPATS - meanATS)/sdATS
      )




    #To fill in missing obs (no fish in stream) with mean value
    d$cBKTStd <- ifelse( is.na(d$cBKTStd), 0, d$cBKTStd )
    d$cBNTStd <- ifelse( is.na(d$cBNTStd), 0, d$cBNTStd )
    d$cATSStd <- ifelse( is.na(d$cATSStd), 0, d$cATSStd )

    d$BKT01DATA <- 1
    d$BNT01DATA <- ifelse( d$river %in% c('west brook','wb jimmy'), 1,0 )

    d$ATS01DATA1 <- ifelse( d$river == 'west brook', 1,0 )
    d$ATS01DATA2 <- ifelse( d$sampleIndex <= 46, 1,0 )
    d$ATS01DATA <- d$ATS01DATA1 * d$ATS01DATA2

    data <- list( encDATA = d$enc,
                  lengthDATA = d$observedLength, #d$lengthDATAStd, #d$lengthDATALnStd,
                  lengthDATAStd = d$lengthDATALnStd,
                  lengthDATAStd_ByYOY_Year = d$lengthDATALnStd_ByYOY_Year,
                  riverDATA = d$riverN,
                  ind = d$tagIndexJags,
                  sampleNumber = d$sampleIndex,
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
                  nMoreThanOneInterval = nMoreThanOneInterval, moreThanOneInterval = moreThanOneInterval,
                  nExactlyOneInterval = nExactlyOneInterval, exactlyOneInterval = exactlyOneInterval,
                  nAllRows = nAllRows,
                  nSeasons = nSeasons,
                  lengthMean = mean(d$observedLength, na.rm = TRUE), #array(means$meanLen, dim = c(nRivers,nSeasons,nSpecies)),
                  lengthSD = sd(d$observedLength, na.rm = TRUE), #array(means$sdLen, dim = c(nRivers,nSeasons,nSpecies)),
                  sampleIntervalMean = array(means$sampleIntervalMean, dim = c(nRivers,nSeasons,nSpecies)),
                  #   cutoffYOYDATA = cutoffYOYDATA,
                  sampleInterval = d$sampleInterval,
                  zForInit = d$zForInit, # z for firstObs gets set to zero in jags. Can't set values in inits for values assigned in jags
                  propSampledDATA = d$proportionSampled, # propSampled$propSampledDATA,
                  countPBySpeciesStd = d$nAllFishBySpeciesPStd_Mean,
                  # countPStdBKT_yoy1 = d$nAllFishBySpeciesPStdBKT_yoy1,
                  # countPStdBKT_yoy2 = d$nAllFishBySpeciesPStdBKT_yoy2,
                  # countPStdBNT_yoy1 = d$nAllFishBySpeciesPStdBNT_yoy1,
                  # countPStdBNT_yoy2 = d$nAllFishBySpeciesPStdBNT_yoy2,
                  # countPStdATS_yoy1 = d$nAllFishBySpeciesPStdATS_yoy1,
                  # countPStdATS_yoy2 = d$nAllFishBySpeciesPStdATS_yoy2,
                  countPStdBKT = d$nAllFishBySpeciesPStdBKT,
                  countPStdBNT = d$nAllFishBySpeciesPStdBNT,
                  countPStdATS = d$nAllFishBySpeciesPStdATS,

                  cBKTStd = d$cBKTStd,
                  cBNTStd = d$cBNTStd,
                  cATSStd = d$cATSStd,

                  BKT01DATA = d$BKT01DATA,
                  BNT01DATA = d$BNT01DATA,
                  ATS01DATA = d$ATS01DATA,

                  countPAllSppStd = d$nAllFishPStd_Mean,
                  countPAllSpp = d$nAllFishP_Mean,
                  tempStd = d$tempStd,
                  flowStd = d$flowStd,
                  tempStd2 = d$tempStd^2,
                  flowStd2 = d$flowStd^2,
                  #               biomassDeltaAllSpp = d$meanBiomassAllSppStdDelta,
                  #              biomassDelta = d$meanBiomassStdDelta,
                  #       logitPhiStd = d$logitPhiStd,
                  #                lForInit = d$lInterp,
                  isYOYDATA = d$isYOYDATA,
                  speciesByInd = speciesByInd$speciesN,
                  grNotUse = d$grLength
    )
  }

  return(list(data,d))
}

