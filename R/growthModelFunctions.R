

#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  #grBetaOutside = array( runif( 11 *2*d$nSpecies*d$nSeasons*d$nRivers,-2,2),c( 11 ,2,d$nSpecies,d$nSeasons,d$nRivers))
  #grBetaOutside[ ,1,,2, ] <- 0

  nBetas <- 18
  nBetasSigma <- 4

  inits <- function(){
    list(
      grInt = array(rnorm(2*d$nSpecies*d$nSeasons*d$nRivers,0.5,0.25),c(2,d$nSpecies,d$nSeasons,d$nRivers)),
      #grBeta[x,1:2,1,1:4,1:4]
      grBeta = array(rnorm(nBetas*2*d$nSpecies*d$nSeasons*d$nRivers,0,0.1),c(nBetas,2,d$nSpecies,d$nSeasons,d$nRivers)),
      #grSigma[ yoy,spp,s,r ]
      grSigma = array(runif(2*d$nSpecies*d$nSeasons*d$nRivers,0,0.05),c(2,d$nSpecies,d$nSeasons,d$nRivers)),
      # sigmaBeta[ b,yoy,spp,s,r ]
      sigmaBeta = array(rnorm(nBetasSigma*2*d$nSpecies*d$nSeasons*d$nRivers,0,0.05),c(nBetasSigma,2,d$nSpecies,d$nSeasons,d$nRivers)),
      grIndRE = rnorm(d$nInd,0,0.1)
    )
  }

  # params <- c('grInt' ,'grIntMu','grIntSigma'
  #             ,'sigmaInt','sigmaIntMu','sigmaIntSigma'
  #             ,'grBeta','grBetaMu','grBetaSigma'
  #             , 'length','expectedGR'#,'expectedGR'
  #    #          , 'grIndRE','grIndREMean','grIndRETau'
  #             )

  params <- c('grInt', 'grBeta', 'grSigma','sigmaBeta', 'lengthExp', 'grIntMu', 'grIntSigma', 'grIndRE', 'sigmaInd', 'grBetaMu', 'grBetaSigma', 'sigmaBetaMu', 'sigmaBetaSigma' )

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel6.jags",
                n.chains = 3,
                n.adapt = 250, #1000
                n.iter = 250,
                n.burnin = 100,
                n.thin = 5,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

#'Add density and biomass data
#'
#'@param ddddGIn growth data
#'@param ddDIn results from the detection model
#'@param ddddGIn detection data
#'@param meanOrIterIn whether the growth model gets mean P from the detection model, or results from an iteration
#'@param sampleToUse if detection results are from an iteration, which iteration
#'@return a data frame with standardized counts, countPStd is standardized by species counts, countPAllSppStd is standardized by the sum of counts of all species in the analysis
#'@export
#'

addDensityData <- function( ddddGIn,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

  if ( meanOrIterIn == "mean") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    #ggplot(filter(den,!is.na(countP)),aes(year,countP, color = species)) + geom_point() + geom_line() + facet_grid(riverOrdered~season)
    #ggplot(filter(den,!is.na(countP)),aes(count,countP, color = species)) + geom_point()
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- meanOrIterIn
    print("in mean")
  }

  if ( meanOrIterIn == "iter") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- sampleToUse

    print("in iter ")
    print(c(meanOrIterIn,sampleToUse))
  }

  denForMerge <- den %>%
    dplyr::select(species, season, riverOrdered, year, countP, meanOrIter, iterIn) %>%
    filter( !is.na(countP) )

  # counts by species
  denForMergeSummary <- denForMerge %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( countPMean = mean( countP, na.rm = T ),
               countPSD = sd( countP, na.rm = T) )

  denForMerge2 <- left_join( denForMerge, denForMergeSummary ) %>%
    mutate( countPStd = ( countP - countPMean )/countPSD )

  ddddGIn <- left_join( ddddGIn,denForMerge2, by = (c("species", "year", "season", "riverOrdered")) )

  # counts are missing for samples with propSampled == 0. For now, fill in mean (0). Need this for when enc==0 and propSampled==0 in the gr model
  ddddGIn$countPStd <- ifelse( is.na(ddddGIn$countPStd), 0, ddddGIn$countPStd )

  #####
  # counts by species in separate columns
  countPStdBySpp <- denForMerge2 %>%
    select(species, season, riverOrdered, year, countPStd) %>%
    spread(species, countPStd, fill=0) %>%
    rename(countPStdBKT = bkt,countPStdBNT = bnt, countPStdATS = ats)
  ddddGIn <- left_join( ddddGIn,countPStdBySpp )

  # get mean masses by species
  massForMergeSummary <- ddddGIn %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( massMean = mean( observedWeight, na.rm = T ),
               massSD = sd( observedWeight, na.rm = T) )

  ddddGIn <- left_join( ddddGIn, massForMergeSummary ) %>%
    mutate( massStd = ( observedWeight - massMean )/massSD,
            biomassStd = massStd * countPStd )

  # biomasses are missing for samples with propSampled == 0. For now, fill in mean (0).
  ddddGIn$biomassStd <- ifelse( is.na(ddddGIn$biomassStd), 0, ddddGIn$biomassStd )

  ###################################################################################
  # counts across species (those in the analysis - e.g. species <- c("bkt", "bnt"))
  # sum over species
  denForMergeAllSpp <- denForMerge %>%
    group_by( season,riverOrdered,year ) %>%
    summarize( countPAllSpp = sum( countP, na.rm = T ) )

  denForMergeSummaryAllSpp <- denForMergeAllSpp %>%
    group_by( season,riverOrdered ) %>%
    summarize( countPAllSppMean = mean( countPAllSpp, na.rm = T ),
               countPAllSppSD = sd( countPAllSpp, na.rm = T) )

  denForMerge2AllSpp <- left_join( denForMergeAllSpp, denForMergeSummaryAllSpp ) %>%
    mutate( countPAllSppStd = ( countPAllSpp - countPAllSppMean )/countPAllSppSD )

  ddddGIn <- left_join( ddddGIn,denForMerge2AllSpp, by = (c("year", "season", "riverOrdered")) )

  # counts are missing for samples with propSampled == 0. For now, fill in mean (0).
  ddddGIn$countPAllSppStd <- ifelse( is.na(ddddGIn$countPAllSppStd), 0, ddddGIn$countPAllSppStd )

  # get mean masses across species
  massForMergeSummaryAllSpp <- ddddGIn %>%
    group_by( season,riverOrdered ) %>%
    summarize( massAllSppMean = mean( observedWeight, na.rm = T ),
               massAllSppSD = sd( observedWeight, na.rm = T) )

  ddddGIn <- left_join( ddddGIn, massForMergeSummaryAllSpp ) %>%
    mutate( massAllSppStd = ( observedWeight - massAllSppMean )/massAllSppSD,
            biomassAllSppStd = massAllSppStd * countPAllSppStd )

  # biomasses are missing for samples with propSampled == 0. For now, fill in mean (0).
  ddddGIn$biomassAllSppStd <- ifelse( is.na(ddddGIn$biomassAllSppStd), 0, ddddGIn$biomassAllSppStd )

  return(ddddGIn)

}

#'Turn observedLength values to NA for a percentage of the observations
#'
#'@param d a dataframe
#'@param runCrossValidation boolean for running cross validation or not
#'@return a data frame with observedLength set to NA for percentLeftOut observations
#'@export
#'
crossValidate <- function(d, runCrossValidationTF){
  d$observedLengthOriginal <- d$observedLength
  if ( runCrossValidationTF ) {
    propFOcc <- sum(d$fOcc) / nrow(d)
    d$leftOut <- ifelse( (runif(nrow(d)) < percentLeftOut/propFOcc/100) & (d$fOcc == 0), T, F ) # adjust percentLeftOut higher to acct for the fOcc that can't be NA
    d$observedLength <- ifelse( d$leftOut, NA, d$observedLength )
  }
  return(d)
}
