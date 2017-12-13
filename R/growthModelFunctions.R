

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
                n.adapt = 100, #1000
                n.iter = 200,
                n.burnin = 100,
                n.thin = 5,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

#'Adjust counts
#'
#'@param ddddGIn growth data
#'@param ddDIn results from the detection model
#'@param ddddGIn detection data
#'@param meanOrIterIn whether the growth model gets mean P from the detection model, or results from an iteration
#'@param sampleToUse if detection results are from an iteration, which iteration
#'@return a data frame with standardized counts, nAllFishBySpeciesPStd is standardized by species counts, nAllFishBySpeciesPAllSppStd is standardized by the sum of counts of all species in the analysis
#'@export
#'

adjustCounts <- function( ddddGIn,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

  if ( meanOrIterIn == "mean") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    #ggplot(filter(den,!is.na(nAllFishBySpeciesP)),aes(yearN,nAllFishBySpeciesP, color = factor(speciesN))) + geom_point() + geom_line() + facet_grid(riverN~seasonN)
    #ggplot(filter(den,!is.na(nAllFishBySpeciesP)),aes(count,nAllFishBySpeciesP, color = species)) + geom_point()
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
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesP, nAllFishP, massAllFishBySpecies, massAllFish, meanOrIter, iterIn) %>%
    filter( !is.na(nAllFishBySpeciesP) )

  # counts by species
  denForMergeSummaryBySpecies <- denForMerge %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( nAllFishBySpeciesPMean = mean( nAllFishBySpeciesP, na.rm = T ),
               nAllFishBySpeciesPSD =     sd( nAllFishBySpeciesP, na.rm = T),
               massAllFishBySpeciesMean = mean( massAllFishBySpecies, na.rm = T ),
               massAllFishBySpeciesSD =     sd( massAllFishBySpecies, na.rm = T))

  # counts for all species
  denForMergeSummary <- denForMerge %>%
    group_by( season,riverOrdered ) %>%
    summarize( nAllFishPMean = mean( nAllFishP, na.rm = T ),
               nAllFishPSD =     sd( nAllFishP, na.rm = T),
               massAllFishMean = mean( massAllFish, na.rm = T ),
               massAllFishSD =     sd( massAllFish, na.rm = T) )


  denForMerge2 <- left_join( denForMerge, denForMergeSummaryBySpecies ) %>%
                  left_join( denForMergeSummary ) %>%
                  mutate( nAllFishBySpeciesPStd = ( nAllFishBySpeciesP - nAllFishBySpeciesPMean )/nAllFishBySpeciesPSD,
                          nAllFishPStd =          ( nAllFishP -          nAllFishPMean )         /nAllFishPSD,
                          massAllFishBySpeciesStd = ( massAllFishBySpecies - massAllFishBySpeciesMean )/massAllFishBySpeciesSD,
                          massAllFishPStd =         ( massAllFish -          massAllFishMean )         /massAllFishSD)

  ddddGIn <- left_join( ddddGIn,denForMerge2, by = (c("species", "year", "season", "riverOrdered")) )

  # counts are missing for samples with propSampled == 0. For now, fill in mean (0). Need this for when enc==0 and propSampled==0 in the gr model
  #should this be a low # instead of 0???
  ddddGIn$nAllFishBySpeciesPStd <-   ifelse( is.na(ddddGIn$nAllFishBySpeciesPStd),   0, ddddGIn$nAllFishBySpeciesPStd )
  ddddGIn$massAllFishBySpeciesStd <- ifelse( is.na(ddddGIn$massAllFishBySpeciesStd), 0, ddddGIn$massAllFishBySpeciesStd )

  #####
  # counts by species in separate columns
  nAllFishBySpeciesPStdBySpp <- denForMerge2 %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesPStd) %>%
    spread(species, nAllFishBySpeciesPStd, fill=0) %>%
    rename(nAllFishBySpeciesPStdBKT = bkt,nAllFishBySpeciesPStdBNT = bnt, nAllFishBySpeciesPStdATS = ats)
  ddddGIn <- left_join( ddddGIn,nAllFishBySpeciesPStdBySpp )


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
