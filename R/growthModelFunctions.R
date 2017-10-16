

#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  inits <- function(){
#    list(grBetaInt = array(rnorm(2*d$nSpecies*d$nSeasons*d$nRivers*d$nYears,0,2.25),c(2,d$nSpecies,d$nSeasons,d$nRivers,d$nYears)))
    list(grInt = array(rnorm(2*d$nSpecies*d$nSeasons*d$nRivers,0,2.25),c(2,d$nSpecies,d$nSeasons,d$nRivers)))
     }

#  params <- c("grBetaInt","muGrBetaInt","sigmaGrBetaInt","grBeta","muGrBeta","grSigmaBeta","sigmaGrSigmaBeta")
  params <- c('grInt','sigmaInt','grBeta',
           #   'sigmaBeta',
              'grIntMu','grIntSigma','sigmaIntMu','sigmaIntSigma','grBetaMu','grBetaSigma'
        #    ,  'sigmaBetaMu',
        #    'sigmaBetaSigma'
 #              , 'grIndRE'
              )

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel2.jags",
                n.chains = 3,
                n.adapt = 1000, #1000
                n.iter = 200,
                n.burnin = 100,
                n.thin = 2,
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

  # get mean masses by species
  massForMergeSummary <- ddddGIn %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( massMean = mean( observedWeight, na.rm = T ),
               massSD = sd( observedWeight, na.rm = T) )

  ddddGIn <- left_join( ddddGIn, massForMergeSummary ) %>%
    mutate( massStd = ( observedWeight - massMean )/massSD,
            biomassStd = massStd * countPStd )

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

  # get mean masses across species
  massForMergeSummaryAllSpp <- ddddGIn %>%
    group_by( season,riverOrdered ) %>%
    summarize( massAllSppMean = mean( observedWeight, na.rm = T ),
               massAllSppSD = sd( observedWeight, na.rm = T) )

  ddddGIn <- left_join( ddddGIn, massForMergeSummaryAllSpp ) %>%
    mutate( massAllSppStd = ( observedWeight - massAllSppMean )/massAllSppSD,
            biomassAllSppStd = massAllSppStd * countPAllSppStd )

  return(ddddGIn)

}

