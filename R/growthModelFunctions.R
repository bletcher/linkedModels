

#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  inits <- function(){
    list(grBetaInt = array(rnorm(2*d$nSpecies*d$nSeasons*d$nRivers*d$nYears,0,2.25),c(2,d$nSpecies,d$nSeasons,d$nRivers,d$nYears)))
  }

  params <- c("grBetaInt","muGrBetaInt","sigmaGrBetaInt","grBeta","muGrBeta","grSigmaBeta","length")

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel.jags",
                n.chains = 3,
                n.adapt = 1000, #1000
                n.iter = 2000,
                n.burnin = 1000,
                n.thin = 4,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

#'Add density data
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
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

  denForMergeSummary <- denForMerge %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( countPMean = mean(countP, na.rm = T),
               countPSD = sd(countP, na.rm = T))

  denForMerge2 <- left_join( denForMerge, denForMergeSummary ) %>%
    mutate( countPStd = ( countP - countPMean )/countPSD )

  ddddGIn <- left_join( ddddGIn,denForMerge2, by = (c("species", "year", "season", "riverOrdered")) )
  return(ddddGIn)

}
