

#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  inits <- function(){
    list(grBetaInt = array(rnorm(2*dddG$nSpecies*dddG$nSeasons*dddG$nRivers*dddG$nYears,0,2.25),c(2,dddG$nSpecies,dddG$nSeasons,dddG$nRivers,dddG$nYears)))
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

addDensityData <- function( ddddG,ddD,ddddD,meanOrIter,sampleToUse ){

  if ( meanOrIter == "mean") {
    den <- getDensities(ddddD, ddD, meanOrIter = meanOrIter, sampleToUse = sampleToUse )
    #ggplot(filter(den,!is.na(countP)),aes(year,countP, color = species)) + geom_point() + geom_line() + facet_grid(riverOrdered~season)
    #ggplot(filter(den,!is.na(countP)),aes(count,countP, color = species)) + geom_point()
  }

  if ( meanOrIter == "iter") {
    den <- getDensities(ddddD, ddD, meanOrIter = meanOrIter, sampleToUse = sampleToUse )
  }

  denForMerge <- den %>%
    dplyr::select(species, season, riverOrdered, year, countP) %>%
    filter( !is.na(countP) )

  ddddG <- left_join( ddddG,denForMerge )
  return(ddddG)

}
