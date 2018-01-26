#'Run the growth model using Nimble
#'
#'@param d, input dataframe
#'@param mcmcInfo, a list containing run info
#'@return a data frame
#'@export

runGrowthModel_Nimble <- function(d,mcmcInfo){

  code <- nimbleCode({
    ##
    for( i in 1:nEvalRows ) {
      ##
      lengthDATA[ evalRows[i] + 1 ] ~ dnorm( lengthDATA[ evalRows[i] ] + gr[ evalRows[i] ], sd = expectedGRSigma[ evalRows[i] ] )
      ##
      gr[ evalRows[i] ] <-
        grInt[         isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
        grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] +
        grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] +
        grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
        grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
        ## grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]]^2 +
        grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]]^2 +
        grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]]^2 +
        grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]]^2 +
        grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
        grBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * countPAllSppStd[evalRows[i]] +
        grBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * flowStd[evalRows[i]] +
        grBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
        ## grBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * countPAllSppStd[evalRows[i]] +
        ## grBeta[13, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * flowStd[evalRows[i]] +
        ## grBeta[14, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * tempStd[evalRows[i]] +
        ## grBeta[16, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
        ## grBeta[17, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * lengthDATA[evalRows[i]] * flowStd[evalRows[i]] +
        ## grBeta[18, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] * lengthDATA[evalRows[i]] +
        grIndRE[ ind[evalRows[i]] ]
      ##
      log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
        sigmaInt[ isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
        sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] +
        sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
        sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
        sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] +
        sigmaBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
        sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * flowStd[evalRows[i]]
    }
    ##
    for( s in 1:nSeasons ) {
      grIntMu[ s ] ~ dnorm(0,0.001)
      grIntSigma[ s ] ~ dunif(0,10)
      ##
      sigmaIntMu[ s ] ~ dnorm(0,0.001)
      sigmaIntSigma[ s ] ~ dunif(0,10)
      ##
      for( b in 1:11 ) {  ############# Beta ##############################################################################
        grBetaMu[ b,s ] ~ dnorm(0,0.001)
        grBetaSigma[ b,s ] ~ dunif(0,100)
      }
      ##
      for( b in 1:6 ) {  ############# Beta ##############################################################################
        sigmaBetaMu[ b,s ] ~ dnorm(0,0.001)
        sigmaBetaSigma[ b,s ] ~ dunif(0,100)
      }
      ##
      for( yoy in 1:2 ) {
        for( r in 1:nRivers ) {
          grInt[ yoy,s,r ] ~ dnorm( grIntMu[ s ], sd = grIntSigma[ s ] )
          sigmaInt[ yoy,s,r ] ~ dnorm( sigmaIntMu[ s ], sd = sigmaIntSigma[ s ] )    #T(0,) #dunif(0,100)
          ##
          for( b in 1:11 ) {  ############# Beta ##################################################################################
            grBeta[ b,yoy,s,r ] ~ dnorm( grBetaMu[ b,s ], sd = grBetaSigma[ b,s ] )
          }
          for( b in 1:6 ) {  ############# Beta ##################################################################################
            sigmaBeta[ b,yoy,s,r ] ~ dnorm( sigmaBetaMu[ b,s ], sd = sigmaBetaSigma[ b,s ] )
          }
        }
      }
    }
    ##
    for(ii in 1:nEvalRows) {
      lengthExp[ evalRows[ii] + 1 ] <- lengthDATA[ evalRows[ii] ] + gr[ evalRows[ii] ]
    }
    ##
    # including this messes up predictions
    # for(ii in 1:nFirstObsRows) {
    #   lengthExp[ firstObsRows[ii]  ] <- lengthDATA[ firstObsRows[ii] ]
    # }
    ##
    ## individual random effect on gr
    for(iii in 1:nInd) {
      grIndRE[ iii ] ~ dnorm( 0, sd = sigmaIndRE )
    }
    ##
    sigmaIndRE ~ dunif(0,100)
  })

  ##tail(ddG[[1]][[2]])
  ##str(d)
  ##names(d)
  ##"encDATA"
  ##"species"
  ##"year"               "yearForCutoff"      "nYears"
  ##"evalRows"           "nFirstObsRows"
  ##"firstObsRows"       "nLastObsRows"       "lastObsRows"
  ##"nAllRows"           "lengthMean"
  ##"lengthSD"           "sampleIntervalMean" "sampleInterval"
  ##"zForInit"           "propSampledDATA"    "countPBySpeciesStd"
  ##"countPStdBKT"       "countPStdBNT"       "countPStdATS"
  ##"speciesByInd"       "grNotUse"

  constants <- list(riverDATA = d$riverDATA,
                    nRivers = d$nRivers,
                    nSpecies = d$nSpecies,
                    ind = d$ind,
                    nInd = d$nInd,
                    season = d$season,
                    nEvalRows = d$nEvalRows,
                    evalRows = d$evalRows,
                    nSeasons = d$nSeasons,
                    countPAllSppStd = d$countPAllSppStd,
                    isYOYDATA = d$isYOYDATA,
                    tempStd = d$tempStd,
                    flowStd = d$flowStd)
  ##
  #data <- list(lengthDATA = d$lengthDATA)
  #data <- list(lengthDATA = d$lengthDATA[1:33518])
  data <- list(lengthDATA = d$lengthDATA[1:34701])
  ##
  nBetas <- 11
  nBetasSigma <- 6
  ##
  inits <- list(grInt = array(rnorm(2 * constants$nSeasons * constants$nRivers, 0.5, 0.25), c(2, constants$nSeasons, constants$nRivers)),
                grBeta = array(rnorm(nBetas * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetas, 2, constants$nSeasons, constants$nRivers)),
                ##grSigma = array(runif(2 * constants$nSeasons * constants$nRivers, 0, 0.05), c(2, constants$nSeasons, constants$nRivers)),
                sigmaBeta = array(rnorm(nBetasSigma * 2 * constants$nSeasons * constants$nRivers, 0, 0.05), c(nBetasSigma, 2, constants$nSeasons, constants$nRivers)),
                grIndRE = rnorm(constants$nInd, 0, 0.1),
                sigmaIndRE = 1,
                grIntMu = rep(0, constants$nSeasons),
                grIntSigma = rep(1, constants$nSeasons),
                grBetaMu = array(0, c(nBetas, constants$nSeasons)),
                grBetaSigma = array(1, c(nBetas, constants$nSeasons)),
                sigmaInt = array(1, c(2, constants$nSeasons, constants$nRivers)),
                sigmaIntMu = rep(0, constants$nSeasons),
                sigmaIntSigma = rep(1, constants$nSeasons),
                sigmaBetaMu = array(0, c(nBetasSigma, constants$nSeasons)),
                sigmaBetaSigma = array(1, c(nBetasSigma, constants$nSeasons))
  )
  ##
  params <- c('grInt', 'grIntMu', 'grIntSigma', 'grBeta', 'grBetaMu',
              'grBetaSigma', 'sigmaInt', 'sigmaIntMu', 'sigmaIntSigma',
              'sigmaBeta', 'sigmaBetaMu', 'sigmaBetaSigma', 'grIndRE', 'sigmaIndRE',
              'lengthExp')



  # Rmodel <- nimbleModel(code, constants, data, inits)
  #
  # #  Rmodel$lengthDATA <- zoo::na.approx(data$lengthDATA[1:33518]) ## length(data$lengthDATA) = 33519, last obs is a fish with a single obs - doesn't get an evalRow
  # Rmodel$lengthDATA <- zoo::na.approx(data$lengthDATA[1:34701]) ## length(data$lengthDATA) = 33519, last obs is a fish with a single obs - doesn't get an evalRow
  #   table(is.na(Rmodel$lengthDATA))
  #
  # #system.time(lp <- Rmodel$calculate())
  # #lp
  #
  # ##Rmodel$getVarNames(nodes = Rmodel$getNodeNames(stochOnly = TRUE))
  # ## [1] "lengthDATA"     "grIntMu"        "grIntSigma"     "sigmaIntMu"
  # ## [5] "sigmaIntSigma"  "grBetaMu"       "grBetaSigma"    "sigmaBetaMu"
  # ## [9] "sigmaBetaSigma" "sigmaIndRE"     "grInt"          "sigmaInt"
  # ##[13] "grBeta"         "sigmaBeta"      "grIndRE"
  #
  #
  # ##for(nn in Rmodel$getVarNames(nodes = Rmodel$getNodeNames(stochOnly = TRUE))) {
  # ##    print(nn)
  # ##    print(any(is.na(Rmodel[[paste0('logProb_', nn)]])))
  # ##}
  # ##
  # ##Rmodel$logProb_lengthDATA
  # ##Rmodel$lengthDATA
  #
  #
  # conf <- configureMCMC(Rmodel)
  # ##(conf <- configureMCMC(Rmodel, useConjugacy = FALSE))
  #
  # ##conf$printSamplers()
  #
  # #conf$removeSamplers('sigmaIntSigma')
  #
  # #for(nn in Rmodel$expandNodeNames('sigmaIntSigma'))
  # #    conf$addSampler(nn, 'RW', control = list(log = TRUE))
  #
  # #conf$addSampler('sigmaIntSigma', 'RW_block')
  #
  # ##conf$printSamplers('sigmaIntSigma')
  #
  # conf$getMonitors()
  # conf$addMonitors(params)
  #
  # ##setdiff(params, conf$getMonitors())
  #
  # Rmcmc <- buildMCMC(conf)
  #
  # Cmodel <- compileNimble(Rmodel)
  #
  # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  #
  # ##args(runMCMC)
  #
  # mcmcInfo$nSamples <- (mcmcInfo$nIter - mcmcInfo$nBurnIn) * mcmcInfo$nChains
  #
  # mcmc <- runMCMC(Cmcmc, nburnin = mcmcInfo$nBurnIn, niter = mcmcInfo$nIter, nchains = mcmcInfo$nChains,
  #                 samples = TRUE, samplesAsCodaMCMC = TRUE,
  #                 summary = TRUE, WAIC = TRUE)

  return(list(code=code,data=data,constants=constants,inits=inits,params=params))
}


#'Run the growth model using jags
#'
#'@param d, the input dataframe
#'@param paralel, boolean for running chains in parallel
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  #grBetaOutside = array( runif( 11 *2*d$nSeasons*d$nRivers,-2,2),c( 11 ,2,d$nSeasons,d$nRivers))
  #grBetaOutside[ ,1,,2, ] <- 0

  nBetas <- 12
  nBetasSigma <- 6

  inits <- function(){
    list(
      grInt = array(rnorm(2*d$nSeasons*d$nRivers,0.5,0.25),c(2,d$nSeasons,d$nRivers)),
      #grBeta[x,1:2,1,1:4,1:4]
      grBeta = array(rnorm(nBetas*2*d$nSeasons*d$nRivers,0,0.1),c(nBetas,2,d$nSeasons,d$nRivers)),
      #grSigma[ yoy,spp,s,r ]
      grSigma = array(runif(2*d$nSeasons*d$nRivers,0,0.05),c(2,d$nSeasons,d$nRivers))
      # sigmaBeta[ b,yoy,spp,s,r ]
    #  sigmaBeta = array(rnorm(nBetasSigma*2*d$nSeasons*d$nRivers,0,0.05),c(nBetasSigma,2,d$nSeasons,d$nRivers))
   #   grIndRE = rnorm(d$nInd,0,0.001)
    )
  }

  # params <- c('grInt' ,'grIntMu','grIntSigma'
  #             ,'sigmaInt','sigmaIntMu','sigmaIntSigma'
  #             ,'grBeta','grBetaMu','grBetaSigma'
  #             , 'length','expectedGR'#,'expectedGR'
  #    #          , 'grIndRE','grIndREMean','grIndRETau'
  #             )

  params <- c('grInt','grIntMu','grIntSigma',
              'grBeta','grBetaMu', 'grBetaSigma',
              'sigmaInt','sigmaIntMu','sigmaIntSigma',
              'sigmaBeta','sigmaBetaMu', 'sigmaBetaSigma',
              'grIndRE', 'sigmaIndRE',
              'lengthExp')

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel7.jags",
                n.chains = 3,
                n.adapt = 500, #1000
                n.iter = 2500,
                n.burnin = 1500,
                n.thin = 5,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

#'Adjust counts
#'
#'@param cdIn core data
#'@param ddDIn results from the detection model
#'@param ddddGIn detection data
#'@param meanOrIterIn whether the growth model gets mean P from the detection model, or results from an iteration
#'@param sampleToUse if detection results are from an iteration, which iteration
#'@return a data frame with standardized counts, nAllFishBySpeciesPStd is standardized by species counts, nAllFishBySpeciesPAllSppStd is standardized by the sum of counts of all species in the analysis
#'@export
#'

adjustCounts <- function( cdIn,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

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
                          massAllFishStd =         ( massAllFish -          massAllFishMean )         /massAllFishSD
                        )

  cdIn <- left_join( cdIn,denForMerge2, by = c("species", "year", "season", "riverOrdered") )

  # counts are missing for samples with propSampled == 0. For now, fill in mean (0). Need this for when enc==0 and propSampled==0 in the gr model
  #should this be a low # instead of 0???
  cdIn$nAllFishBySpeciesPStd <-   ifelse( is.na(cdIn$nAllFishBySpeciesPStd),   0, cdIn$nAllFishBySpeciesPStd )
  cdIn$nAllFishPStd <-            ifelse( is.na(cdIn$nAllFishPStd),            0, cdIn$nAllFishPStd )
  cdIn$massAllFishBySpeciesStd <- ifelse( is.na(cdIn$massAllFishBySpeciesStd), 0, cdIn$massAllFishBySpeciesStd )
  cdIn$massAllFishStd <-          ifelse( is.na(cdIn$massAllFishStd),          0, cdIn$massAllFishStd )

  #####
  # counts by species in separate columns
  nAllFishBySpeciesPStdBySpp <- denForMerge2 %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesPStd) %>%
    spread(species, nAllFishBySpeciesPStd, fill=0) %>%
    rename(nAllFishBySpeciesPStdBKT = bkt,nAllFishBySpeciesPStdBNT = bnt, nAllFishBySpeciesPStdATS = ats)
  cdIn <- left_join( cdIn,nAllFishBySpeciesPStdBySpp )


  return(cdIn)

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


#'Remove fish with many intermediate NAs - these cause problems with indRE estimates
#'
#'@param d a dataframe
#'@return a data frame with problem fish removed
#'@export
#'

removeFishWithManyIntermediateNAs <- function(d){

  tagsToRemove <- c('00088cc02a','1c2d6c51d3','00088cc364')
  d <- d %>% filter( !(tag %in% tagsToRemove) )

  return(d)
}
