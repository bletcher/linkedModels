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
      mu[ evalRows[i] ] <- lengthDATA[ evalRows[i] ] + ( gr[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
      lengthDATA[ evalRows[i] + 1 ] ~ dnorm( mu[ evalRows[i] ],
                                             sd = expectedGRSigma[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
      ##

      gr[ evalRows[i] ] <-
        grInt[         isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +

          grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] +
          grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] +
          grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] +
          grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] +
          grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
          grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      # #
          grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]]^2 +
          grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]]^2 +
          grBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]]^2 +
          grBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]]^2 +
          grBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]]^2 +
          grBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]]^2 +
      #
          grBeta[13, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[14, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[15, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[16, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[17, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * countPStdBKT[evalRows[i]] +
          grBeta[18, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * countPStdBNT[evalRows[i]] +
          grBeta[19, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * countPStdATS[evalRows[i]] +
      #
          grBeta[20, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[21, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
          grBeta[22, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

          grIndRE[ ind[evalRows[i]] ]


      # gr[ evalRows[i] ] <-
      #   grInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
      #   grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthStd[evalRows[i]] +
      #   grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] +
      #   grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      #   grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      #   ## grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]]^2 +
      #   grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]]^2 +
      #   grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]]^2 +
      #   grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]]^2 +
      #   grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
      #   grBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * countPAllSppStd[evalRows[i]] +
      #   grBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * flowStd[evalRows[i]] +
      #   grBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
      #   ## grBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * countPAllSppStd[evalRows[i]] +
      #   ## grBeta[13, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * flowStd[evalRows[i]] +
      #   ## grBeta[14, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * tempStd[evalRows[i]] +
      #   ## grBeta[16, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
      #   ## grBeta[17, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * lengthDATA[evalRows[i]] * flowStd[evalRows[i]] +
      #   ## grBeta[18, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] * lengthDATA[evalRows[i]] +
      #   grIndRE[ ind[evalRows[i]] ]
      ##
      log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
        sigmaInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
        sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] +
        sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] +
        sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] +
        sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
        sigmaBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
        sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] * tempStd[evalRows[i]] +
        sigmaBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] * tempStd[evalRows[i]] +
        sigmaBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] * tempStd[evalRows[i]] +
        sigmaBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
        sigmaBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBKT[evalRows[i]] * flowStd[evalRows[i]] +
        sigmaBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdBNT[evalRows[i]] * flowStd[evalRows[i]] +
        sigmaBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPStdATS[evalRows[i]] * flowStd[evalRows[i]]


      # log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
      #   sigmaInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
      #   sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] +
      #   sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      #   sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      #   sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * tempStd[evalRows[i]] +
      #   sigmaBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
      #   sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * countPAllSppStd[evalRows[i]] * flowStd[evalRows[i]]



    }
    ##
    for( s in 1:nSeasons ) {
      grIntMu[ s ] ~ dnorm(0,sd = 30) #changed from 0.001
      grIntSigma[ s ] ~ dunif(0,10)
      ##
      sigmaIntMu[ s ] ~ dnorm(0,sd = 30)
      sigmaIntSigma[ s ] ~ dunif(0,10)
      ##
      for( b in 1:22 ) {  ############# Beta ##############################################################################
        grBetaMu[ b,s ] ~ dnorm(0,sd = 30)
        grBetaSigma[ b,s ] ~ dunif(0,100)
      }
      ##
      for( b in 1:12 ) {  ############# Beta ##############################################################################
        sigmaBetaMu[ b,s ] ~ dnorm(0,sd = 30)
        sigmaBetaSigma[ b,s ] ~ dunif(0,100)
      }
      ##
      for( yoy in 1:2 ) {
        for( r in 1:nRivers ) {
          grInt[ yoy,s,r ] ~ dnorm( grIntMu[ s ], sd = grIntSigma[ s ] )
          sigmaInt[ yoy,s,r ] ~ dnorm( sigmaIntMu[ s ], sd = sigmaIntSigma[ s ] )    #T(0,) #dunif(0,100)
          ##
          for( b in 1:22 ) {  ############# Beta ##################################################################################
            grBeta[ b,yoy,s,r ] ~ dnorm( grBetaMu[ b,s ], sd = grBetaSigma[ b,s ] )
          }
          for( b in 1:12 ) {  ############# Beta ##################################################################################
            sigmaBeta[ b,yoy,s,r ] ~ dnorm( sigmaBetaMu[ b,s ], sd = sigmaBetaSigma[ b,s ] )
          }
        }
      }
    }
    ##
    for(ii in 1:nEvalRows) {
      lengthExp[ evalRows[ii] + 1 ] <- lengthDATA[ evalRows[ii] ] + gr[ evalRows[ii] ] * sampleInterval[ evalRows[ii] ]
    }

    ############ standardized length for sizeBetas
    for( i in 1:(nEvalRows) ){
      lengthStd[ evalRows[i] ] <-  ( lengthDATA[ evalRows[i] ] - lengthMean ) / lengthSD
    }


   # for(ii in 1:nEvalRows) {
   #  grExp[ evalRows[ii] ] <- gr[ evalRows[ii] ] / sampleInterval[ evalRows[ii] ]
  #  }

    # for(y in 1:2) {
    #   for(s in 1:nSeasons){
    #     for(r in 1:nRivers){
    #       for( len in 1:7 ){
    #         for( count in 1:7 ){
    #           for( temp in 1:7 ){
    #             for( flow in 1:7 ){
    #
    #               predGR[ y, s, r, len,count,temp,flow ] <-
    #                   grInt[ y,s,r ] +
    #                   grBeta[ 1, y,s,r ] * (-2 + len*0.5) +
    #                   grBeta[ 2, y,s,r ] * (-2 + count*0.5) +
    #                   grBeta[ 3, y,s,r ] * (-2 + temp*0.5) +
    #                   grBeta[ 4, y,s,r ] * (-2 + flow*0.5) +
    #                   grBeta[ 5, y,s,r ] * (-2 + count*0.5)^2 +
    #                   grBeta[ 6, y,s,r ] * (-2 + temp*0.5)^2 +
    #                   grBeta[ 7, y,s,r ] * (-2 + flow*0.5)^2 +
    #                   grBeta[ 8, y,s,r ] * (-2 + temp*0.5) * (-2 + flow*0.5) +
    #                   grBeta[ 9, y,s,r ] * (-2 + temp*0.5) * (-2 + count*0.5) +
    #                   grBeta[10, y,s,r ] * (-2 + count*0.5) * (-2 + flow*0.5) +
    #                   grBeta[11, y,s,r ] * (-2 + count*0.5) * (-2 + temp*0.5) * (-2 + flow*0.5)
    #
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }

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

  # Need to trim off fish at the end of the data that aer only observed once and aren't in evalRows


  constants <- list(riverDATA = d$riverDATA,
                    nRivers = d$nRivers,
                    nSpecies = d$nSpecies,
                    ind = d$ind,
                    nInd = d$nInd,
                    season = d$season,
                    nEvalRows = d$nEvalRows,
                    evalRows = d$evalRows,
                    nSeasons = d$nSeasons,
                  #  countPAllSppStd = d$countPAllSppStd,
                    countPStdBKT = d$countPStdBKT,
                    countPStdBNT = d$countPStdBNT,
                    countPStdATS = d$countPStdATS,
                    isYOYDATA = d$isYOYDATA,
                    tempStd = d$tempStd,
                    flowStd = d$flowStd,
                    #lengthDATAStd = d$lengthDATAStd,
                    lengthMean = d$lengthMean,
                    lengthSD = d$lengthSD,
                    sampleInterval = d$sampleInterval[1:(max(d$evalRows)+1)]
              #      predRange = seq(-15,15,5)
                    )
  ##
  data <- list(lengthDATA = d$lengthDATA[1:(max(constants$evalRows)+1)]) # so don't have trailing single obs fish at end of df
  print(paste("Trimmed", length(d$lengthDATA) - length(data$lengthDATA), "fish that had single observation(s) at end of df"))
  ##
  nBetas <- 22
  nBetasSigma <- 12
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
#'@return Std n's for missing occasion samples (e.g. some winters) are interpolated
#'@export
#'

adjustCounts <- function( cdIn,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

  if ( meanOrIterIn == "mean") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    #ggplot(filter(den,!is.na(nAllFishBySpeciesP)),aes(yearN,nAllFishBySpeciesP, color = factor(speciesN))) + geom_point() + geom_line() + facet_grid(riverN~seasonN)
    #ftable(den$yearN,den$speciesN,is.na(den$nAllFishBySpeciesP))
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
  #  filter( !is.na(nAllFishBySpeciesP) ) %>%
    arrange(species,riverOrdered,year,season)

  # stats for counts by species
  denForMergeSummaryBySpecies <- denForMerge %>%
    group_by( species,season,riverOrdered ) %>%
    summarize( nAllFishBySpeciesPMean = mean( nAllFishBySpeciesP, na.rm = T ),
               nAllFishBySpeciesPSD =     sd( nAllFishBySpeciesP, na.rm = T),
               massAllFishBySpeciesMean = mean( massAllFishBySpecies, na.rm = T ),
               massAllFishBySpeciesSD =     sd( massAllFishBySpecies, na.rm = T))

  # stats for counts for all species
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

  # create template for all possible occasions
  possibleOccasions <- data.frame(table(denForMerge$season,denForMerge$year,denForMerge$species,denForMerge$riverOrdered)) %>%
    mutate( FreqUp = lead(Freq), FreqDown = lag(Freq),
            occFilled = ifelse(FreqUp * FreqDown + Freq > 0,1,0) ) %>%
    filter(occFilled == 1) %>%
    rename( season=Var1, year=Var2, species=Var3, riverOrdered=Var4 ) %>%
    mutate( season = as.numeric(season), year = as.numeric(year) + min(denForMerge2$year, na.rm = TRUE) - 1, species = as.character(species),
            riverOrdered = factor(riverOrdered,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)
          )
  # use this as merge template, then interpolate std values

  denForMerge_possibleOccasions <- left_join( possibleOccasions,denForMerge2 ) %>%
    mutate( nAllFishBySpeciesPStd = zoo::na.approx(nAllFishBySpeciesPStd),
            nAllFishPStd = zoo::na.approx(nAllFishPStd),
            massAllFishBySpeciesStd = zoo::na.approx(massAllFishBySpeciesStd),
            massAllFishStd = zoo::na.approx(massAllFishStd)) %>%
    dplyr::select(-c(Freq,FreqUp,FreqDown))

  cdIn <- left_join( cdIn,denForMerge_possibleOccasions, by = c("species", "year", "season", "riverOrdered") )


  #####
  # counts by species in separate columns
  nAllFishBySpeciesPStdBySpp <- denForMerge_possibleOccasions %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesPStd) %>%
    spread(key=species, value=nAllFishBySpeciesPStd, fill = -2.5) %>%
    rename(nAllFishBySpeciesPStdBKT = bkt,
           nAllFishBySpeciesPStdBNT = bnt,
           nAllFishBySpeciesPStdATS = ats
           )

  cdIn <- left_join( cdIn,nAllFishBySpeciesPStdBySpp )

 #  #####
 #  # counts by species and yoy in separate columns
 #  # yoy1 and yoy2 are too highly correlated - don't use
 #  nAllFishBySpeciesPStdBySppYOY <- denForMerge_possibleOccasions %>%
 #    dplyr::select(isYOYN,species, season, riverOrdered, year, nAllFishBySpeciesPStd) %>%
 #    unite(yoySpp, species, isYOYN, sep = "_") %>%
 #    spread(key=yoySpp, value=nAllFishBySpeciesPStd, fill = -2.5) %>%
 #    rename(nAllFishBySpeciesPStdBKT_yoy1 = bkt_1,
 #           nAllFishBySpeciesPStdBKT_yoy2 = bkt_2,
 #           nAllFishBySpeciesPStdBNT_yoy1 = bnt_1,
 #           nAllFishBySpeciesPStdBNT_yoy2 = bnt_2,
 #           nAllFishBySpeciesPStdATS_yoy1 = ats_1,
 #           nAllFishBySpeciesPStdATS_yoy2 = ats_2)
 #
 #  cdIn <- left_join( cdIn,nAllFishBySpeciesPStdBySppYOY )
 # # ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdBKT_yoy1,nAllFishBySpeciesPStdBKT_yoy2)) + geom_point()
 # # ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdBNT_yoy1,nAllFishBySpeciesPStdBNT_yoy2)) + geom_point()
 # #  ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdATS_yoy1,nAllFishBySpeciesPStdATS_yoy2)) + geom_point()
 #  nAllFishBySpeciesPStdBySppYOY[nAllFishBySpeciesPStdBySppYOY == -2.5] <- NA
 #  round(cor(nAllFishBySpeciesPStdBySppYOY[,c('nAllFishBySpeciesPStdBKT_yoy1','nAllFishBySpeciesPStdBKT_yoy2','nAllFishBySpeciesPStdBNT_yoy1','nAllFishBySpeciesPStdBNT_yoy2','nAllFishBySpeciesPStdATS_yoy1','nAllFishBySpeciesPStdATS_yoy2')], use="complete.obs"),2)

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

  tagsToRemove <- c('00088cc02a','1c2d6c51d3','00088cc364','1c2c582218')
  d <- d %>% filter( !(tag %in% tagsToRemove) )

  return(d)
}
