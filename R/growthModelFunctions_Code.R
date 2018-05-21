library(nimble)

codeSpp <- list()

## BKT
codeSpp[[1]] <- nimbleCode({
  ##
  for( i in 1:nEvalRows ) {
    ##
    mu[ evalRows[i] ] <- lengthDATA[ evalRows[i] ] + ( gr[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    lengthDATA[ evalRows[i] + 1 ] ~ dnorm( mu[ evalRows[i] ],
                                           sd = expectedGRSigma[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    ##

    gr[ evalRows[i] ] <-
      grInt[        isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +

      grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] +

      grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] +
      grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd2[evalRows[i]] +
      grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd2[evalRows[i]] +

      grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]]^2 +

      grBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
  #    grBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]]^2 +
  #    grBetaBNT[3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +

      grBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] +

      grIndRE[ ind[evalRows[i]] ]

    ##
    log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
      sigmaInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
      sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] +
      sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      sigmaBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
      #      sigmaBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +
      sigmaBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]]
    # sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * tempStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * tempStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * tempStd[evalRows[i]] * ATS01DATA[evalRows[i]] +
    # sigmaBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
    # sigmaBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * flowStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * flowStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * flowStd[evalRows[i]] * ATS01DATA[evalRows[i]]

  }
  ##
  for( s in 1:nSeasons ) {
    grIntMu[ s ] ~ dnorm(0,sd = 30) #changed from 0.001
    grIntSigma[ s ] ~ dunif(0,10)
    ##
    sigmaIntMu[ s ] ~ dnorm(0,sd = 30)
    sigmaIntSigma[ s ] ~ dunif(0,10)
    ##
    for( b in 1:nBetas ) {  ############# Beta ##############################################################################
      grBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      grBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( b in 1:nBetasSigma ) {  ############# Beta ##############################################################################
      sigmaBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      sigmaBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( yoy in 1:2 ) {
      for( r in 1:nRivers ) {
        grInt[ yoy,s,r ] ~ dnorm( grIntMu[ s ], sd = grIntSigma[ s ] )
        sigmaInt[ yoy,s,r ] ~ dnorm( sigmaIntMu[ s ], sd = sigmaIntSigma[ s ] )    #T(0,) #dunif(0,100)
        ##
        for( b in 1:nBetas ) {  ############# Beta ##################################################################################
          grBeta[ b,yoy,s,r ] ~ dnorm( grBetaMu[ b,s ], sd = grBetaSigma[ b,s ] )
        }
        for( b in 1:nBetasSigma ) {  ############# Beta ##################################################################################
          sigmaBeta[ b,yoy,s,r ] ~ dnorm( sigmaBetaMu[ b,s ], sd = sigmaBetaSigma[ b,s ] )
        }
      }
    }
  }

  # # nATS priors
  # for( b in 1:nBetasATS ) {
  #   for( yoy in 1:2 ) {
  #     for( s in 1:nSeasons ) {
  #       grBetaATS[ b,yoy,s ] ~ dnorm( grBetaATSMu[ b,s ], sd = grBetaATSSigma[ b,s ] )
  #       sigmaBetaATS[ b,yoy,s ] ~ dnorm( sigmaBetaATSMu[ b,s ], sd = sigmaBetaATSSigma[ b,s ] )
  #     }
  #   }
  # }
  # for( b in 1:nBetasATS ) {
  #   for( s in 1:nSeasons ) {
  #     grBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
  #     grBetaATSSigma[ b,s ] ~ dunif( 0,100 )
  #     sigmaBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
  #     sigmaBetaATSSigma[ b,s ] ~ dunif( 0,100 )
  #   }
  # }

  #nATS priors
  for( b in 1:nBetasATS ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {

        grBetaATS[ b,yoy,s,1 ] ~ dnorm( grBetaATSMu[ b,s ], sd = grBetaATSSigma[ b,s ] )
        grBetaATS[ b,yoy,s,2 ] <- 0
        grBetaATS[ b,yoy,s,3 ] <- 0
        grBetaATS[ b,yoy,s,4 ] <- 0

        sigmaBetaATS[ b,yoy,s,1 ] ~ dnorm( sigmaBetaATSMu[ b,s ], sd = sigmaBetaATSSigma[ b,s ] )
        sigmaBetaATS[ b,yoy,s,2 ] <- 0
        sigmaBetaATS[ b,yoy,s,3 ] <- 0
        sigmaBetaATS[ b,yoy,s,4 ] <- 0
      }
    }
  }
  for( b in 1:nBetasATS ) {
    for( s in 1:nSeasons ) {
      grBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaATSSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaATSSigma[ b,s ] ~ dunif( 0,100 )
    }
  }

  #nBNT priors
  for( b in 1:nBetasBNT ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {

        grBetaBNT[ b,yoy,s,1 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
        grBetaBNT[ b,yoy,s,2 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
                  grBetaBNT[ b,yoy,s,3 ] <- 0
                 grBetaBNT[ b,yoy,s,4 ] <- 0

        sigmaBetaBNT[ b,yoy,s,1 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
        sigmaBetaBNT[ b,yoy,s,2 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
                  sigmaBetaBNT[ b,yoy,s,3 ] <- 0
                  sigmaBetaBNT[ b,yoy,s,4 ] <- 0
      }
    }
  }
  for( b in 1:nBetasBNT ) {
    for( s in 1:nSeasons ) {
      grBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
    }
  }

  ##
  # for( yoy in 1:2 ) {
  #   for( s in 1:nSeasons ) {
  #     grATSOffset[ yoy,s ] ~ dnorm( 0, sd = 30 )
  #   }
  # }
  for(ii in 1:nEvalRows) {
    lengthExp[ evalRows[ii] + 1 ] <- lengthDATA[ evalRows[ii] ] + gr[ evalRows[ii] ] * sampleInterval[ evalRows[ii] ]
  }

  ############ standardized length for sizeBetas
  for( i in 1:(nEvalRows) ){
    lengthStd[ evalRows[i] ] <-  ( lengthDATA[ evalRows[i] ] - lengthMean ) / lengthSD
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


## BNT
codeSpp[[2]] <- nimbleCode({
  ##
  for( i in 1:nEvalRows ) {
    ##
    mu[ evalRows[i] ] <- lengthDATA[ evalRows[i] ] + ( gr[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    lengthDATA[ evalRows[i] + 1 ] ~ dnorm( mu[ evalRows[i] ],
                                           sd = expectedGRSigma[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    ##

    gr[ evalRows[i] ] <-
      grInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +

      grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] +

      grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]]+# * BKT01DATA[evalRows[i]] +
      grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd2[evalRows[i]] +
      grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd2[evalRows[i]] +

      #    grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * cBNTStd[evalRows[i]] +

      grBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
      #  grBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +

      grBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]] ] * cATSStd[evalRows[i]] +

      grIndRE[ ind[evalRows[i]] ]

    ##
    log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
      sigmaInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
      sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] +
      sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      sigmaBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
      #      sigmaBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +
      sigmaBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]] ] * cATSStd[evalRows[i]]
    # sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * tempStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * tempStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * tempStd[evalRows[i]] * ATS01DATA[evalRows[i]] +
    # sigmaBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
    # sigmaBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * flowStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * flowStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * flowStd[evalRows[i]] * ATS01DATA[evalRows[i]]

  }
  ##
  for( s in 1:nSeasons ) {
    grIntMu[ s ] ~ dnorm(0,sd = 30) #changed from 0.001
    grIntSigma[ s ] ~ dunif(0,10)
    ##
    sigmaIntMu[ s ] ~ dnorm(0,sd = 30)
    sigmaIntSigma[ s ] ~ dunif(0,10)
    ##
    for( b in 1:nBetas ) {  ############# Beta ##############################################################################
      grBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      grBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( b in 1:nBetasSigma ) {  ############# Beta ##############################################################################
      sigmaBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      sigmaBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( yoy in 1:2 ) {
      for( r in 1:nRivers ) {
        grInt[ yoy,s,r ] ~ dnorm( grIntMu[ s ], sd = grIntSigma[ s ] )
        sigmaInt[ yoy,s,r ] ~ dnorm( sigmaIntMu[ s ], sd = sigmaIntSigma[ s ] )    #T(0,) #dunif(0,100)
        ##
        for( b in 1:nBetas ) {  ############# Beta ##################################################################################
          grBeta[ b,yoy,s,r ] ~ dnorm( grBetaMu[ b,s ], sd = grBetaSigma[ b,s ] )
        }
        for( b in 1:nBetasSigma ) {  ############# Beta ##################################################################################
          sigmaBeta[ b,yoy,s,r ] ~ dnorm( sigmaBetaMu[ b,s ], sd = sigmaBetaSigma[ b,s ] )
        }
      }
    }
  }

  # nATS priors
  for( b in 1:nBetasATS ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {
        grBetaATS[ b,yoy,s ] ~ dnorm( grBetaATSMu[ b,s ], sd = grBetaATSSigma[ b,s ] )
        sigmaBetaATS[ b,yoy,s ] ~ dnorm( sigmaBetaATSMu[ b,s ], sd = sigmaBetaATSSigma[ b,s ] )
      }
    }
  }
  for( b in 1:nBetasATS ) {
    for( s in 1:nSeasons ) {
      grBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaATSSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaATSSigma[ b,s ] ~ dunif( 0,100 )
    }
  }


  #nBNT priors
  for( b in 1:nBetasBNT ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {

        grBetaBNT[ b,yoy,s,1 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
        grBetaBNT[ b,yoy,s,2 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
        #          grBetaBNT[ b,yoy,s,3 ] <- 0
        #         grBetaBNT[ b,yoy,s,4 ] <- 0

        sigmaBetaBNT[ b,yoy,s,1 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
        sigmaBetaBNT[ b,yoy,s,2 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
        #          sigmaBetaBNT[ b,yoy,s,3 ] <- 0
        #          sigmaBetaBNT[ b,yoy,s,4 ] <- 0
      }
    }
  }
  for( b in 1:nBetasBNT ) {
    for( s in 1:nSeasons ) {
      grBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
    }
  }

  ##
  # for( yoy in 1:2 ) {
  #   for( s in 1:nSeasons ) {
  #     grATSOffset[ yoy,s ] ~ dnorm( 0, sd = 30 )
  #   }
  # }
  for(ii in 1:nEvalRows) {
    lengthExp[ evalRows[ii] + 1 ] <- lengthDATA[ evalRows[ii] ] + gr[ evalRows[ii] ] * sampleInterval[ evalRows[ii] ]
  }

  ############ standardized length for sizeBetas
  for( i in 1:(nEvalRows) ){
    lengthStd[ evalRows[i] ] <-  ( lengthDATA[ evalRows[i] ] - lengthMean ) / lengthSD
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

## ATS
codeSpp[[3]] <- nimbleCode({
  ##
  for( i in 1:nEvalRows ) {
    ##
    mu[ evalRows[i] ] <- lengthDATA[ evalRows[i] ] + ( gr[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    lengthDATA[ evalRows[i] + 1 ] ~ dnorm( mu[ evalRows[i] ],
                                           sd = expectedGRSigma[ evalRows[i] ] * sampleInterval[ evalRows[i] ] )
    ##

    gr[ evalRows[i] ] <-
      grInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +

      grBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * lengthDATA[evalRows[i]] +

      grBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]]+# * BKT01DATA[evalRows[i]] +
      grBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      grBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      grBeta[ 5, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      grBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd2[evalRows[i]] +
      grBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd2[evalRows[i]] +

      #    grBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * cBNTStd[evalRows[i]] +

      grBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
      #  grBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +

      grBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]] ] * cATSStd[evalRows[i]] +

      grIndRE[ ind[evalRows[i]] ]

    ##
    log( expectedGRSigma[ evalRows[i] ] ) <- #grSigma
      sigmaInt[     isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] +
      sigmaBeta[ 1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] +
      sigmaBeta[ 2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] +
      sigmaBeta[ 3, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * flowStd[evalRows[i]] +
      sigmaBeta[ 4, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +

      sigmaBetaBNT[1, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] +
      #      sigmaBetaBNT[2, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * cBKTStd[evalRows[i]] +
      sigmaBetaATS[1, isYOYDATA[evalRows[i]], season[evalRows[i]] ] * cATSStd[evalRows[i]]
    # sigmaBeta[ 6, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * tempStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[ 7, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * tempStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[ 8, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * tempStd[evalRows[i]] * ATS01DATA[evalRows[i]] +
    # sigmaBeta[ 9, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * tempStd[evalRows[i]] * flowStd[evalRows[i]] +
    # sigmaBeta[10, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBKTStd[evalRows[i]] * flowStd[evalRows[i]] * BKT01DATA[evalRows[i]] +
    # sigmaBeta[11, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cBNTStd[evalRows[i]] * flowStd[evalRows[i]] * BNT01DATA[evalRows[i]] +
    # sigmaBeta[12, isYOYDATA[evalRows[i]], season[evalRows[i]], riverDATA[evalRows[i]] ] * cATSStd[evalRows[i]] * flowStd[evalRows[i]] * ATS01DATA[evalRows[i]]

  }
  ##
  for( s in 1:nSeasons ) {
    grIntMu[ s ] ~ dnorm(0,sd = 30) #changed from 0.001
    grIntSigma[ s ] ~ dunif(0,10)
    ##
    sigmaIntMu[ s ] ~ dnorm(0,sd = 30)
    sigmaIntSigma[ s ] ~ dunif(0,10)
    ##
    for( b in 1:nBetas ) {  ############# Beta ##############################################################################
      grBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      grBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( b in 1:nBetasSigma ) {  ############# Beta ##############################################################################
      sigmaBetaMu[ b,s ] ~ dnorm(0,sd = 30)
      sigmaBetaSigma[ b,s ] ~ dunif(0,100)
    }
    ##
    for( yoy in 1:2 ) {
      for( r in 1:nRivers ) {
        grInt[ yoy,s,r ] ~ dnorm( grIntMu[ s ], sd = grIntSigma[ s ] )
        sigmaInt[ yoy,s,r ] ~ dnorm( sigmaIntMu[ s ], sd = sigmaIntSigma[ s ] )    #T(0,) #dunif(0,100)
        ##
        for( b in 1:nBetas ) {  ############# Beta ##################################################################################
          grBeta[ b,yoy,s,r ] ~ dnorm( grBetaMu[ b,s ], sd = grBetaSigma[ b,s ] )
        }
        for( b in 1:nBetasSigma ) {  ############# Beta ##################################################################################
          sigmaBeta[ b,yoy,s,r ] ~ dnorm( sigmaBetaMu[ b,s ], sd = sigmaBetaSigma[ b,s ] )
        }
      }
    }
  }

  # nATS priors
  for( b in 1:nBetasATS ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {
        grBetaATS[ b,yoy,s ] ~ dnorm( grBetaATSMu[ b,s ], sd = grBetaATSSigma[ b,s ] )
        sigmaBetaATS[ b,yoy,s ] ~ dnorm( sigmaBetaATSMu[ b,s ], sd = sigmaBetaATSSigma[ b,s ] )
      }
    }
  }
  for( b in 1:nBetasATS ) {
    for( s in 1:nSeasons ) {
      grBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaATSSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaATSMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaATSSigma[ b,s ] ~ dunif( 0,100 )
    }
  }


  #nBNT priors
  for( b in 1:nBetasBNT ) {
    for( yoy in 1:2 ) {
      for( s in 1:nSeasons ) {

        grBetaBNT[ b,yoy,s,1 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
        #grBetaBNT[ b,yoy,s,2 ] ~ dnorm( grBetaBNTMu[ b,s ], sd = grBetaBNTSigma[ b,s ] )
        #          grBetaBNT[ b,yoy,s,3 ] <- 0
        #         grBetaBNT[ b,yoy,s,4 ] <- 0

        sigmaBetaBNT[ b,yoy,s,1 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
        #sigmaBetaBNT[ b,yoy,s,2 ] ~ dnorm( sigmaBetaBNTMu[ b,s ], sd = sigmaBetaBNTSigma[ b,s ] )
        #          sigmaBetaBNT[ b,yoy,s,3 ] <- 0
        #          sigmaBetaBNT[ b,yoy,s,4 ] <- 0
      }
    }
  }
  for( b in 1:nBetasBNT ) {
    for( s in 1:nSeasons ) {
      grBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      grBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
      sigmaBetaBNTMu[ b,s ] ~ dnorm( 0,sd = 30 )
      sigmaBetaBNTSigma[ b,s ] ~ dunif( 0,100 )
    }
  }

  ##
  # for( yoy in 1:2 ) {
  #   for( s in 1:nSeasons ) {
  #     grATSOffset[ yoy,s ] ~ dnorm( 0, sd = 30 )
  #   }
  # }
  for(ii in 1:nEvalRows) {
    lengthExp[ evalRows[ii] + 1 ] <- lengthDATA[ evalRows[ii] ] + gr[ evalRows[ii] ] * sampleInterval[ evalRows[ii] ]
  }

  ############ standardized length for sizeBetas
  for( i in 1:(nEvalRows) ){
    lengthStd[ evalRows[i] ] <-  ( lengthDATA[ evalRows[i] ] - lengthMean ) / lengthSD
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

