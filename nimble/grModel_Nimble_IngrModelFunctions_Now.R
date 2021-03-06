

library(nimble)
library(arrayhelpers)
library(zoo)
library(jagsUI)
library(coda)
library(tidyverse)

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


load('./data/out/dG_1515095839_bkt1997.RData')
jd <- ddG[[1]][[1]]

##tail(ddG[[1]][[2]])
##str(jd)
##names(jd)
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

constants <- list(riverDATA = jd$riverDATA,
                  nRivers = jd$nRivers,
                  nSpecies = jd$nSpecies,
                  ind = jd$ind,
                  nInd = jd$nInd,
                  season = jd$season,
                  nEvalRows = jd$nEvalRows,
                  evalRows = jd$evalRows,
                  nSeasons = jd$nSeasons,
                  countPAllSppStd = jd$countPAllSppStd,
                  isYOYDATA = jd$isYOYDATA,
                  tempStd = jd$tempStd,
                  flowStd = jd$flowStd)
##
#data <- list(lengthDATA = jd$lengthDATA)
data <- list(lengthDATA = jd$lengthDATA[1:33518])
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



Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$lengthDATA <- zoo::na.approx(data$lengthDATA[1:33518]) ## length(data$lengthDATA) = 33519, last obs is a fish with a single obs - doesn't get an evalRow
table(is.na(Rmodel$lengthDATA))

#system.time(lp <- Rmodel$calculate())
#lp

##Rmodel$getVarNames(nodes = Rmodel$getNodeNames(stochOnly = TRUE))
## [1] "lengthDATA"     "grIntMu"        "grIntSigma"     "sigmaIntMu"
## [5] "sigmaIntSigma"  "grBetaMu"       "grBetaSigma"    "sigmaBetaMu"
## [9] "sigmaBetaSigma" "sigmaIndRE"     "grInt"          "sigmaInt"
##[13] "grBeta"         "sigmaBeta"      "grIndRE"


##for(nn in Rmodel$getVarNames(nodes = Rmodel$getNodeNames(stochOnly = TRUE))) {
##    print(nn)
##    print(any(is.na(Rmodel[[paste0('logProb_', nn)]])))
##}
##
##Rmodel$logProb_lengthDATA
##Rmodel$lengthDATA


conf <- configureMCMC(Rmodel)
##(conf <- configureMCMC(Rmodel, useConjugacy = FALSE))

##conf$printSamplers()

#conf$removeSamplers('sigmaIntSigma')

#for(nn in Rmodel$expandNodeNames('sigmaIntSigma'))
#    conf$addSampler(nn, 'RW', control = list(log = TRUE))

#conf$addSampler('sigmaIntSigma', 'RW_block')

##conf$printSamplers('sigmaIntSigma')

conf$getMonitors()
conf$addMonitors(params)

##setdiff(params, conf$getMonitors())

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##args(runMCMC)
mcmcInfo <- list()
mcmcInfo$nChains <- 3
mcmcInfo$nIter <- 3000
mcmcInfo$nBurnIn <- 1500
mcmcInfo$nSamples <- (mcmcInfo$nIter - mcmcInfo$nBurnIn) * mcmcInfo$nChains

mcmc <- runMCMC(Cmcmc, nburnin = mcmcInfo$nBurnIn, niter = mcmcInfo$nIter, nchains = mcmcInfo$nChains,
                               samples = TRUE, samplesAsCodaMCMC = TRUE,
                               summary = TRUE, WAIC =TRUE)

# Use jagsUI functions to get output ionto jagsUI format for analysis functions
source("./R/jagsUIFunctions.R")
mcmcProcessed <- process.output(mcmc$samples, DIC=FALSE, params.omit=FALSE)#,params.omit=("lengthExp")) #jagsUI function, DIC = FALSE because it requires 'deviance'

obsPred <- getRMSE_Nimble(mcmcProcessed,0.6)
obsPred$rmse
obsPred$outliers%>% as.data.frame()

#####
source("./R/plotFunctions.R") #will put into functions later
plotInt_Nimble(mcmcProcessed)
plotBetas_Nimble(mcmcProcessed,1:2)
plotBetas_Nimble(mcmcProcessed,3:4)
plotBetas_Nimble(mcmcProcessed,5:18)

plotSigmaInt_Nimble(mcmcProcessed)
plotSigmaBetas_Nimble(mcmcProcessed,1:2)
plotSigmaBetas_Nimble(mcmcProcessed,3:4)

#########
limits <- 1.5 # -/+ limits on standardized range of input variable
nPoints <- 5

nItersForPred <- 100
itersForPred <- sample( 1:mcmcInfo$nSamples, nItersForPred )

# predictions across the grid
p <- getPrediction( mcmcProcessed, limits, nPoints, itersForPred, constants, c("len", "temp", "flow","count") )
#######################
# graph function, in analyzeOutputFunctions.R

plotPred(p, "len", 1, "bkt") #spp is just a pass-through for title until I combine the species results into one df
plotPred(p, "temp", 1, "bkt")
plotPred(p, "flow", 1, "bkt")
plotPred(p, "count", 1, "bkt")
plotPred(p, c("flow", "temp"), 1, "bkt")
plotPred(p, c("temp", "flow"), 0, "bkt")
plotPred(p, c("temp","count"), 1, "bkt")
plotPred(p, c("flow","count"), 1, "bkt")
plotPred(p, c("len","count"), 1, "bkt")
plotPred(p, c("flow","len"), 1, "bkt")




























extractSamples <- function(samples, varIn, thin = 3) {
  outDFAll <- data.frame()

  # more than 1 chain
  if(is.list(samples)) { nChains <- length(samples)
  ## process samples
  var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', varIn))   ## add \\ before any '[' or ']' appearing in var
  var <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples[[1]]) , value=TRUE)))  ## expanded any indexing

  for( i in 1:nChains ){
    out <- samples[[i]][,var]
    outDF <- array2df(out) %>% rename(estimate = out, rowNumber = d1, variable = d2) %>% mutate( chain = i )
    outDFAll <- bind_rows(outDFAll,outDF)
  }

  # 1 chain
  } else {
    nChains <- 1
    ## process samples
    var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', varIn))   ## add \\ before any '[' or ']' appearing in var
    var <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
    out <- samples[,var]
    outDFAll <- array2df(out) %>% rename(estimate = out, rowNumber = d1, variable = d2) %>% mutate( chain = 1 )
  }
  outDFAll <- outDFAll %>% filter( rowNumber %% thin == 0)
  return(outDFAll)
}


mvSamples <- samples2$mvSamples
samplesMatrix <- as.matrix(mvSamples)

str(samples2)
dim(samples)
colnames(samples)
samplesPlot(samples$samples, 'sigmaIntSigma')






###############################
# or, do all in one step
# system.time(
#   mcmcOut <- nimbleMCMC(code = code, constants = constants,
#                         data = data, inits = inits,
#                         nchains = 2, niter = 1000,
#                         summary = TRUE, WAIC = TRUE)
# )


library(coda)
ess <- apply(mcmcProcessed, 2, effectiveSize)

ess5 <- ess * 5



summary(ess5)

samplesPlot(samples, 'sigmaIntSigma')













