
#install.packages("devtools")
# library(devtools)
#devtools::install_github('bletcher/linkedModels')
#devtools::install_github(repo = "Conte-Ecology/westBrookData", subdir = "getWBData")
library(arm)
library(zoo)
library(arrayhelpers)
library(linkedModels)
library(coda)
library(jagsUI)
library(nimble)
library(getWBData)
library(stringr)
library(lubridate)
library(dbplyr)
library(tidyverse)
######################################################
# selection criteria

reconnect()

# run to here to start so can reconnect

drainage <- "west" # ==

speciesDet <- c("bkt", "bnt","ats") #keep as all three spp
speciesInDet <- factor(speciesDet, levels = c('bkt','bnt','ats'), ordered = T)

speciesGr <- "ats"
#speciesGr = c("bkt", "bnt","ats")
speciesInGr <- factor(speciesGr, levels = c('bkt','bnt','ats'), ordered = T)

riverOrderedIn <- factor(c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)

minCohort <- 1997#1995 # >=
maxSampleInterval <- 200 # <

recreatecdFileBeforeDetMod_TF <- FALSE
runDetectionModel_TF <- FALSE
recreateCD_TF <- FALSE
runCrossValidation_TF <- TRUE ############################
percentLeftOut <- 10

meanOrIter <- "mean"; iter <- 1
####### or ########
#meanOrIter = "iter"

######################################################
# Get data from database
# Only need to run this when data or functions have changed
# or you want to change drainages

cdFileBeforeDetMod <- paste0('./data/cd_',drainage,"_",minCohort,'_BeforeDetMod.RData')
if ( recreatecdFileBeforeDetMod_TF ) {
  # core data - tagged fish
  cdBeforeDetMod <- getCoreData(drainage) %>%
    cleanData(drainage) %>%
    mergeSites(drainage) %>%
    mutate(drainage = drainage,
           countP = NA)  # placeholder so prepareDataForJagsNimble() works for detection model

  cdBeforeDetMod$isYOY <- ifelse( cdBeforeDetMod$ageInSamples <= 3, 1, 2 )

  cdBeforeDetMod <- cdBeforeDetMod %>%

   # addRawCounts(drainage, filteredAreas = c("inside","trib")) %>% # counts and summed masses of all fish, tagged and untagged. Not adjusted for P (happens in addDensities())
   # do this in detection model, removeUnsampledRows(drainage, removeIncomplete = T) %>% # removes enc=0 rows for samples with no or very few() sections sampled (p=0 otherwise for those samples)
    removeLowAbundanceRivers(drainage) # removes ats,jimmy fish and bnt,mitchell from drainage=='west' - too few fish for estimates
  #
  #ftable(cdBeforeDetMod$isYOY,cdBeforeDetMod$species,cdBeforeDetMod$river,cdBeforeDetMod$season,cdBeforeDetMod$year)

  save(cdBeforeDetMod, file = cdFileBeforeDetMod)
} else {
  load(cdFileBeforeDetMod)
}
#################################
# Detection model
# all species at once
#################################
#
#
dModelName <- paste0(paste0(speciesDet,collapse = ''),"_",minCohort)

(start <- Sys.time())

ddddD <- cdBeforeDetMod %>%
  filter(  species %in% speciesInDet,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval) %>%  # this removes the later yearly samples. Want to stick with seasonal samples
           removeUnsampledRows(drainage, removeIncomplete = TRUE) # removes enc=0 rows for samples with no or very few() sections sampled (p=0 otherwise for those samples)

dddD <- ddddD %>% prepareDataForJags_Nimble('detection')

save(ddddD,dddD, file = paste0('./data/out/ddddD_', dModelName,'.RData'))

#################################################################################################
if (runDetectionModel_TF) {
  ddD <- dddD[[1]] %>% runDetectionModel(parallel = TRUE) ###################### currently for jags
  save(ddD, file = paste0('./data/out/ddD_', dModelName,'.RData'))
  #whiskerplot(ddD, parameters="pBetaInt")
} else {
  load( file = paste0('./data/out/ddD_', dModelName,'.RData') )
}

cdFile <- paste0('./data/cd_',drainage,"_",minCohort,'.RData')
if (recreateCD_TF) {

  ##########################################################################
  # get counts of all fish, tagged and untagged

  allFishBySpeciesYOY <- getRawCounts(drainage, filteredAreas = c("inside","trib"))
  allFishBySpeciesYOY$riverOrdered <- factor(allFishBySpeciesYOY$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)

  allFishBySpeciesYOY <- adjustCounts_allFish( allFishBySpeciesYOY,ddD,ddddD,meanOrIter,iter )

  # nAllFishBySpeciesPATS (and other spp) have NA when no fish were caught in a river. Set to 0 for ATS in WB so we can test effect of 0 fish there
  allFishBySpeciesYOY$nAllFishBySpeciesPATS <- ifelse( is.na(allFishBySpeciesYOY$nAllFishBySpeciesPATS) & allFishBySpeciesYOY$river == "west brook",
                                                          0, allFishBySpeciesYOY$nAllFishBySpeciesPATS)

#  ftable(allFishBySpeciesYOY$isYOY,allFishBySpeciesYOY$species,allFishBySpeciesYOY$river,allFishBySpeciesYOY$season,allFishBySpeciesYOY$year)
#  ggplot(allFishBySpeciesYOY, aes(year,nAllFishBySpeciesPATS)) + geom_point() + facet_grid(riverOrdered~season)
#  ggplot(allFishBySpeciesYOY, aes(year,nAllFishBySpeciesPBNT)) + geom_point() + facet_grid(riverOrdered~season)
  ###########################################################################

  cd <- left_join( cdBeforeDetMod,allFishBySpeciesYOY )

  # add counts derived from counts/p to cdBeforeDetMod. ddD and ddddD need to have all three spp to create spp-specific abundances
  #cd <- adjustCounts( cdBeforeDetMod,ddD,ddddD,meanOrIter,iter ) # function in growthModelFunctions.R
  save(cd, file = cdFile)
} else {
  load( file = cdFile )
}


# with full intercepts, effect of nPasses is insignificant (massive overlap). Leaving, nPasses out. Just going with pBetaInt.
#whiskerplot(ddD, parameters = "pBetaInt")
#whiskerplot(ddD, parameters = "pBeta")

#################################
# Movement model
#################
# To start, will just assume fish don't move
# until they are seen somewhere else.
# Accomplished with %>% fillSizeLocation(size = F) in getCoreData()

##############################################################################################################################
# Growth model
# One species at a time
###############

#########################################
# Filter raw data (cd)
ddddG <- cd %>%
  filter(  species %in% speciesGr,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval, # this removes the later yearly samples. Want to stick with seasonal samples
           #enc == 1 #change to btw first and last so we can impute missing obs?
           knownZ == 1

           #riverN == 1 ############
  )
nSeasons <- n_distinct(ddddG$season, na.rm=T)

# merge in density data, either overall means or single iterations at a time
# meanOrIter ="mean" uses means of all iterations
# meanOrIter ="iter" uses the sample that is the combo of sampleToUse and chainToUse (ignored if meanOrIter="mean")
#  for numOfItersToUse samples

#ggplot(ddddG, aes(sampleNumber,nAllFishBySpeciesPStdATS)) + geom_point() + facet_grid(riverOrdered~season)
#ggplot(ddddG, aes(sampleNumber,nAllFishBySpeciesPStdBNT)) + geom_point() + facet_grid(riverOrdered~season)
#ggplot(ddddG, aes(sampleNumber,nAllFishBySpeciesPStdBKT)) + geom_point() + facet_grid(riverOrdered~season)

######################################
######################################

chainToUse <- 1
numItersToUse <- 2
if (meanOrIter == "iter") {
  itersToUse <- sort(sample(((dd$mcmc.info$n.samples/dd$mcmc.info$n.chains) * (chainToUse - 1)):
                              ((dd$mcmc.info$n.samples/dd$mcmc.info$n.chains) * (chainToUse - 0)),
                            numItersToUse))
} else {
  itersToUse <- 1 # only run one loop over iter
}

######################################
######################################
### loop over iters

modelName <- paste0(paste0(speciesGr,collapse = ''),minCohort)

dddG <- list()
ddG <- list()
dG <- list()
print(elapsed)
elapsed <- list()
# run the growth model for detection model iterations in itersToUse
ii <- 0
iter=1
#######for (iter in itersToUse) {
  ii <- ii + 1

  #########################################
  # Prepare data for model run

  # dddG[[ii]] <- addBiomassDeltas( dddG[[ii]] )
  dddG[[ii]] <- addSurvivals( ddddG,ddD,meanOrIter,iter ) %>% removeFishWithManyIntermediateNAs()
  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags_Nimble("growth") # returns a list, list for running model in [[ii]][[1]], input df in [[ii]][[2]]
    #save(dG,ddG,dddG, file = paste0('./data/out/dG_', modelName, '_forLmer.RData')) # for Lmer model with all spp


  # #########################################
  # # mcmc run data
  # mcmcInfo <- list()
  # mcmcInfo$nChains <- 3
  # mcmcInfo$nIter <- 1000 #25000
  # mcmcInfo$thinRate <- 10 # nIter gets thinned
  # mcmcInfo$nSamples <- 25 # Number of thinned samples
  #
  # try( if(mcmcInfo$nSamples > mcmcInfo$nIter/mcmcInfo$thinRate) stop( 'nSamples > thinned nIter' ) )
  #
  # mcmcInfo$nBurnIn <- mcmcInfo$nIter/mcmcInfo$thinRate - mcmcInfo$nSamples
  # #mcmcInfo$nSamples <- round( (mcmcInfo$nIter - mcmcInfo$nBurnIn) * mcmcInfo$nChains / mcmcInfo$thinRate )

  #########################################
  # mcmc run data
  mcmcInfo <- list()
  mcmcInfo$nChains <- 3
  mcmcInfo$nIter <- 50000 #25000
  mcmcInfo$nBurnIn <- 25000
  mcmcInfo$AvailForSampling <- mcmcInfo$nIter - mcmcInfo$nBurnIn
  mcmcInfo$thinRate <- 10
  mcmcInfo$nSamples <- mcmcInfo$AvailForSampling / mcmcInfo$thinRate

  #########################################
  # Get data for model run
  source('./R/growthModelFunctions_Code.R') #################################################
  code <- codeSpp[[as.numeric(speciesInGr)]]

  #####
  print(c("in loop",meanOrIter,ii,iter))
  nB <- list()

  if(speciesGr == 'bkt'){
    nB$nBetas <- 7
    nB$nBetasBNT <- 2
    nB$nBetasATS <- 1
    nB$nBetasSigma <- 3
  }

  if(speciesGr == 'bnt'){
    nB$nBetas <- 7
    nB$nBetasBNT <- 2
    nB$nBetasATS <- 1
    nB$nBetasSigma <- 3
  }

  if(speciesGr == 'ats'){
    nB$nBetas <- 7
    nB$nBetasBNT <- 2
    nB$nBetasATS <- 1
    nB$nBetasSigma <- 3
  }

  dG[[ii]] <- ddG[[ii]][[1]] %>% runGrowthModel_Nimble(mcmcInfo,code,speciesGr,nB)

  #########################################
  # Nimble model run

  rm(Rmodel); rm(conf); rm(Rmcmc); rm(Cmodel); rm(Cmcmc); rm(mcmc); rm(mcmcProcessed); rm(obsPred)

  start <- Sys.time()
  print(start)

  Rmodel <- nimbleModel(dG[[ii]]$code, dG[[ii]]$constants, dG[[ii]]$data, dG[[ii]]$inits)
  #Rmodel$lengthDATA <- zoo::na.approx(dG[[ii]]$data$lengthDATA) ### check this is ok

  conf <- configureMCMC(Rmodel)
  conf$setThin(mcmcInfo$thinRate)
  conf$addMonitors(dG[[ii]]$params)

####################################### Added with Daniel 6/26/18 ######

  # sample sigma variables on log scale. Helps with very small values
   sigmaNodes <- c("grIntSigma",'sigmaIntSigma','grBetaSigma','sigmaBetaSigma')
   conf$removeSamplers(sigmaNodes)
   for( node in Rmodel$expandNodeNames(sigmaNodes) ){
     conf$addSampler(node, "RW", control = list(logScale = TRUE))
   }

  # remove samplers for last obs for each fish, for speed
#  nodeNames <- paste0("lengthDATA[", dG[[1]]$constants$lastObsRows, "]")
#  conf$removeSamplers(nodeNames)


  # Block sample - check cor, and block sample correlated vars
#e.g.  corNodeNames <- c('grBeta[1,2,3,4]','grBeta[2,2,3,4]')
#  conf$removeSamplers(corNodeNames)
#  conf$addSampler(corNodeNames, 'RW_block')


########################################

  Rmcmc <- buildMCMC(conf)#, enableWAIC = TRUE)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  mcmc <- runMCMC(Cmcmc,
                  nburnin = mcmcInfo$nBurnIn,
                  niter = mcmcInfo$nIter,
                  nchains = mcmcInfo$nChains,
                  samples = TRUE, samplesAsCodaMCMC = TRUE,
                  summary = FALSE, WAIC = FALSE)

  done <- Sys.time()
  elapsed[[ii]] <- done - start
  print(paste("Elapsed =", round(as.numeric(elapsed),2)))

  #########################################
  # Use jagsUI functions to get output into jagsUI format for analysis functions

  source("./R/jagsUIFunctions.R")
  mcmcProcessed <- process.output(mcmc, DIC=FALSE, params.omit=FALSE)#,params.omit=("lengthExp")) #jagsUI function, DIC = FALSE because it requires 'deviance'

  # this saves /figures/obsPred_speciesGr.png
  obsPred <- getRMSE_Nimble(mcmcProcessed,20)
  obsPred$rmse
  #obsPred$outliers%>% as.data.frame()

  getPropOverlap0Betas(mcmcProcessed,nB)

  #hist(mcmcProcessed$mean$grIndRE[mcmcProcessed$mean$grIndRE != 0] ,breaks=1000)
  mcmcProcessed$mean$sigmaIndRE

  # rHats
  mcmcProcessed$Rhat$grIntMu
  mcmcProcessed$Rhat$grIntSigma
  mcmcProcessed$Rhat$grInt

  mcmcProcessed$Rhat$grBetaMu
  mcmcProcessed$Rhat$grBetaSigma
  mcmcProcessed$Rhat$grBeta
##
  mcmcProcessed$Rhat$sigmaIntMu
  mcmcProcessed$Rhat$sigmaIntSigma
  mcmcProcessed$Rhat$sigmaInt

  mcmcProcessed$Rhat$sigmaBetaMu
  mcmcProcessed$Rhat$sigmaBetaSigma
  mcmcProcessed$Rhat$sigmaBeta
##
  mcmcProcessed$Rhat$sigmaIndRE
  plot(mcmcProcessed$sims.list$sigmaIndRE)

  #########################################
  # save data to file
  save(mcmcInfo,mcmcProcessed,dG,ddG,nB,speciesGr,riverOrderedIn, file = paste0('./data/out/dG_', as.integer(Sys.time()),'_crossVal',runCrossValidation_TF,"_", modelName, '_Nimble.RData'))

  ######################
  # Plot traces

  source("./R/plotFunctions.R") #will put into functions later
  plotInt_Nimble(mcmcProcessed)

  plotBetas_Nimble(mcmcProcessed,1:4)
  plotBetas_Nimble(mcmcProcessed,5:7)

  plotBetasBNT_Nimble(mcmcProcessed,1:3)
  plotBetasATS_Nimble(mcmcProcessed,1:2)


  plotSigmaInt_Nimble(mcmcProcessed)
  plotSigmaBetas_Nimble(mcmcProcessed,1:3)

  #plotSigmaBetasBNT_Nimble(mcmcProcessed,1)
  #plotSigmaBetasATS_Nimble(mcmcProcessed,1)

  #########################
  # Predictions
  ### temporary variables needed if do predictions after loading output files
#  ii <- 1
#  riverOrderedIn <- factor(c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)
#  speciesGr <- 'ats'
  ###

  limits <- 1.5 # -/+ limits on standardized range of input variable
  nPoints <- 5

  nItersForPred <- 100
  itersForPred <- sample( 1:mcmcInfo$nSamples, nItersForPred )

  # predictions across the grid
  p <- getPrediction( mcmcProcessed, limits, nPoints, itersForPred, dG[[1]]$constants, ddG[[ii]][[1]]$sampleIntervalMean,
                      c("len","temp", "flow","cBKT", "cBNT", "cATS") )#######################
  #p$variable <- "mean"
  # predictions of sigma across the grid

  #### need to comment out betas for ATS
  pSigma <- getPredictionSigma( mcmcProcessed, limits, nPoints, itersForPred, dG[[1]]$constants, ddG[[ii]][[1]]$sampleIntervalMean,
                                c("temp", "flow") )

  #pSigma$variable <- "sd"
  # combine predicted growth rates and sigma in predicted growth rates to get cv in predicted growth rates
  pBoth <- left_join( p,pSigma, by = c('flow','temp','iter','isYOY','season','river')) %>%
    mutate(predCV = predGrSigma/predGr)
  #pBoth$variable <- "cv"
  pBoth$predCV <- ifelse(pBoth$predCV > 10, NA, pBoth$predCV) # to avoid very large CVs when gr is very small
  pBoth$predCV <- ifelse(pBoth$predCV < 0, NA, pBoth$predCV) # to avoid negative CVs when gr is very small and negative

  save(pBoth,file = paste0("./data/out/P_ForMike_",speciesGr,".RData"))

  #
  # graph function, in analyzeOutputFunctions.R

  plotPred(pBoth, "predGr", "len", 1, speciesGr)
  plotPred(pBoth, "predGr", "temp", 1, speciesGr)
  plotPred(pBoth, "predGr", "flow", 1, speciesGr)
  plotPred(pBoth, "predGr", "cBKT", 1, speciesGr)
  plotPred(pBoth, "predGr", "cBNT", 1, speciesGr)
  plotPred(pBoth, "predGr", "cATS", 1, speciesGr)
  plotPred(pBoth, "predGr", c("temp", "flow"), 1, speciesGr)
  plotPred(pBoth, "predGr", c("temp","cBKT"), 1, speciesGr)
  plotPred(pBoth, "predGr", c("flow","cBKT"), 1, speciesGr)
  plotPred(pBoth, "predGr", c("cBKT","cBNT"), 1, speciesGr)

  plotPred(pBoth, "predGr", c("len","cBNT"), 1, speciesGr)
  #######################
  # graph function, in analyzeOutputFunctions.R

#  plotPred(pBoth, "predGrSigma", "len", 1, speciesGr) #shows means
  plotPred(pBoth, "predGrSigma", "temp", 1, speciesGr)
  plotPred(pBoth, "predGrSigma", "flow", 1, speciesGr)
  plotPred(pBoth, "predGrSigma", c("temp", "flow"), 1, speciesGr)

  plotPred(pBoth, "predCV", "len", 1, speciesGr)
  plotPred(pBoth, "predCV", "temp", 1, speciesGr)
  plotPred(pBoth, "predCV", "flow", 1, speciesGr)
  plotPred(pBoth, "predCV", c("temp", "flow"), 1, speciesGr)


#To do:


