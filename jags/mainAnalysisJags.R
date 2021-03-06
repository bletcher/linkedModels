# next -
# lit review on soze-dep growth
# take out phi

# compare countPStd (by species) with countPAllSppStd (all spp in analysis)

# species-, age, location-specific rates of size-dep growth
# estimate spp-sp abundances for each river for each sample
# size effect ~ 1/spp/river/season/year
# run models:
# size
# size*spp
# size*spp*river
# size*spp*river*season
# then add in variable for 'harshness', either production (biomass) over interval or abundance
# (overall and spp-specific) or env (flow, temp). Test how size-dep gr varies across harshness for best model from moedl selection.


#install.packages("devtools")
# library(devtools)
#devtools::install_github('bletcher/linkedModels')
#devtools::install_github(repo = "Conte-Ecology/westBrookData", subdir = "getWBData")
library(arm)
library(zoo)
library(arrayhelpers)
library(linkedModels)
library(jagsUI)
library(getWBData)
library(stringr)
library(lubridate)
library(dbplyr)
library(tidyverse)
######################################################
# selection criteria

reconnect()

drainage <- "west" # ==

speciesDet <- c("bkt", "bnt","ats") #keep as all three spp
speciesInDet <- factor(speciesDet, levels = c('bkt','bnt','ats'), ordered = T)

speciesGr <- "bkt"
speciesInGr <- factor(speciesGr, levels = c('bkt','bnt','ats'), ordered = T)

riverOrderedIn <- factor(c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)

minCohort <- 2010#1995 # >=
maxSampleInterval <- 200 # <
runDetectionModelTF <- F
runCrossValidationTF <- F
percentLeftOut <- 10

meanOrIter <- "mean"; iter <- 1
####### or ########
#meanOrIter = "iter"

#make sure species are always in order and indexed correctly for arrays
#riverOrderedIn <- factor(1:3,levels=c('mainstem', 'west', 'east'),labels = c('mainstem', 'west', 'east'), ordered=T)

# update yoy cutoffs as get new data using
# getYOYCutoffs(cd,drainage) - make sure to call with cd so all data enter and minYear==1997 for 'west'
# which is called within prepareDataForJags() and is saved in
# getAndPrepareDataWB.R
# will need to update code for stanley


######################################################
# Get data from database
# Only need to run this when data or functions have changed
# or you want to change drainages

cdFileBeforeDetMod <- paste0('./data/cd_',drainage,"_",minCohort,'_BeforeDetMod.RData')

if ( file.exists(cdFileBeforeDetMod) ) {
  load(cdFileBeforeDetMod)
} else {

  # core data
  cdBeforeDetMod <- getCoreData(drainage) %>%
    cleanData(drainage) %>%
    mergeSites(drainage) %>%
    mutate(drainage = drainage,
           countP = NA) %>% # placeholder so prepareDataForJags() works for detection model
    addRawCounts(drainage, filteredAreas = c("inside","trib")) %>% # counts and summed masses of all fish, tagged and untagged. Not adjusted for P (happens in addDensities())
   # do this in detection model, removeUnsampledRows(drainage, removeIncomplete = T) %>% # removes enc=0 rows for samples with no or very few() sections sampled (p=0 otherwise for those samples)
    removeLowAbundanceRivers(drainage) # removes ats,jimmy fish and bnt,mitchell from drainage=='west' - too few fish for estimates

  save(cdBeforeDetMod, file = cdFileBeforeDetMod)
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
           sampleInterval < maxSampleInterval %>%  # this removes the later yearly samples. Want to stick with seasonal samples
           removeUnsampledRows(drainage, removeIncomplete = T) # removes enc=0 rows for samples with no or very few() sections sampled (p=0 otherwise for those samples)
  )

dddD <- ddddD %>%
  prepareDataForJags('detection')

if (runDetectionModelTF) {
  ddD <- dddD[[1]] %>% runDetectionModel(parallel = TRUE)
  save(ddD,dddD,ddddD, file = paste0('./data/out/ddD_', dModelName,'.RData'))
  #whiskerplot(ddD, parameters="pBetaInt")

  done <- Sys.time()
  (elapsed <- done - start)

  # add counts derived from counts/p to cdBeforeDetMod. ddD and ddddD need to have all three spp to create spp-specific abundances
  cd <- adjustCounts( cdBeforeDetMod,ddD,ddddD,meanOrIter,iter )
  cdFile <- paste0('./data/cd_',drainage,"_",minCohort,'.RData')
  save(cd, file = cdFile)

} else {
  load( file = paste0('./data/out/ddD_', dModelName,'.RData') )
  load( file = paste0('./data/cd_',drainage,"_",minCohort,'.RData') )
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

##########################################################################################################################################
# Growth model
# One species at a time
###############
#
#
#speciesInGr <- factor(speciesGr, levels = c('bkt','bnt','ats'), ordered = T)

#cd <- cd %>% filter(!(tag %in% tmp$tag)) #### for bnt test

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

  # saving into a list for now, could also map()

  # dddG[[ii]] <- addBiomassDeltas( dddG[[ii]] )
  dddG[[ii]] <- addSurvivals( ddddG,ddD,meanOrIter,iter ) %>% removeFishWithManyIntermediateNAs()

  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth") # returns a list, list for running model in [[ii]][[1]], input df in [[ii]][[2]]
    #save(dG,ddG,dddG, file = paste0('./data/out/dG_', modelName, '_forLmer.RData')) # for Lmer model with all spp
  #####
  start <- Sys.time()
  print(start)
  print(c("in loop",meanOrIter,ii,iter))

  dG[[ii]] <- ddG[[ii]][[1]] %>% runGrowthModel( parallel = FALSE )

  save(dG,ddG,dddG, file = paste0('./data/out/dG_', as.integer(Sys.time()),'_', modelName, '.RData'))

  done <- Sys.time()
  elapsed[[ii]] <- done - start
  print(paste("Elapsed =",elapsed))

  obsPred <- getRMSE(0.6)
  obsPred$rmse
  obsPred$outliers%>% as.data.frame()
  #####
source("./R/plotFunctions.R") #will put into functions later
  plotInt()
  plotBetas(1:2)
  plotBetas(3:4)
  plotBetas(5:18)

  plotSigma()
  plotSigmaBetas(1:2)
  plotSigmaBetas(3:4)


##### }

ggplot(ddG[[1]][[2]], aes(nAllFishBySpeciesPStd,grLength)) + geom_point() + facet_grid(species+river~season)

whiskerplot(dG[[1]], parameters = "lengthExp[1:20]")
traceplot(dG[[1]], parameters = "lengthExp[7:8]")

whiskerplot(dG[[1]], parameters = "grInt")
whiskerplot(dG[[1]], parameters = "grIntMu")
whiskerplot(dG[[1]], parameters = "grSigma")
jagsUI::traceplot(dG[[1]], parameters = "grSigma")

jagsUI::traceplot(dG[[1]], parameters = "grInt[2,1,1,3]")
jagsUI::traceplot(dG[[1]], parameters = "grIntSigma")

whiskerplot(dG[[1]], parameters = "grSigma")
whiskerplot(dG[[1]], parameters = "grSigmaMu")
whiskerplot(dG[[1]], parameters = "sigmaBetaSigma")
whiskerplot(dG[[1]], parameters = "sigmaBeta")
traceplot(dG[[1]], parameters = "sigmaIntSigma")

whiskerplot(dG[[1]], parameters = "grBeta")
whiskerplot(dG[[1]], parameters = "grBeta[2,,,,]")
whiskerplot(dG[[1]], parameters = "grBetaMu")
whiskerplot(dG[[1]], parameters = "grBetaSigma")
jagsUI::traceplot(dG[[1]], parameters = "grBeta")

whiskerplot(dG[[1]], parameters = "grIndRE[1:100]")
jagsUI::traceplot(dG[[1]], parameters = "grIndRE[1:10]")

whiskerplot(dG[[1]], parameters = "grBetaMu")
jagsUI::traceplot(dG[[1]], parameters = "grBetaMu")


######################################
# Explore predictions

limits <- 1.5 # -/+ limits on standardized range of input variable
nPoints <- 5

nItersForPred <- 100
itersForPred <- sample( 1:dG[[1]]$mcmc.info$n.samples,nItersForPred )

# predictions across the grid
p <- getPrediction( dG[[1]], limits, nPoints, itersForPred, c("len", "temp", "flow","count") )
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

# predictions of sigma across the grid
pSigma <- getPredictionSigma( dG[[1]], limits, nPoints, itersForPred, c("len", "temp", "flow","count") )
#######################
# graph function, in analyzeOutputFunctions.R

plotPred(pSigma, "len", 1, "bkt")
plotPred(pSigma, "temp", 1, "bkt")
plotPred(pSigma, "flow", 1, "bkt")
plotPred(pSigma, "count", 1, "bkt")



pairs(ddG[[1]][[1]][c('countPStd','tempStd','flowStd')])


#To do:


