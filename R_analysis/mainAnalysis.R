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
library(tidyverse)
######################################################
# selection criteria

drainage <- "west" # ==

species <- c("bkt", "bnt","ats") #
speciesIn <- factor(species, levels = c('bkt','bnt','ats'), ordered = T)
riverOrderedIn <- factor(c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)

minCohort <- 1997#1995 # >=
maxSampleInterval <- 200 # <
runDetectionModelTF <- T
runCrossValidationTF <- F
percentLeftOut <- 10

reconnect()

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

cdFile <- paste0('./data/cd_',drainage,"_",paste0(paste0(species,collapse = ''),minCohort),'.RData')

if ( file.exists(cdFile) ) {
  load(cdFile)
} else {

  # core data
  cd <- getCoreData(drainage) %>%
    cleanData(drainage) %>%
    mergeSites(drainage) %>%
    mutate(drainage = drainage,
           countP = NA) %>% # placeholder so prepareDataForJags() works for detection model
    addCounts(drainage) # counts and summed masses of all fish, tagged and untagged. Not adjusted for P (happens in addDensities())

  save(cd, file = cdFile)
}



#################################
# Detection model
##################
#
#
dModelName <- paste0(paste0(species,collapse = ''),minCohort)

(start <- Sys.time())

ddddD <- cd %>%
  filter(  species %in% speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval  # this removes the later yearly samples. Want to stick with seasonal samples
  )

dddD <- ddddD %>% prepareDataForJags('detection')

if (runDetectionModelTF) {
  ddD <- dddD[[1]] %>% runDetectionModel(parallel = TRUE)
  save(ddD,dddD,ddddD, file = paste0('./data/out/ddD_', dModelName,'.RData'))
  #whiskerplot(ddD, parameters="pBetaInt")
}
done <- Sys.time()
(elapsed <- done - start)

# with full intercepts, effect of nPasses is insignificant (massive overlap). Leaving, nPasses out. Just going with pBetaInt.
#whiskerplot(ddD, parameters = "pBetaInt")
#whiskerplot(ddD, parameters = "pBeta")

#################################
# Movement model
#################
# To start, will just assume fish don't move
# until they are seen somewhere else.
# Accomplished with %>% fillSizeLocation(size = F) in getCoreData()

#################################
# Growth model
###############
#
#

ddddG <- cd %>%
  filter(  species %in% speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval, # this removes the later yearly samples. Want to stick with seasonal samples
           #enc == 1 #change to btw first and last so we can impute missing obs?
           knownZ == 1
  )
nSeasons <- n_distinct(ddddG$season, na.rm=T)
load(file = paste0('./data/out/ddD_', dModelName,'.RData'))

# merge in density data, either overall means or single iterations at a time
# meanOrIter ="mean" uses means of all iterations
# meanOrIter ="iter" uses the sample that is the combo of sampleToUse and chainToUse (ignored if meanOrIter="mean")
#  for numOfItersToUse samples

######################################
######################################

meanOrIter <- "mean"
####### or ########
#meanOrIter = "iter"

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

modelName <- paste0(paste0(species,collapse = ''),minCohort)

dddG <- list()
ddG <- list()
dG <- list()
print(elapsed)
elapsed <- list()
# run the growth model for detection model iterations in itersToUse
ii <- 0
for (iter in itersToUse) {
  ii <- ii + 1

  # saving into a list for now, could also map()

  dddG[[ii]] <- addDensityData( ddddG,ddD,ddddD,meanOrIter,iter )
  # dddG[[ii]] <- addBiomassDeltas( dddG[[ii]] )
  dddG[[ii]] <- addSurvivals( dddG[[ii]],ddD,meanOrIter,iter )
  #dddG[[ii]] <- crossValidate( dddG[[ii]],runCrossValidationTF ) # Moved inside prepareDataforJags()

  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth") # returns a list, list for running model in [[1]], input df in [[2]]

  #####
  start <- Sys.time()
  print(start)
  print(c("in loop",meanOrIter,ii,iter))

  dG[[ii]] <- ddG[[ii]][[1]] %>% runGrowthModel( parallel = TRUE )

  done <- Sys.time()
  elapsed[[ii]] <- done - start
  print(paste("Elapsed =",elapsed))

  getRMSE()
  #####

  plotInt()
  plotBetas(1:2)

  plotSigma()
  plotSigmaBetas(1:2)
  plotSigmaBetas(3:4)

  save(dG,ddG,dddG, file = paste0('./data/out/dG_', as.integer(Sys.time()),'_', modelName, '.RData'))
}

whiskerplot(dG[[1]], parameters = "lengthExp[1:20]")
traceplot(dG[[1]], parameters = "lengthExp[7:8]")

whiskerplot(dG[[1]], parameters = "grInt")
whiskerplot(dG[[1]], parameters = "grIntMu")
whiskerplot(dG[[1]], parameters = "grSigma")
jagsUI::traceplot(dG[[1]], parameters = "grSigma")

jagsUI::traceplot(dG[[1]], parameters = "grInt[2,1,1,3]")
jagsUI::traceplot(dG[[1]], parameters = "grIntSigma")

whiskerplot(dG[[1]], parameters = "sigmaInt")
whiskerplot(dG[[1]], parameters = "sigmaBetaMu")
whiskerplot(dG[[1]], parameters = "sigmaBetaSigma")
whiskerplot(dG[[1]], parameters = "sigmaBeta[4,,,,]")
traceplot(dG[[1]], parameters = "sigmaIntSigma")

whiskerplot(dG[[1]], parameters = "grBeta")
whiskerplot(dG[[1]], parameters = "grBeta[2,,,,]")
whiskerplot(dG[[1]], parameters = "grBetaMu")
whiskerplot(dG[[1]], parameters = "grBetaSigma")
jagsUI::traceplot(dG[[1]], parameters = "grBeta")

whiskerplot(dG[[1]], parameters = "grIndRE[1:10]")
jagsUI::traceplot(dG[[1]], parameters = "grIndRE[1:10]")

whiskerplot(dG[[1]], parameters = "grBetaMu")
jagsUI::traceplot(dG[[1]], parameters = "grBetaMu")


########################################################################################
# traceplots

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotInt <- function(){
  ggGrInt <- array2df(dG[[1]]$sims.list$grInt, label.x = "est")

  ggGrInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-1.5,1.5) + facet_grid(d2+d4~d5+d3)
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotBetas <- function(b){
  ggGrBeta <- array2df(dG[[1]]$sims.list$grBeta, label.x = "est")

  ggGrBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

  gg <- list()
  numBetas <- 18
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + ylim(-1,1) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigma <- function(){
  ggSigmaInt <- array2df(dG[[1]]$sims.list$grSigma, label.x = "est")

  ggSigmaInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggSigmaInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggSigmaInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + facet_grid(d2+d4~d5+d3)
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigmaBetas <- function(b){
  ggGrBeta <- array2df(dG[[1]]$sims.list$sigmaBeta, label.x = "est")

  ggGrBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

  gg <- list()
  numBetas <- 4
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + ylim(-1,1) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}

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

plotPred(p, "len", 1, "bkt")
plotPred(p, "temp", 1, "bkt")
plotPred(p, "flow", 1, "bkt")
plotPred(p, "count", 1, "bnt")
plotPred(p, c("temp", "flow"), 1, "bkt")
plotPred(p, c("temp", "flow"), 0, "bkt")
plotPred(p, c("temp","count"), 1, "bkt")
plotPred(p, c("flow","count"), 1, "bkt")
plotPred(p, c("len","count"), 1, "bkt")




########################################################################################
# traceplots

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
ggGrInt <-array2df(dG[[1]]$sims.list$grInt, label.x = "est")

ggGrInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggGrInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + facet_grid(d2+d4~d5+d3)+ylim(-1,1)

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
ggSigmaInt <- array2df(dG[[1]]$sims.list$sigmaInt, label.x = "est")

ggSigmaInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggSigmaInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggplot(filter(ggSigmaInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + facet_grid(d2+d4~d5+d3)

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
ggGrBeta <- array2df(dG[[1]]$sims.list$grBeta, label.x = "est")

ggGrBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggGrBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

gg <- list()
numBetas <- 11
for (i in 1:numBetas){
  gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + ylim(-1,1) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
}


ggSigmaBeta <- array2df(dG[[1]]$sims.list$sigmaBeta, label.x = "est")

ggSigmaBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggSigmaBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

ggSigma <- list()
numBetas <- 7
for (i in 1:numBetas){
  ggSigma[[i]] <- ggplot(filter(ggSigmaBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
}

# isYOY[ evalRows[i] ],season[ evalRows[i] ] ]
# [1:675, 1:2, 1:4]
gg2 <- array2df(dG[[1]]$sims.list$grIntMu, label.x = "est")
gg2 <- array2df(dG[[1]]$sims.list$grIntSigma, label.x = "est")
gg2 <- array2df(dG[[1]]$sims.list$sigmaIntMu, label.x = "est")
gg2 <- array2df(dG[[1]]$sims.list$sigmaIntSigma, label.x = "est")

gg2$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
gg2$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggplot(gg2, aes(iter,est)) + geom_point( aes(color = factor(chain)), size = 0.1 ) + facet_grid(d2~d3)

gg3 <- array2df(dG[[1]]$sims.list$grBetaMu, label.x = "est")
gg3 <- array2df(dG[[1]]$sims.list$grBetaSigma, label.x = "est") # seems uninformative

gg3$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
gg3$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggplot(gg3, aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + facet_grid(d2~d4+d3)


pairs(ddG[[1]][c('countPStd','tempStd','flowStd')])


#To do:
# Line up ddd$YOY with actual year in model run

