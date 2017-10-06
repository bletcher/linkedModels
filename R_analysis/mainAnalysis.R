# next -


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
#devtools::install_github('bletcher/linkedModels')
#install.packages("devtools")
#devtools::install_github('Conte-Ecology/getWBData')
library(arm)
library(arrayhelpers)
library(linkedModels)
library(jagsUI)
library(getWBData)
library(lubridate)
library(tidyverse)
######################################################
# selection criteria

drainage <- "west" # ==
species <- c("bkt", "bnt") #"bkt" # ==
minCohort <- 2002 # >=
maxSampleInterval <- 200 # <
runDetectionModelTF <- TRUE

#make sure species are always in order and indexed correctly for arrays
speciesIn <- factor(species, levels = c('bkt','bnt','ats'), ordered = T)
riverOrderedIn <- factor(c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T)
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

cdFile <- paste0('./data/cd_',drainage,'.RData')

if ( file.exists(cdFile) ) {
  load(cdFile)
} else {
  cd <- getCoreData(drainage) %>%
    cleanData(drainage) %>%
    mergeSites(drainage) %>%
    mutate(drainage = drainage,
           countP = NA) # placeholder so prepareDataForJags() works for detection model

  cd <- addNPasses(cd,drainage)

  save(cd, file = cdFile)
}


#################################
# Detection model
##################
#
#

(start <- Sys.time())
ddddD <- cd %>%
  filter(  species %in% speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval # this removes the later yearly samples. Want to stick with seasonal samples
  )

dddD <- ddddD %>% prepareDataForJags('detection')

if (runDetectionModelTF) {
  ddD <- dddD %>% runDetectionModel(parallel = TRUE)
  save(ddD, file = './data/out/ddD.RData')
}
done <- Sys.time()
(elapsed <- done - start)


#whiskerplot(ddD, parameters = "pBetaInt")

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

# need to check grBeta estimates

ddddG <- cd %>%
  filter(  species %in% speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval, # this removes the later yearly samples. Want to stick with seasonal samples
           enc == 1
  )

load(file = './data/out/ddD.RData')

# merge in density data, either overall means or single iterations at a time
# meanOrIter ="mean" uses means of all iterations
# meanOrIter ="iter" uses the sample that is the combo of sampleToUse and chainToUse (ignored if meanOrIter="mean")
#  for numOfItersToUse samples

######################################
######################################

meanOrIter = "mean"
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

dddG <- list()
ddG <- list()
dG <- list()
elapsed <- list()
# run the growth model for detection model iterations in itersToUse
ii <- 0
for (iter in itersToUse) {
  ii <- ii + 1
  start <- Sys.time()
  print(c("in loop",Sys.time(),meanOrIter,ii,iter))

  # saving into a list for now, could also map()

  dddG[[ii]] <- addDensityData( ddddG,ddD,ddddD,meanOrIter,iter )

  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth")

  dG[[ii]] <- ddG[[ii]] %>% runGrowthModel( parallel = TRUE )

  done <- Sys.time()
  elapsed[[ii]] <- done - start
  print(paste("Elapsed =",elapsed))
  save(dG, file = './data/out/dG.RData')
}

######################################

# Explore predictions
limits <- 2
nPoints <- 5
iterForPred <- 100
preds <- getPrediction(dG[[1]], limits, nPoints, iterForPred)

p <- preds %>% filter(flow == 0,temp==0,isYOY==1,count==0)

ggplot(p, aes(len,predGr,color=species)) +
  geom_point() +
 # geom_line() +
  geom_smooth(method="lm") +
  facet_grid(river~season)


# 2, isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
whiskerplot(dG[[1]], parameters = "grBeta[4,2,,,]")


#  [isYOY,species,season,riverDATA]
whiskerplot(dG[[1]], parameters = "muGrBetaInt[2,,1,]")
#  [1:2,isYOY,species,season,riverDATA]
whiskerplot(dG[[1]], parameters = "grBeta[2,,,,]")
#  [isYOY,species,season,riverDATA,year]
whiskerplot(dG[[1]], parameters = "grBetaInt[2,,,1]")


whiskerplot(dG[[1]], parameters = "muGrBeta[,]")



ggplot(filter(dddG[[1]], grLength<0.95 & grLength > (-0.5)), aes(countPStd,grLength)) + geom_point() +
#ggplot(filter(dddG[[1]], grLength<0.95 & grLength > (-0.5)), aes(countPAllSppStd,grLength)) + geom_point() +
  geom_smooth(method='lm') +
  facet_grid(riverOrdered~season)

ggplot(filter(dddG[[1]], grLength<0.95 & grLength > (-0.5)), aes(countPAllSppStd,countPStd)) + geom_point() +
  geom_smooth(method='lm') +
  facet_grid(riverOrdered~season)

pairs(ddG[[1]][c('countPStd','tempStd','flowStd')])


#To do:
# Line up ddd$YOY with actual year in model run

