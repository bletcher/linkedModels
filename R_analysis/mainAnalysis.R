# next - add species to det and gr models (everywhere)


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
library(tidyverse)
library(linkedModels)
library(jagsUI)
library(getWBData)
library(lubridate)

######################################################
# selection criteria

drainage <- "west" # ==
speciesIn <- c("bkt", "bnt") #"bkt" # ==
minCohort <- 2002 # >=
maxSampleInterval <- 200 # <


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

dddD <- ddddD %>% prepareDataForJags()

ddD <- dddD %>% runDetectionModel(parallel = TRUE)
done <- Sys.time()
(elapsed <- done - start)

save(ddD, file = './data/out/ddD.RData')
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
# add species variable category

ddddG <- cd %>%
  filter(  species %in% speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval, # this removes the later yearly samples. Want to stick with seasonal samples
           enc == 1
  )

# merge in density data, either overall means or single iterations at a time
# meanOrIter ="mean" uses means of all iterations
# meanOrIter ="iter" uses the sample that is the combo of sampleToUse and chainToUse (ignored if meanOrIter="mean")
#  for numOfItersToUse samples

######################################
######################################

meanOrIter = "mean"
####### or ########
meanOrIter = "iter"

chainToUse <- 1
numItersToUse <- 2
if (meanOrIter == "iter") {
  itersToUse <- sort(sample(((dd$mcmc.info$n.samples/dd$mcmc.info$n.chains) * (chainToUse - 1)):
                            ((dd$mcmc.info$n.samples/dd$mcmc.info$n.chains) * (chainToUse - 0)),
                          numItersToUse))
} else {
  itersToUse == 1 # only run one loop over iter
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
  print(c("in loop",as.POSIXct(start),meanOrIter,ii,iter))

  # saving into a list for now, could also map()

  dddG[[ii]] <- addDensityData( ddddG,ddD,ddddD,meanOrIter,iter )

  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags()

  dG[[ii]] <- ddG[[ii]] %>% runGrowthModel( parallel = TRUE )

  done <- Sys.time()
  (elapsed[[ii]] <- done - start)
}





# get cutoffYOY right

whiskerplot(ddG, parameters = "muGrBetaInt")
traceplot(ddG, parameters = "muGrBetaInt")

whiskerplot(ddG, parameters = "muGrBeta")




#To do:
# Line up ddd$YOY with actual year in model run
# add intervalMeans back in to gr[]
