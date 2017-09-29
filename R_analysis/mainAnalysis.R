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
    mutate(drainage = drainage)

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

whiskerplot(ddD, parameters = "pBetaInt")

den <- getDensities(ddddD, ddD, meanOrIter = "mean", sampleToUse = 30, chainToUse = 1 )
ggplot(filter(den,!is.na(countP)),aes(year,countP, color = species)) + geom_point() + geom_line() + facet_grid(riverOrdered~season)
ggplot(filter(den,!is.na(countP)),aes(count,countP, color = species)) + geom_point()
# save output, calc densities

st <- 160
end <- 180

data.frame(
  len = dddD$lengthDATA[st:end],
  sample = dddD$sample[st:end],
  ind = dddD$ind[st:end],
  enc = dddD$encDATA[st:end],
  zForInit = dddD$zForInit[st:end],
  season = dddD$season[st:end],
  river = dddD$riverDATA[st:end],
  year = dddD$year[st:end],
  yearForCutoff = dddD$yearForCutoff[st:end]
)

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

propSampled <- 1
(start <- Sys.time())
dddG <- cd %>%
  filter(  species == speciesIn,
           cohort >= minCohort,
           sampleInterval < maxSampleInterval, # this removes the later yearly samples. Want to stick with seasonal samples
           enc == 1
           # tag %in% sample(unique(tag), propSampled*length(unique(tag)))
  )  %>% #,distMoved < 48, distMoved > 0, enc == 1)
  prepareDataForJags()

ddG <- dddG %>% runGrowthModel(parallel = TRUE)
done <- Sys.time()
(elapsed <- done - start)

# get cutoffYOY right

whiskerplot(ddG, parameters = "muGrBetaInt")
traceplot(ddG, parameters = "muGrBetaInt")

whiskerplot(ddG, parameters = "muGrBeta")



st <- 9355
end <- 9370
data.frame(
  occ = c(st:end),
  dddG$lengthDATA[st:end],
  dddG$ind[st:end],
  dddG$occ[st:end],
  dddG$season[st:end],
  dddG$riverDATA[st:end]
)

#To do:
# Line up ddd$YOY with actual year in model run
# add intervalMeans back in to gr[]
