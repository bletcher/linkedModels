#install.packages("devtools")
#devtools::install_github('bletcher/linkedModels')
#install.packages("devtools")
#devtools::install_github('Conte-Ecology/getWBData')
library(tidyverse)
library(linkedModels)
library(jagsUI)
library(getWBData)
library(lubridate)

dr <- "west"
cdFile <- paste0('./data/cd_',dr,'.RData')

if ( file.exists(cdFile) ) {
  load(cdFile)
} else {
  cd <- getCoreData(dr) %>%
    cleanData(dr) %>%
    mergeSites(dr) %>%
    mutate(drainage = dr)

  save(cd, file = cdFile)
}

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

propSampled <- 1
start <- Sys.time()
ddd <- cd %>%
  filter(  species == "bkt",
           cohort %in% c(2006),
           tag %in% sample(unique(tag), propSampled*length(unique(tag)))
  )  %>% #,distMoved < 48, distMoved > 0, enc == 1)
  prepareDataForJags()

dd <- ddd %>% runGrowthModel()
done <- Sys.time()
(elapsed <- done - start)

st <- 9355
end <- 9370
data.frame(
  occ = c(st:end),
  ddd$lengthDATA[st:end],
  ddd$ind[st:end],
  ddd$occ[st:end],
  ddd$season[st:end],
  ddd$riverDATA[st:end]
)

#To do:
# Line up ddd$YOY with actual year in model run
# add intervalMeans back in to gr[]
