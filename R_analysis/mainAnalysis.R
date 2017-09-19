#install.packages("devtools")
#devtools::install_github('bletcher/linkedModels')
library(tidyverse)
library(linkedModels)
library(jagsUI)

dr <- "west"
cd <- getCoreData(dr) %>%
        cleanData(dr) %>%
        mergeSites(dr) %>%
        mutate(drainage = dr)

save(cd, file = paste0('/home/ben/linkedModels/data/cd_',dr,'.RData'))


# specific analyses
# limit to bkt
if (!exists("cd")) load(paste0('/home/ben/linkedModels/data/cd_',dr,'.RData'))

#################################
# Movement model
# To start, will just assume fish don't move
# until they are seen somewhere else.
# Accomplished with %>% fillSizeLocation(size = F) in getCoreData()

#################################
# Growth model
#
#
#
start <- Sys.time()
ddd <- cd %>%
         filter(species == "bkt", cohort %in% c(2006)) %>% #,distMoved < 48, distMoved > 0, enc == 1) %>%
         prepareDataForJags() %>%
         runGrowthModel()
elapsed <- Sys.time() - start
