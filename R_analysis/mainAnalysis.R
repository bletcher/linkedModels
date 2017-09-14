#install.packages("devtools")
#devtools::install_github('bletcher/linkedModels')
library(linkedModels)

dr <- "west"
cd <- getCoreData(dr) %>%
        cleanData(dr) %>%
        mergeSites(dr) %>%
        mutate(drainage = dr)

save(cd, file = paste0('/home/ben/linkedModels/data/cd_',dr,'.RData'))


# specific analyses
# limit to bkt
if (!exists("cd")) load(paste0('/home/ben/linkedModels/data/cd_',dr,'.RData'))
dat <- cd %>% filter(species == "bkt") #,distMoved < 48, distMoved > 0, enc == 1)

