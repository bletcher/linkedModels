#install.packages("devtools")
#devtools::install_github('bletcher/linkedModels')
library(linkedModels)

dr <- "west"
cd <- getCoreData(dr) %>%
        cleanData(.,dr) %>%
        mutate(drainage = dr)

sites <- getSites(dr)
# merge in riverMeter for sections
cd <- left_join(cd, sites, by = c("river","section","area"))
cd$riverMeter <- ifelse( cd$survey == "shock" | cd$survey == "portableAntenna", cd$river_meter, cd$riverMeter )

save(cd, file = '/home/ben/linkedModels/data/cdWB.RData')


# specific analyses
# limit to bkt
dat <- cd %>% filter(species == "bkt",distMoved < 48, distMoved > 0, enc == 1)
