
library(tidyverse)
library(jagsUI)



load('ddG_forSim.RData')
dd <- ddG[[1]]
load('dddG_forSim.RData')
ddd <- dddG[[1]]

nInd <- n_distinct(ddd$tagIndex, na.rm = T)
nOcc <- n_distinct(ddd$sampleIndex, na.rm = T)
minOcc <- min(ddd$sampleIndex)
nRivers <- n_distinct(ddd$riverOrdered, na.rm = T)
nSeasons <- n_distinct(ddd$season, na.rm = T)
nYears <- n_distinct(ddd$year, na.rm = T)
nSpecies <- n_distinct(ddd$species, na.rm = T)

ddd$riverN <- as.numeric(ddd$riverOrdered)
ddd$speciesN <- as.numeric(as.factor(ddd$species))

means <- ddd %>%
  group_by( speciesN,season,riverN ) %>%
  summarize( meanLen = mean(observedLength, na.rm = T),
             sdLen = sd(observedLength, na.rm = T),
             sampleIntervalMean = mean(sampleInterval, na.rm = T)
  )

ddd$observedLengthStd <- (ddd$observedLength - mean(ddd$observedLength,na.rm=T)) / sd(ddd$observedLength,na.rm=T)

ddd <- ddd %>%
  group_by(tag) %>%
  mutate( minOcc = min(sampleNumber), fOcc = (sampleNumber == minOcc)*1,
          maxOcc = max(sampleNumber), lOcc = (sampleNumber == maxOcc)*1,
          nObs = n()
  ) %>%
  ungroup()

ddd <- ddd %>% filter(tagIndex < 501) ##############################

evalRows2 <- which( ddd$lOcc == 0 )
firstObsRows2 <- which( ddd$fOcc == 1 )
lastObsRows2 <- which( ddd$lOcc == 1 )


nEvalRows2 <- length(evalRows2)
nFirstObsRows2 <- length(firstObsRows2)
nLastObsRows2 <- length(lastObsRows2)

nAllRows2 <- nEvalRows2 + nLastObsRows2

ddd$rowNumber <- 1:nrow(ddd)


dataIn2 <- list(
  nEvalRows = nEvalRows2, evalRows = evalRows2,
  nFirstObsRows = nFirstObsRows2, firstObsRows = firstObsRows2,
  nLastObsRows = nLastObsRows2, lastObsRows = lastObsRows2,
  lengthMean = array(means$meanLen, dim = c(nRivers,nSeasons,nSpecies)),
  lengthSD = array(means$sdLen, dim = c(nRivers,nSeasons,nSpecies)),
 # lengthMean = unique(simddd$meanLen),
#  lengthSD = unique(simddd$sdLen),
  lengthDATA = ddd$observedLengthStd,
  season = ddd$season,
  riverDATA = ddd$riverN,
  species = as.numeric(factor(ddd$species, levels = c('bkt','bnt','ats'), ordered = T)),
  ind = ddd$tagIndex,
  nInd = max(ddd$tagIndex) #length(unique(ddd$tagIndex))
)

params <- c('grBeta','grSigma','lengthDATA','lengthExp','grIntMean')

start <- Sys.time()
print(start)

s <- jags(data = dataIn2,
          inits = NULL,
          parameters.to.save = params,
          model.file = "./grSimWData3.txt",
          n.chains = 3,
          n.adapt = 5000, #1000
          n.iter = 2000,
          n.burnin = 500,
          n.thin = 4,
          parallel = T
)

done <- Sys.time()
elapsed <- done - start
print(paste("Elapsed =",elapsed))

whiskerplot(s, parameters = "lengthExp[1:80]")
whiskerplot(s, parameters = "lengthDATA[1:80]")
whiskerplot(s, parameters = "grInt[1:80]")
traceplot(s, parameters = "length[7]")
traceplot(s, parameters = "grInt")
traceplot(s, parameters = "grBeta")
traceplot(s, parameters = "grSigma")
traceplot(s, parameters = "obsErr")


ggplot(ddd,aes(sampleNumber,observedLength)) + geom_line() + geom_point() +facet_wrap(~tagIndex)
ggplot(ddd %>% filter(rowNumber %in% evalRows2),aes(ageInSamples,observedLength, color=factor(tagIndex))) +geom_line() + geom_point()+facet_wrap(~tagIndex)




ddd <- ddd %>%
  group_by(tag) %>%
  mutate( minOcc = min(sampleNumber), fOcc = (sampleNumber == minOcc)*1,
          maxOcc = max(sampleNumber), lOcc = (sampleNumber == maxOcc)*1
  ) %>%
  ungroup()

evalRowsS <- which( ddd$lOcc == 0 )
firstObsRowsS <- which( ddd$fOcc == 1 )
lastObsRowsS <- which( ddd$lOcc == 1 )
