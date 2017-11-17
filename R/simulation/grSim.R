
library(tidyverse)
library(jagsUI)

propMissing <- 0.5
numSamps <- 15
nInd <- 100
initSize <- 100

trimOcc <- TRUE
binomProb <- 0.45
  
simD <- data.frame(samp = rep(1:numSamps), ind = rep(1:nInd,each=numSamps))
simD$len <- ifelse( simD$samp == 1, round(rnorm(nInd,initSize,10)), NA )


gr <- function(l){
  out <- l*0.25 + rnorm(1,0,2)
  return(out)
}

if (trimOcc){
  lastOcc <- data.frame( ind = 1:nInd, lastOcc = rbinom(nInd,numSamps-1,binomProb) + 1 )
  simD <- left_join( simD,lastOcc ) %>% filter( samp <= lastOcc )
}
hist(simD$samp)


simD <- simD %>%
  group_by(ind) %>%
  mutate( first = min(samp),
          firstObs = ifelse( samp == first,TRUE,FALSE),
          last = max(samp),
          lastObs = ifelse( samp == last,TRUE,FALSE)) %>%
  ungroup()

for(i in 1:nrow(simD)){
  if (simD$firstObs[i] == 0) {
    simD$len[i] <- simD$len[ i-1 ] + gr(simD$len[ i-1 ])
  }
}

simD$random <- runif(nrow(simD))

simD$missing <- simD$firstObs == 0 & simD$random < propMissing
simD$lenMissing <- ifelse( simD$missing, NA, simD$len )


simD <- simD %>%
  mutate( meanLen = mean( lenMissing, na.rm=T ),
          sdLen = sd( lenMissing, na.rm=T ),
          stdLen = (lenMissing - meanLen)/sdLen
  )

# ggplot(simD, aes(samp,stdLen, color=factor(ind))) +
#   geom_point() +
#   geom_line()

str(simD)

evalRows <- which( !simD$lastObs )
firstObsRows <- which(simD$firstObs)
lastObsRows <- which(simD$lastObs)

simD$rowNumber <- 1:nrow(simD)
#ggplot(simD ,aes(samp,lenMissing)) +geom_line() + geom_point() +facet_wrap(~ind)

# simD$lenMissing[4] <- NA
# simD$lenMissing[5] <- NA
# simD$lenMissing[6] <- NA
simD$lenMissingStd <- (simD$lenMissing - mean(simD$lenMissing,na.rm=T))/sd(simD$lenMissing,na.rm=T)

dataIn <- list(
  nEvalRows = length(evalRows), evalRows = evalRows,
  nFirstObsRows = length(firstObsRows), firstObsRows = firstObsRows,
  nLastObsRows = length(lastObsRows), lastObsRows = lastObsRows,
  lengthMean = unique(simD$meanLen),
  lengthSD = unique(simD$sdLen),
  lengthDATA = simD$lenMissingStd,
  
  nInd = nInd,
  ind = simD$ind
  )

params <- c('grInt','grBeta','grSigma','obsErr','length')

start <- Sys.time()
print(start)

ss <- jags(data = dataIn,
          inits = NULL,
          parameters.to.save = params,
          model.file = "./grSim.txt",
          n.chains = 3,
          n.adapt = 5000, #1000
          n.iter = 5000,
          n.burnin = 3000,
          n.thin = 4,
          parallel = T
)

done <- Sys.time()
elapsed <- done - start
print(paste("Elapsed =",elapsed))

whiskerplot(ss, parameters = "length[1:80]")
traceplot(ss, parameters = "length[3:7]")
traceplot(ss, parameters = "grInt")
traceplot(ss, parameters = "grBeta")
traceplot(ss, parameters = "grSigma")
traceplot(ss, parameters = "obsErr")

