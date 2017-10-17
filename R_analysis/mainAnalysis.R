# next -

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
runDetectionModelTF <- FALSE

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
print(elapsed)
elapsed <- list()
# run the growth model for detection model iterations in itersToUse
ii <- 0
for (iter in itersToUse) {
  ii <- ii + 1
  start <- Sys.time()
  print(start)
  print(c("in loop",meanOrIter,ii,iter))

  # saving into a list for now, could also map()

  dddG[[ii]] <- addDensityData( ddddG,ddD,ddddD,meanOrIter,iter )
  dddG[[ii]] <- addBiomassDeltas( dddG[[ii]] )
  dddG[[ii]] <- addSurvivals( dddG[[ii]],ddD,meanOrIter,iter )

  ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth")

  dG[[ii]] <- ddG[[ii]] %>% runGrowthModel( parallel = TRUE )

  done <- Sys.time()
  elapsed[[ii]] <- done - start
  print(paste("Elapsed =",elapsed))
  save(dG, file = paste0('./data/out/dG_', as.integer(Sys.time()), '.RData'))
}

######################################

# Explore predictions
limits <- 2 # -/+ limits on standatrdized range of input variable
nPoints <- 5

nItersForPred <- 100
itersForPred <- sample( 1:dG[[1]]$mcmc.info$n.samples,nItersForPred )
preds <- getPrediction( dG[[1]], limits, nPoints, itersForPred )
#predsSigma <- getPredictionSigma( dG[[1]], limits, nPoints, itersForPred )

#######################
# length graph
p <- preds %>% filter(flow == 0,temp==0,isYOY==1,count==0,species=="bkt")

ggplot(p, aes(len,predGr,group = iter)) +
#  geom_point() +
  geom_line( color="lightgrey", alpha=0.25 ) +
 # geom_smooth(method="lm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  facet_grid(river~season)

####################################################################################
# len by biomass graph
yoyIn <- 0
speciesIn <- "bkt"
p <- preds %>% filter(isYOY==yoyIn,species==speciesIn,biomass %in% c(-2,2))
p$iterBiomass <- paste0(p$iter,p$biomass)

ggplot(p, aes(len,predGr,group = iterBiomass)) +
  geom_line( aes(color=biomass), alpha=0.25 ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  ggtitle(paste0('isYOY = ',yoyIn,', species = ',speciesIn)) +
  facet_grid(river~season)
####################################################################################

# flow graph
p <- preds %>% filter(len == 0,temp==0,isYOY==1,count==0,species=="bkt")
#p <- predsSigma %>% filter(temp==0,isYOY==1,count==0,species=="bkt")

ggplot(p, aes(flow,predGr,group=iter)) +
  #  geom_point() +
  geom_line( color="lightgrey", alpha=0.25 ) +
  # geom_smooth(method="lm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  facet_grid(river~season)

# temperature*flow graph
p <- preds %>% filter(len == 0,isYOY==1,count==0,species=="bkt")
#p <- predsSigma %>% filter(isYOY==1,count==0,species=="bkt")
p$iterFlow <- paste0(p$iter,p$flow)

ggplot(p, aes(temp,predGr,group=iterFlow, color=flow)) +
  #  geom_point() +
  geom_line( alpha=0.25 ) +
  # geom_smooth(method="lm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  facet_grid(river~season)

# count graph
p <- preds %>% filter(len == 0,temp==0,isYOY==1,flow==0,species=="bkt")
p <- preds %>% filter(temp==0,isYOY==1,flow==0,species=="bkt")
p$iterLen <- paste0(p$iter,p$len)

#p <- predsSigma %>% filter(temp==0,isYOY==1,flow==0,species=="bkt")

ggplot(p, aes(count,predGr,group = iterLen,color = (len))) +
  #  geom_point() +
  geom_line(   alpha=0.25 ) +
  # geom_smooth(method="lm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  facet_grid(river~season)



  meanGRs <- ddddG %>%
    group_by(species,riverOrdered,season) %>%
    summarize(mean=mean(grLength*sampleInterval,na.rm=T),
              int=mean(sampleInterval,na.rm=T))

########################################################################################
# traceplots

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
ggGrInt <-array2df(dG[[1]]$sims.list$grInt, label.x = "est")

ggGrInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggGrInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + facet_grid(d2+d4~d5+d3)

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
numBetas <- 5
for (i in 1:numBetas){
  gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
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

###########################################################
# [yoy,s]
whiskerplot(dG[[1]], parameters = "grIntMu[,]")
traceplot(dG[[1]], parameters = "grIntMu")
whiskerplot(dG[[1]], parameters = "grIntSigma[,]")

# sigmaInt[ yoy,spp,s,r ]
whiskerplot(dG[[1]], parameters = "sigmaInt[2,1,,1]")

# sigmaIntMu[ yoy,s ]
whiskerplot(dG[[1]], parameters = "sigmaIntMu[,]")
whiskerplot(dG[[1]], parameters = "sigmaIntSigma[,]")

#  [isYOY,species,season,riverDATA]
whiskerplot(dG[[1]], parameters = "grInt[2,,,1]")

#grBetaMu[ i,yoy,s ]
whiskerplot(dG[[1]], parameters = "grBetaMu[,2,]")


# 2, isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
whiskerplot(dG[[1]], parameters = "grBeta[,2,,,]")
whiskerplot(dG[[1]], parameters = "sigmaBeta[,2,,,]")


#  [isYOY,species,season,riverDATA]
whiskerplot(dG[[1]], parameters = "muGrBetaInt[2,,4,]")
#  [1:2,isYOY,species,season,riverDATA]
whiskerplot(dG[[1]], parameters = "grBeta[3,2,1,1,]")


# yoy,spp,s,r,y
whiskerplot(dG[[1]], parameters = "grSigmaBeta[2,1,2,1,]")



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

