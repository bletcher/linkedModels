
library(plyr)
library(rjags)
library(ggplot2)
library(abind)
library(parallel)
#rjags::load.module("dic")
load.module("dic")
#load.module("glm")

dMData$length[dMData$tagNumberCH=='1BF1FF6207' & dMData$season == 3 & dMData$year == 2005] <- NA
dMData$length[dMData$tagNumberCH=='1BF1FF6521' & dMData$season == 3 & dMData$year == 2005] <- NA
dMData$length[dMData$tagNumberCH=='1BF18CE7ED' & dMData$season == 2 & dMData$year == 2006] <- NA
dMData$length[dMData$tagNumberCH=='1BF20FF1B9' & dMData$season == 3 & dMData$year == 2005] <- NA
dMData$length[dMData$tagNumberCH=='257C67CA48' ] <- NA
dMData$length[dMData$tagNumberCH=='1BF20EB7A4' & dMData$season == 4 & dMData$year == 2008] <- NA
dMData$length[dMData$tagNumberCH=='00088D215D' & dMData$season == 4 & dMData$year == 2010] <- NA



dMData$riverOrdered <- factor(dMData$river,levels=c('WEST BROOK','WB JIMMY','WB MITCHELL','WB OBEAR'), ordered=T)

# means for standardizing
#####################################################################
stdBySeasonRiverYear <- ddply( dMData, .(riverOrdered,riverN,season,year), summarise,
                               lengthMn=mean(length, na.rm=TRUE),
                               tempMn=mean(fullMeanT, na.rm=TRUE),
                               tempMnP=mean(temperatureForP, na.rm=TRUE),
                               flowMn=mean(fullMeanD, na.rm=TRUE),
                               dischMnP=mean(dischargeForP,na.rm=T) )
stdBySeasonRiverYear<-stdBySeasonRiverYear[!is.na(stdBySeasonRiverYear$riverN),]


stdBySeasonRiver <- ddply( stdBySeasonRiverYear, .(riverOrdered,riverN,season), summarise,
                           lengthMean=mean(lengthMn, na.rm=TRUE),
                           lengthSd=sd(lengthMn, na.rm=TRUE),
                           lengthLo = quantile(lengthMn,c(0.025), na.rm=TRUE),
                           lengthHi = quantile(lengthMn,c(0.975), na.rm=TRUE),
                           tempMean=mean(tempMn, na.rm=TRUE),
                           tempMeanP=mean(tempMnP, na.rm=TRUE),
                           tempSd=sd(tempMn, na.rm=TRUE),
                           tempSdP=sd(tempMnP, na.rm=TRUE),
                           tempLo = quantile(tempMn,c(0.025), na.rm=TRUE),
                           tempHi = quantile(tempMn,c(0.975), na.rm=TRUE),
                           flowMean=mean(flowMn, na.rm=TRUE),
                           flowSd=sd(flowMn, na.rm=TRUE),
                           dischMeanP=mean(dischMnP,na.rm=T),
                           dischSdP=sd(dischMnP,na.rm=T),
                           flowLo = quantile(flowMn,c(0.025), na.rm=TRUE),
                           flowHi = quantile(flowMn,c(0.975), na.rm=TRUE) )
############# To get rid of NA Rivers
stdBySeasonRiver<-stdBySeasonRiver[!is.na(stdBySeasonRiver$riverN),]

#####################################################################
stdBySeason <- ddply( dMData, .(season), summarise,
                      lengthMean=mean(length, na.rm=TRUE),
                      lengthSd=sd(length, na.rm=TRUE),
                      lengthLo = quantile(length,c(0.025), na.rm=TRUE),
                      lengthHi = quantile(length,c(0.975), na.rm=TRUE),
                      tempMean=mean(fullMeanT, na.rm=TRUE),
                      tempMeanP=mean(temperatureForP, na.rm=TRUE),
                      tempSd=sd(fullMeanT, na.rm=TRUE),
                      tempSdP=sd(temperatureForP, na.rm=TRUE),
                      tempLo = quantile(fullMeanT,c(0.025), na.rm=TRUE),
                      tempHi = quantile(fullMeanT,c(0.975), na.rm=TRUE),
                      flowMean=mean(fullMeanD, na.rm=TRUE),
                      flowSd=sd(fullMeanD, na.rm=TRUE),
                      dischMeanP=mean(dischargeForP,na.rm=T),
                      dischSdP=sd(dischargeForP,na.rm=T),
                      flowLo = quantile(fullMeanD,c(0.025), na.rm=TRUE),
                      flowHi = quantile(fullMeanD,c(0.975), na.rm=TRUE) )

# standardize by river  - for age0 fall lengths
stdByRiver <- ddply( dMData, .(riverOrdered,riverN), summarise,
                     lengthSd0 = sd  (subset( length, age == 0 & season == 3 ), na.rm=TRUE),
                     lengthMean0 = mean(subset( length, age == 0 & season == 3 ), na.rm=TRUE) )

stdByRiver <- stdByRiver[!is.na(stdByRiver$riverN),]
stdByRiver$river <- as.numeric(stdByRiver$riverOrdered)

#stdBySeasonRiver<-rbind(stdBySeasonRiver,c('zRiv1','0',rep(NA,ncol(stdBySeasonRiver)-2)))

#####
# fdDATA is flood and drought frequencies and durations
fdDATA$year <- as.numeric( fdDATA$year )
fdDATA$year2 <- fdDATA$year
fdDATA$year <- fdDATA$year-min(fdDATA$year) + 1

floodDur <- matrix(0,max(fdDATA$season),max(fdDATA$year))
droughtDur <- matrix(0,max(fdDATA$season),max(fdDATA$year))
floodFreq <- matrix(0,max(fdDATA$season),max(fdDATA$year))
for ( i in 1:nrow(fdDATA) ){
  floodDur[fdDATA$season[i],fdDATA$year[i]] <- fdDATA$floodDur[i]
  droughtDur[fdDATA$season[i],fdDATA$year[i]] <- fdDATA$droughtDur[i]
  floodFreq[fdDATA$season[i],fdDATA$year[i]] <- fdDATA$floodFreq[i]

}
#####


# function to add dummy rows and columns for zRiv=1
addRowColMeans <- function(m){
  m <- cbind( rowMeans(m),m )
  m <- rbind( colMeans(m),m )
  return ( m )
}
# function to add dummy columns for zRiv=1
addColMeans <- function(m){
  m <- cbind( rowMeans(m),m )
  return ( m )
}

meanForTotalN <- matrix(NA,nrow=4,ncol=5)
sdForTotalN <- matrix(NA,nrow=4,ncol=5)
for(s in 1:4){
  for(r in 1:5){
    meanForTotalN[s,r] <- mean(estTotalN[s,,r])
    sdForTotalN[s,r] <- sd(estTotalN[s,,r])
  }
}

meanForN <- array(NA,dim=c(4,5,2))
sdForN <- array(NA,dim=c(4,5,2))
for(s in 1:4){
  for(r in 1:5){
    for(yoy in 1:2){
      meanForN[s,r,yoy] <- mean(estN[s,,r,yoy])
      sdForN[s,r,yoy] <- sd(estN[s,,r,yoy])
    }
  }
}

meanForBNTTotalN <- matrix(NA,nrow=4,ncol=5)
sdForBNTTotalN <- matrix(NA,nrow=4,ncol=5)
for(s in 1:4){
  for(r in 1:5){
    meanForBNTTotalN[s,r] <- mean(estBNTTotalN[s,,r])
    sdForBNTTotalN[s,r] <- sd(estBNTTotalN[s,,r])
  }
}

meanForBNTN <- array(NA,dim=c(4,5,2))
sdForBNTN <- array(NA,dim=c(4,5,2))
for(s in 1:4){
  for(r in 1:5){
    for(yoy in 1:2){
      meanForBNTN[s,r,yoy] <- mean(estBNTN[s,,r,yoy])
      sdForBNTN[s,r,yoy] <- sd(estBNTN[s,,r,yoy])
    }
  }
}
countForAllTotalN <- matrix(NA,nrow=4,ncol=5)
meanForAllTotalN <- matrix(NA,nrow=4,ncol=5)
sdForAllTotalN <- matrix(NA,nrow=4,ncol=5)
for(s in 1:4){
  for(r in 1:5){
    countForAllTotalN[s,r] <- (estTotalN[s,,r] + estBNTTotalN[s,,r])
    meanForAllTotalN[s,r] <- mean(estTotalN[s,,r] + estBNTTotalN[s,,r])
    sdForAllTotalN[s,r] <- sd(estTotalN[s,,r] + estBNTTotalN[s,,r])
  }
}

countForAllN <- array(NA,dim=c(4,5,2))
meanForAllN <- array(NA,dim=c(4,5,2))
sdForAllN <- array(NA,dim=c(4,5,2))
for(s in 1:4){
  for(r in 1:5){
    for(yoy in 1:2){
      countForAllN[s,r] <- (estN[s,,r,yoy] + estBNTN[s,,r,yoy])
      meanForAllN[s,r,yoy] <- mean(estN[s,,r,yoy] + estBNTN[s,,r,yoy])
      sdForAllN[s,r,yoy] <- sd(estN[s,,r,yoy] + estBNTN[s,,r,yoy])
    }
  }
}

tempForN<- array(NA,dim=c(4,5,max(dMData$year-min(dMData$year) + 1)))
for(s in 1:4){
  for(y in 1:max(dMData$year-min(dMData$year) + 1)){
    tempForN[s,1,y]<-0
    for(r in 1:4){
      tempForN[s,r+1,y]<-(mean(dMData$fullMeanT[dMData$season==s&as.numeric(dMData$riverOrdered)==r&(dMData$year-min(dMData$year) + 1)==y],na.rm=T)
                          - stdBySeasonRiver$tempMean[ 4*(r-1)+s ] ) / stdBySeasonRiver$tempSd[ 4*(r-1)+s ]
      if(tempForN[s,r+1,y]=='NaN')tempForN[s,r+1,y]<-(stdBySeasonRiver$tempMean[4*(r-1)+s]- stdBySeasonRiver$tempMean[ 4*(r-1)+s] ) / stdBySeasonRiver$tempSd[ 4*(r-1)+s ]
    }
  }
}


flowForN<- array(NA,dim=c(4,5,max(dMData$year-min(dMData$year) + 1)))
for(s in 1:4){
  for(y in 1:max(dMData$year-min(dMData$year) + 1)){
    flowForN[s,1,y]<-0
    for(r in 1:4){
      flowForN[s,r+1,y]<-(mean(dMData$fullMeanD[dMData$season==s&as.numeric(dMData$riverOrdered)==r&(dMData$year-min(dMData$year) + 1)==y],na.rm=T)
                          - stdBySeasonRiver$flowMean[4*(r-1)+s] ) / stdBySeasonRiver$flowSd[4*(r-1)+s]
      if(flowForN[s,r+1,y]=='NaN')flowForN[s,r+1,y]<-(stdBySeasonRiver$flowMean[4*(r-1)+s]- stdBySeasonRiver$flowMean[4*(r-1)+s] ) / stdBySeasonRiver$flowSd[4*(r-1)+s]
    }
  }
}




############ Predictors that are in a matrix have season in rows and river in columns
d <- within(
  data = list(),
  expr = {

    encDATA = as.numeric(dMData$enc) #$msEnc
    riverDATA = as.numeric(dMData$riverOrdered) #-3
    nRivers = length(unique(dMData$riverN))-1 #may need to add one for unobs
    lengthDATA = dMData$length
    availableDATA = dMData$available01
    ind = as.numeric(dMData$tagNumberCH)
    # For standardizing length
    lengthMean = addColMeans( matrix(stdBySeasonRiver$lengthMean,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    lengthSd =   addColMeans( matrix(stdBySeasonRiver$lengthSd,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )

    lengthMean0 = stdByRiver$lengthMean0
    lengthSd0 = stdByRiver$lengthSd0
    # environmental covariates pertaining to intervals.  These are
    # covariates of growth and survival

    # For standardizing env predictors of growth and surv
    tempMean = addColMeans( matrix(stdBySeasonRiver$tempMean,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    tempSd =   addColMeans( matrix(stdBySeasonRiver$tempSd,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    flowMean = addColMeans( matrix(stdBySeasonRiver$flowMean,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    flowSd =   addColMeans( matrix(stdBySeasonRiver$flowSd,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )

    ## Predictors of phi for correcting N1 where countForN ==0
    tempForN = tempForN
    flowForN = flowForN

    # not used anymore, see below
    #  tempDATA = ( as.numeric(dMData$fullMeanT) - stdBySeason$tempMean[ as.numeric(dMData$season)] ) / stdBySeason$tempSd[ as.numeric(dMData$season) ]
    #  flowDATA = ( as.numeric(dMData$fullMeanD) - stdBySeason$flowMean[ as.numeric(dMData$season)] ) / stdBySeason$flowSd[ as.numeric(dMData$season) ]


    # Now doing standardization in the bugs code to make sure that unobserved fish get observed std env data
    tempDATA = as.numeric(dMData$fullMeanT)
    flowDATA = as.numeric(dMData$fullMeanD)

    flowMeanDATA = flowMean
    flowSDDATA =   flowSd
    tempMeanDATA = tempMean
    tempSDDATA =   tempSd

    # emPermNA, used to censor likelihood for permanent emigrants
    # 1 on line before last observation with subsequent bottom of the study site antenna hit. 0's before and after if em, NAs otherwise
    # trying emPerm without the NAs
    emPermDATA = dMData$emPerm

    intervalDays = as.numeric(dMData$fullMeanIntLen )

    # Environmental covariates for p
    flowP = as.numeric(dMData$dischargeForP)
    temperatureP = as.numeric(dMData$temperatureForP)
    #For standardizing env predictors of p
    flowMeanP = addRowColMeans( matrix(stdBySeasonRiver$dischMeanP,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    flowSdP = addRowColMeans( matrix(stdBySeasonRiver$dischSdP,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    tempMeanP = addRowColMeans( matrix(stdBySeasonRiver$tempMeanP,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )
    tempSdP = addRowColMeans( matrix(stdBySeasonRiver$tempSdP,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1) )

    # , growthSd = sd(((dMData$lagLength - dMData$length)/(as.numeric(dMData$intervalLength)))*365/4, na.rm=TRUE)
    ######## NEVER!!!! #########  gr = (dMData$lagLength - dMData$length)/(as.numeric(dMData$intervalLength))
    # indexing of the input and state vectors
    year = dMData$year-min(dMData$year) + 1
    nYears = max(dMData$year)-min(dMData$year)+1
    season = as.numeric(as.character(dMData$season))
    nAllRows = length(dMData[,1])
    nFirstObsRows = evalList$nFirstObsRows
    firstObsRows = evalList$firstObsRows
    nOcc = length(unique(dMData$sampleNum))
    occ = dMData$sampleNum-min(dMData$sampleNum)-1
    nEvalRows = evalList$nEvalRows # rows that will matter if we start using JS, and
    evalRows = evalList$evalRows   # that matter now for the growth model
    lastPossibleRows = subset( 1:nrow(dMData),dMData$lastAIS==dMData$ageInSamples ) # need to put this in makedMData
    nLastPossibleRows = evalList$nFirstObsRows

    lastObsRows = evalList$lastObsRows
    nLastObsRows = evalList$nLastObsRows

    lastRows = lastPossibleRows
    nLastRows = nLastPossibleRows

    nOut = evalList$nEvalRows # evalRows to output for each trace

    #create variables that hold information on counts - data held in statsForN (made in makeDMData.R - based on pheno2Long, so has all cohorts. need to throw away years before dMData's first cohort)
    minYear <- min(dMData$year)
    firstYearIndex <- minYear-statsForN$minYear + 1
    # countForN has dummy river 1 in it

    countForTotalNDATA <- estTotalN
    meanForTotalN <- meanForTotalN
    sdForTotalN <- sdForTotalN

    countForNDATA <- estN
    meanForN <- meanForN
    sdForN <- sdForN

    countForBNTTotalNDATA <- estBNTTotalN
    meanForBNTTotalN <- meanForBNTTotalN
    sdForBNTTotalN <- sdForBNTTotalN

    countForBNTNDATA <- estBNTN
    meanForBNTN <- meanForBNTN
    sdForBNTN <- sdForBNTN

    countForAllTotalNDATA <- countForAllTotalN
    meanForAllTotalN <- meanForAllTotalN
    sdForAllTotalN <- sdForAllTotalN

    countForAllNDATA <- countForAllN
    meanForAllN <- meanForAllN
    sdForAllN <- sdForAllN

    #  dMDataF <- dMData[ dMData$first == dMData$sampleNum, ]
    #  nTagged1 <- table(dMDataF$season,dMDataF$year,dMDataF$riverOrdered)

    #Fill in random #s for zRiv=1
    #  nTagged <- abind(matrix(round(runif(4*nYears,10,50)), nrow=4,ncol=nYears),nTagged1)
    floodDurDATA <- floodDur
    droughtDurDATA <- droughtDur
    floodFreqDATA <- floodFreq

    # mean intervaldays by season and river for interval boundaries [ s,r ]
    dIntDays <- data.frame(enc=encDATA, int=intervalDays, river=riverDATA, season=season)
    dIntDaysMean <- ddply(dIntDays[dIntDays$enc==1,], .(season,river), colMeans)
    intervalMeans <- addColMeans( matrix(dIntDaysMean$int,nrow=length(unique(dMData$season)),ncol=length(unique(as.numeric(dMData$riverN)-0))-1, byrow=T) )
    rm(dIntDays)

    # propSamp - proportion of each season,river,year combo that got sampled (proportion of sctions sampled)
    # this doesn't work indexed by evalRows becasue there's no river value when fish aren't captured
    #propSampledDATA = dMData$propSampled

    # check data for propSampledDATA
    #tmp <- data.frame(cbind(dMData$year,dMData$season,dMData$river,dMData$enc,dMData$section))
    #names(tmp ) <- c('year','season','river','enc','section')
    #tmp2 <- tmp[tmp$enc==1 ,]  #& tmp$section %in% 1:47
    #f <- ftable(tmp2$year,tmp2$river,tmp2$season,tmp2$section)
    #cbind(1:nrow(f),rowSums(f))

    propSampledDATA =  array( 1, c(4,5,11) ) #season, river year

    propSampledDATA[ c(1,4),2:5,1 ] <- 0     #all spring and winter samples in 2002
    propSampledDATA[ 2,3:4,1 ] <- 0          #J and M summer samples in 2002
    #propSampledDATA[ 4,2:5,1 ] <- 0          #all winter samples in 2002
    propSampledDATA[ 4,2,2 ] <- 30/47        #WB winter sample in 2003
    propSampledDATA[ 4,2,3 ] <- 3/47         #WB winter sample in 2004
    propSampledDATA[ 4,2,4 ] <- 0            #WB winter sample in 2005
    propSampledDATA[ 4,2,6 ] <- 0            #WB winter sample in 2007


    # zeroSectionsDATA - completely unsampled winter samples. not including samples before season 4,year 1 because we didn't want to rewrite the meanPhiS34 indexing [mostly noise ni these estimates anyway]

    zeroSectionsDATA <- array( 0, c(4,5,11) ) #season, river year

    zeroSectionsDATA[ 3:4,2:5,1 ] <- 1          #all winter samples in 2002
    zeroSectionsDATA[ 3:4,2,4 ] <- 1            #WB winter sample in 2005
    zeroSectionsDATA[ 3:4,2,6 ] <- 1            #WB winter sample in 2007



    ####################################
    # create a data frame of max lengths for YOYs from Matt with some visual fixes,

    # including fall,winter 0+ fish in YOY category
    #  cutoffYOYDATA <- cutoffYOYDATA

    # including fall,winter 0+ and spring 1+ fish in YOY category
    cutoffYOYDATA <- cutoffYOYInclSpring1DATA
    #
    ########################################
  }
)


# function to make initial z matrix, with 1s when known alive and NAs otherwise
encInit<-function(sN, first, last){#, river){
  z.iv <- array(NA, dim=length(first))
  z.iv[(sN>first)&(sN<=(last))] <- 1 #z[first] gets set to 1 in bugs code-new version of jags doesn't like
  return(z.iv)
}


emPermInit <- function(e){
  eOut <- array(NA, dim=length(e))
  eOut <- ifelse( is.na(e), 0, e )
  return(eOut)
}

encInitMS<-function(sN, first, last, river){
  for (i in 1:(length(first))){
    river[i] <- river[i] - 0
    if ((sN[i] >= first[i]) & (sN[i] <= (last[i]))) {
      if( is.na(river[i]) ) river[i] <- river[i-1]
    }
    else river[i] <- NA
  }

  for (i in 1:(length(first))){
    if(sN[i] == first[i]) river[i] <- NA
  }
  return(river + 1)
}

inits <- function(){

  within(
    data = list(),
    expr = {

      ####Need to provide random initial values for at least one param besides z, otherwise all chains are the same

      z = as.numeric(encInit(dMData$sampleNum, dMData$first, dMData$last))
      #       isYOY1 = (dMData$length > 90) +1 # could make this s,r,y-specific
      #       zRiv = as.numeric(encInitMS(dMData$sampleNum, dMData$first, dMData$last, as.numeric(dMData$riverOrdered)))

      #Changed from 4 to add DD parm
      grBeta = array( runif(5*2*4*5,-2,2),c(5,2,4,5))  #array(rnorm(11*2*4*5),c(11,2,4,5))
      grBeta[ ,1,1:2, ] <- 0
      grSigmaBeta = array( abs(runif(2*4*5*11,0,5) ),c(2,4,5,11)) #T[0,]

      muGrSigmaBeta = array( abs(runif(2*4*5,0,5) ),c(2,4,5))
      sigmaGrSigmaBeta= array(runif(2*4*5,0,5),c(2,4,5))

      # let Jags assign initials because priors are season-dependent and we're lazy
      #grBetaInt = array( ( rnorm(2*11*4*5) ), c(2,4,5,11) )
      #muGrBetaInt= array( (rnorm(2*4*5) ),c(2,4,5))
      #    sigmaGrBetaInt= array(runif(2*4*5,0,10),c(2,4,5))

      muPhiBeta = array( (runif(5*2, -0.5,0.5) ),c(5,2))
      sigmaPhiBeta = array(runif(5*2,0,1),c(5,2))

      pBetaInt = array(runif(2*4*11*5, -2.5, 2),c(2,4,11,5))    #array(rnorm(5*4*5),c(5,4,5))    #Includes is YOY
      pBeta =    array(runif(2*4*11*5, -0.5, 0.5),c(2,4,11,5))    #array(rnorm(5*4*5),c(5,4,5))    #Includes is YOY
      muPBeta = array( (runif(2, -0.5,0.5) ),c(2))
      sigmaPBeta = array(runif(2,0,1),c(2))
      #Changed from 4 to add DD parm
      phiBeta = array(runif(5*2*4*5, -0.5, 0.5),c(5,2,4,5))
      phiBetaInt = array( ( runif(2*4*5, -0.5, 0.5) ), c(2,4,5) )
      psiBeta = array(runif(4*4*4,-5.5,-2.5),c(4,4,4))

      # grSigmaBeta[1,,] = runif(4)

    })

}

# MCMC settings
nAdapt <- 500 #500

nIter <- 500000 #25000                      # total num of iters
nIterChunk <- 100 #round(nIter/nSave)   # num of iters per chunk
nSave <- round(nIter/nIterChunk)                       # num of saves

nThin <- 5
nChains <- 1

#top('Pause')

varsToMonitor<-c(


  'pBeta'
  , 'pBetaInt'
  , 'muPBeta'
  , 'sigmaPBeta'

  , 'phiBeta'
  #  , 'phiBetaSize'
  , 'phiBetaInt'
  #  , 'meanPhiS34'
  , 'sigmaPhiBeta'
  , 'muPhiBeta'

  #  , 'muPhiYear'
  #  , 'sigmaPhiYear'

  #  , 'xi.phi'
  #  , 'mu.phi.raw'
  #  , 'sd.phi'


  , 'psiBeta'

  , 'grBeta'
  , 'muGrBeta'
  , 'sigmaGrBeta'

  , 'grSigmaBeta'
  , 'muGrSigmaBeta'
  , 'sigmaGrSigmaBeta'

  , 'grBetaInt'
  , 'muGrBetaInt'
  , 'sigmaGrBetaInt'


  , 'deviance'  # take out if dic module is not loaded

  #  , 'zOut'

  #  , 'lengthOut'
  #  , 'stdLengthOut'
  #  , 'zRivOut'

  #  , 'emPermDATAOut'
  #  , 'probEmPerm'

  #  , 'grLengthOut'
  #  , 'grTempOut'
  #  , 'grFlowOut'
  #  , 'grNOut'

  #  ,'phiLengthOut'
  #  ,'phiTempOut'
  #  ,'phiFlowOut'

  #  ,'pLengthOut'
  #  ,'pFlowOut'
  #   , 'stdN'
  #   , 'N1'
  #   , 'meanN'
  #   , 'sdN'

)

# out1=out; fileOutName <- 'outMSRiver2.RData'  ## for running a second set of iters



cl <- makeCluster(5)                                # Request 3 cores - equivalent to # of chains
clusterExport(cl, c("d", "inits", "varsToMonitor", 'bugsName', 'nChains', 'nAdapt','nThin', 'nIterChunk',
                    'encInit', 'dMData','phiBetaIntMeans','abind')) # Make these available -need to send any functions used in inits()
clusterSetRNGStream(cl = cl)#, 4345)


( beforeJags <- Sys.time() )
print( beforeJags )
print( 'Adapting, no progress bar available when running chains in parallel' )
beforeAdapt <- Sys.time()
###################################
# adaptation
####################################

adapt <- clusterEvalQ(cl, {
  library(rjags)

  jm <- jags.model(
    file = bugsName,
    data = d,
    inits = inits,
    n.chains = nChains,
    n.adapt = nAdapt,
  )
  return( jm )
})

afterAdapt <- Sys.time()
print(paste('adaptation time:',afterAdapt - beforeAdapt,sep=" "))
adaptTime = afterAdapt - beforeAdapt

####################################
# sampling
####################################

# run 1 iter just to get structure of outP
outP <- clusterEvalQ(cl, {

  library(rjags)
  load.module("dic")

  js <- jags.samples(  ##############coda.samples( #
    model = jm,
    variable.names = varsToMonitor,
    n.iter = 1,
    thin = nThin,
    progress.bar = 'text'
  )
  return( js )  ###############as.mcmc
})


print( paste( Sys.time(), Sys.time() - beforeJags, sep=" " ) )

####################################
# run separate chunks of iterations, nIterChunk at a time
# each chunk gets appended to outP and saved as 'fileOutName'
####################################
sampleTime <- array( NA,nSave )

for( i in 1:nSave ){

  holdTime <- Sys.time()
  print( paste( 'Before js, i = ',i,'; Start time = ',holdTime, sep="" ) )

  outP[[i]] <- clusterEvalQ(cl, {

    library(rjags)
    load.module("dic")

    js <- jags.samples( ##############coda.samples( #
      model = jm,
      variable.names = varsToMonitor,
      n.iter = nIterChunk,
      thin = nThin,
      progress.bar = 'text'
    )

    return( (js) )  ##############as.mcmc
  })

  print( paste( 'After js, Loop = ',i,' out of ',nSave, ', Saved iters done = ',nIterChunk*i, sep="") )
  sampleTime[i] = Sys.time() - holdTime
  print( paste( 'Chunk end time =',Sys.time(), '; Chunk run time =', round( Sys.time() - holdTime, 2 ), sep=" " ) )
  print( getwd() )
  print( paste( 'mean sampleTime = ', round(mean(sampleTime, na.rm=T),2), '; adaptTime = ', round( adaptTime,2 ), sep='' ) )

  print( "##########" )

  save(d, outP, i, adaptTime, sampleTime, file = fileOutName)
}


( done <- Sys.time() )

print(afterAdapt - beforeAdapt)
print(done - beforeJags)
