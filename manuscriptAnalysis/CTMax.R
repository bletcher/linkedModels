#From Evan's Ecology paper'
library(tidyverse)
library(ggplot2)

d <- data.frame(temp = seq(1,20.5,0.1))

# bkt first, bnt second

# Perf curve ###############
tOpt <- c(11.9,13.3)
cTMax <- c(19.5,19.3)
sigma <- c(3.7,4.1)

i <- 1
d$pT_BKT <- ifelse( d$temp <= tOpt[i],
                 exp(-((d$temp - tOpt[i])/(2*sigma[i]))^2), #there's an error in the paper
                 1 - ((d$temp - tOpt[i])/(tOpt[i] - cTMax[i]))^2
              )

i <- 2
d$pT_BNT <- ifelse( d$temp <= tOpt[i],
                    exp(-((d$temp - tOpt[i])/(2*sigma[i]))^2), #there's an error in the paper
                    1 - ((d$temp - tOpt[i])/(tOpt[i] - cTMax[i]))^2
)

ggplot(d, aes(temp,pT_BKT)) + geom_line() + geom_line(aes(temp,pT_BNT), color = "brown")
#################


getPerf <- function(tempIn,sppIn=1,ctMaxIn=1,sigmaIn=1){

  perf <- ifelse( tempIn <= tOpt[sppIn],
                      exp(-((tempIn - tOpt[sppIn])/(2*sigma[sppIn]))^2), #there's an error in the paper
                      1 - ((tempIn - tOpt[sppIn])/(tOpt[sppIn] - cTMax[sppIn]))^2
  )
  return(perf)
}

## gOpt part
lInf <- 354
k <- 0.0013
flowBeta <- 0.098
bktBeta <- -0.051
bntBeta <- -0.039

bktDen <- 7.01
bntDen <- 13.00
flow <- 0.14

len <- 100

getGOpt <- function(lenIn,flowIn=flow,bktDenIn=bktDen,bntDenIn=bntDen){

  gOpt <- k*(lInf-lenIn) + flowBeta*flowIn + bktBeta*bktDenIn + bntBeta*bntDenIn

  return(gOpt)
}

getGrowth <- function(lenIn,tempIn=10,intervalLength=70,flowIn=flow,bktDenIn=bktDen,bntDenIn=bntDen,sppIn=1,ctMaxIn=1,sigmaIn=1){

  growth <- getGOpt(lenIn) * getPerf(tempIn,sppIn,ctMaxIn=1,sigmaIn=1)
  len2 <- lenIn + growth * intervalLength

  return(list(growth=growth,len2=len2))
}

getGrowth(60:100,10)
getGrowth(100,10:18)




