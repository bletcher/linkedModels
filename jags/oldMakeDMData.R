Skip to content
This repository
Search
Pull requests
Issues
Marketplace
Explore
@bletcher
Sign out
Unwatch 1
Star 0
Fork 0 bletcher/makeDMData
Code  Issues 0  Pull requests 0  Projects 0  Wiki  Settings Insights
Branch: master Find file Copy pathmakeDMData/makeDMData.R
4b14080  on Nov 6, 2013
@bletcher bletcher Revert "test3"
1 contributor
RawBlameHistory
1748 lines (1380 sloc)  78.5 KB
## as of 11/6/2013
# file location that was on felek got deleted during a switch
# saving this on git now.
# have not updated directories and file sources etc..







###

##########################################################################
# File location on felek: /var/lib/data/share
#
# Files to update on felek before running this file:
#
# fourRiversAllSpp.txt
# file created using access file query  fourRiversAllSpp:fourRivers query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp\fourRiversAllSpp.mdb
# export query as text file without formatting, include names on first row
# copy fourRiversAllSpp.txt to felek (home/ben/allSpp)

# adjSampNums.csv
# right now, median dates are based on WB only. probably should incude tribs in dates
#      make sure to udpate alignedDates.csv and medianSampleDates.csv in C:\PITTAGMAIN\CMR Analyses\adjustedDates
#        1) update by copying new sample rows from fourRiversAllSpp:Median Date for Each Sample Table and pasting into medianSampleDates.csv
#           will need to copy and paste and delete the river columns before pasting
#        2) rerun sampleDates.r
#        3) copy adjSampNums.csv to felek (home/ben/allSpp/adjustedDates)

# envData.txt
# file created using access file query  fourRiversAllSpp:envData query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
# export query as text file without formatting, include names on first row
# copy to felek home/ben/allSpp

# smolttrapFourRiversAllSpp.txt
# file created using access file query  fourRiversAts:smoltTrapfourRiversAts query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
#export query as text file without formatting

# antennaTribs.txt
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp fourRiversAllSpp:antennaTribs query
# export query as text file without formatting
# then copy to felek home/ben/allSpp

# antennaDuda.txt
# file created using access file query  fourRiversAts:antennaDuda query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
#export query as text file without formatting
# then copy to felek home/ben/allSpp

# yoySizeLimits
# file made with query of yoy size limits
# file created using access file query  fourRiversAllSpp:stockYearBySizeClasses query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp\fourRiversAllSpp.mdb
# run query [stockYearBySizeClasses], export and copy to
#'/data/projects/jwpcenter/bayesian/allSpp/stockYearBySizeClasses.txt'
# update any visual fixes to the size limits in callMSRiverP.R on Denver computer

# unsampled rivers
# update propSampledDATA and zeroSectionsDATA in callMSRiverP.R for unsampled rivers/years

#
#############################################################################

#############################################################################
# to do: code a newSection variable that includes above and below
#############################################################################

#############################################################################
# Things to be careful of
# 1. Make sure all the samples that need to included in dMdata are in adjDates

# 2. fill in yearSeasonListAnt for non-salmon samples
# 3. if have any samples that were not attempted (propSampledATSWest ==0), update the df 'temp' before it gets rbind to sampleDateList
# 4. Update new rivers (West and Shorey done) for propSampled when any sample was incomplete. This is just for ATS. need to ignorre or change for bkt
# 4. When writing code, make sure to order on tagNumber and sampleNumber after a merge if the order matters (e.g. evalRows)
#############################################################################

#############################################################################
# Major steps:
#   pheno: raw data from access
#     access query subsets on area = "inside" Or "trib" Or "above" Or "below", excludes antennas
#     basic data prep, cleans up inconsistencies, etc.
#     incoroprates adjDates from spreadsheet to line up sample nums among rivers
#
#   pheno2: [prepare to create long format]
#     1st step, subset pheno based on species and river
#
#   pheno2Long: [long format (augmented), one observation per row and one row
#               per sample (gets subsetted later)
#
#   dMData: long format, in addition to other subsetting, dMData only contains
#           available == 1 rows (truncated for fish that are too old - determined
#           by subsetDMdataAgeInSamples [ or fish that have emigrated permanently - no more like this ].
#           All environmental data get merged into dMData as well. Data for P are
#           on the occasion and data for phi and growth are over the interval.
#           Missing observations have either data from the average day during the
#           occ for P or data averaged (etc) over the median date interval.
#
#############################################################################

#############################################################################
#   Metadata for objects saved in the .RData output file

# "evalJSRows"          rows in dMData to evaluate when estimating phi in JS context. from ageInSamples==1 to subsetDMdataAgeInSamples
# "evalRows"            rows in dMData to evaluate when estimating phi in CJS context. from first to subsetDMdataAgeInSamples
# "firstObsRows"        rows in dMData that represent the first observation for each fish
# "lastObsRows"         rows in dMData that represent the last observation for each fish

# "nEvalJSRows"         number of rows in evalJSRows
# "nEvalRows"           number of rows in evalRows
# "nFirstObsRows"       number of rows in firstObsRows
# "nLastObsRows"        number of rows in lastObsRows

# "nonSummerObsRows"    rows in dMData that are not summer occasions (used for estimating p(maturation)
# "nNonSummerObsRows"   number of rows in nonSummerObsRows
# "nonSummerAIS"        ageInSample values that are not summer occasions
# "nSummerObsRows"      number of rows in nonSummerAIS
# "summerAIS"           ageInSample values that are summer occasions
# "summerObsRows"       number of rows in summerAIS

# "dMData"              This is the core data input file. it is vectorized so that each row represents a potential fish observation
#                       Around the actual observation occasoins for each fish, it is augmented backwards to ageInSamples == 1
#                       and forward to subsetDMdataAgeInSamples

# variables in dMData

# $ tagNumber             : Factor: PIT tag number
# $ ageInSamples          : num: increments by one for each combination of season and age. 1 = season 3 (Fall), age 0
# $ year                  : num:  year of sample occasion
# $ season                : Factor: season of sample occasion. 1=spring, 2=summer, 3=fall, 4=winter
# $ sampleNum             : num:  sample number occasion. Scaled to start with 1 for first occasion in subset of data in dMData.
#                                 May not match original sample numbers. yearSeasonList2 for original
#                                (for WB, adjusted for other rivers-see adjSampNums.csv) sample # (sampleNumAdj) and new sample # sampleNumAdjConsec
# $ drainage              : Factor: in pheno, "CATAMARAN" "SAWMILL"   "SHOREY"    "WEST"
# $ river                 : Factor: in pheno, "CATAMARAN BROOK" "SAWMILL RIVER" "SHOREY BROOK" "WB JIMMY" "WB MITCHELL" "WB OBEAR"  "WEST BROOK"
# $ date2                 : POSIXct: date of capture if captured
# $ length                : int: fork lenght in mm
# $ mature01              : num: numeric indicator of male parr maturation (0=not observed mature, 1=observed mature)
# $ section               : int: stream section number (usually 20-m lengths)
# $ enc                   : num: numeric indicator for capture (0=not captured, 1=captured)
# $ species               : chr: in pheno,  "ATS" "BKT" "BNT"
# $ age                   : Factor: age in years at capture - be careful, this is a factor so use as.numeric(as.character(pheno$age)) to get meaningful #s
# $ gtFirstOcc            : num: numeric indicator for an occasion greater than first capture (0 = first, 1 > first)
# $ first                 : int: sample number of first capture
# $ last                  : int:  sample number of last capture
# $ cohort                : num: year that a fish is estimated to be age-0
# $ lagDate2              : POSIXct: date of next capture occasion
# $ medianDate2           : POSIXct: median date of the sample occasion
# $ everEmWest            : num: numeric indicator for permanent emigration from the West Brook system (includes Jimmy and Mitchell).
#                                (1= yes, 0=no). Does not indcate occasion of last observation, just that a fish ever emigrated
# $ fromTo                : int: not very useful here. indicator of movement from one stream segement in the previous occasion to the current occ
# $ everMat01             : chr: indicator for whether a fish was ever observed as a mature male parr (1=yes, 0=no)
# $ area                  : Factor: location indicator, "ABOVE"  "BELOW"  "INSIDE" "TRIB", inisde is sections 1-47 in the WB
# $ riverM                : num: river meter from a fixed downstream location.
# $ tempDuringOccSampDays : num: average water temperature during the sampling occasion
# $ dischDuringOccSampDays: num: average water dishcarge during the sampling occasion
# $ temperature           : num: average water temperature on the day of capture
# $ discharge             : num: average water discharge on the day of capture
# $ temperatureForP       : num: 'temperature' if captured, 'tempDuringOccSampDays' if not captured
# $ dischargeForP         : num  'discharge' if captured, 'dischDuringOccSampDays' if not captured
# $ medianDateSampleNum   : POSIXct: median date of the sampling occasion
# $ propSampled           : num: proportion of the 47 sections that were sampled during occasions. Accounts for incomplete samples. Samples were all complete in the other rivers except for West Brook except for 2002 winter.
# $ propDaysOnAntEff      : num: proportion of days over the following interval we estimated (visually from a graph) that there were no reads on Duda's antennas - just for WEST now
# $ dateForEnv            : POSIXct: median sample date for non-capture occsasions, capture date (date2) for capture occasions. Used to fill in environmental data using addEnvironmentalData2
# $ lagDateForEnv         : POSIXct: next occasion
# $ matNA                 : num: numeric indicator for observed male parr maturity (1=yes, 0=no)
#                            1= matured in that maturity year (suummer-spring)
#                            0 = observed immature in fall
#                            NA = not observed mature and not captured in fall
# $ intervalLength        :Class 'difftime': interval length in days
# $ emPerm                : num: numeric indicator for occasion of permanent emigration (1=yes, 0=no)
# $ emPermNA              : num: numeric indicator for occasion of permanent emigration (1=yes, 0=no), with NAs after last occasion that a fish was observed
# $ maxT                  : num: maximun temperature during the interval following the current occasion
# $ meanT                 : num:  mean temperature during the interval following the current occasion
# $ medianT               : num:  median temperature during the interval following the current occasion
# $ minT                  : num:  minimum temperature during the interval following the current occasion
# $ sdT                   : num:  sd of temperature during the interval following the current occasion
# $ skewT                 : num:  skew of temperature during the interval following the current occasion
# $ maxD                  : num  maximun discharge during the interval following the current occasion
# $ meanD                 : num  mean discharge during the interval following the current occasion
# $ medianD               : num  median discharge during the interval following the current occasion
# $ minD                  : num  minimum discharge during the interval following the current occasion
# $ sdD                   : num  sd of discharge during the interval following the current occasion
# $ skewD                 : num  skew of discharge during the interval following the current occasion
#


#############################################################################


library(ecoPiffle)
library(lattice)
library(latticeExtra)
library(ggplot2)
#library(reshape)
library(moments)
#library(multicore) #linux only
library(date) #to get month/day from julian date using date.mdy

# data subsetting
#############################################################################

#choose a species
speciesSubset <- c('ATS')
speciesSubset <- 'BNT'

#############################################################################
if (speciesSubset == 'ATS'  ) {
  riverSubset <- c('WEST BROOK')#,'WB JIMMY') #,'WB MITCHELL',"WB OBEAR")
  #riverSubset <- c("SHOREY BROOK") # "SAWMILL RIVER","WB JIMMY","WB OBEAR", "WB MITCHELL", "CATAMARAN BROOK",
  #riverSubset <- c("SAWMILL RIVER")
  areaSubset <- c('INSIDE', 'ABOVE', 'BELOW')#, 'TRIB' )

  # subset variables for dMData
  # if include fish from cohorts < 1997, they will have ageInSamples that
  #don't go back to 1. for now we are leaving out cohort < 1997
  #when we want to include them, we'll need to augment back to ageInsamples 1
  #bay adding in negative smample numbers
  subsetDMdataCohortMin <- 1996 # >=
  subsetDMdataCohortMax <- 2050 # <=
  subsetDMdataAgeInSamples <- 15 # <

  # exclude fish that were captured for the first time after the following ageInSamples
  # set equal to subsetDMdataAgeInSamples - 1 to have no effect
  maxAgeInSamplesFirstCapt <- 4 #subsetDMdataAgeInSamples - 1 #4
}

if (speciesSubset == 'BKT' | speciesSubset == 'BNT'  ) {
  riverSubset <- c('WEST BROOK','WB JIMMY','WB MITCHELL',"WB OBEAR")
  areaSubset <- c('INSIDE', 'ABOVE', 'BELOW', 'TRIB','ABOVE ABOVE','aboveabove','BELOW BELOW' )

  subsetDMdataCohortMin <- 2002 # >=
  subsetDMdataCohortMax <- 2012 # <=

  # DON'T do this (by Year) because it messes up firstObsRows etc'
  #subsetDMdataYearMin <- 2002 # >=
  #subsetDMdataYearMax <- 2006 # <=

  subsetDMdataAgeInSamples <- 15 # <

  maxAgeInSamplesFirstCapt <- subsetDMdataAgeInSamples - 1 #4
}

# columns to include in  pheno2LongList
pheno2LongList <- c('tagNumber','sampleNumAdj','sampleNumAdjConsec',
                    'drainage','river','cohort','species',
                    'season','age','length','weight',
                    'section','riverM',
                    'sex','mature01','everMat01',
                    'date2','medianDate2','fromTo',
                    'everEmWest','everEmMain','area','riverM', 'riverN'
)
# columns to include in  pheno2LongList - need to update 'names' on next line
dMDataList <-  c('tagNumber','sampleNumAdj','sampleNumAdjConsec','length','weight','mature01','section',
                 'enc', 'river', 'speciesConsec',
                 'seasonConsec','drainageConsec','ageConsec','ageInSamplesConsec',
                 'yearConsec',
                 'gtFirstOcc','firstConsec','lastConsec','cohortConsec',
                 'date2','medianDate2','everEmConsec01','fromTo',
                 'everMat01Consec','area','riverM', 'riverN','sex'
)

dMDataNames <- c('tagNumber','sampleNumOrig','sampleNum','length','weight','mature01','section','enc', 'river',
                 'species','season','drainage','age','ageInSamples',
                 'year',
                 'gtFirstOcc','first','last','cohort',
                 'date2','medianDate2','everEmWest','fromTo',
                 'everMat01','area','riverM', 'riverN','sex'
)

#############################################################################

setwd("/var/lib/data/share")
#setwd('C:/PITTAGMAIN/CMR Analyses/Hierach_Bugs/allSpp')

pheno <- as.data.frame(read.csv(file="fourRiversAllSpp.txt", header=TRUE))
# file created using access query fourRiversAllSpp
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp\fourRiverAllSpp.accdb
# export query as text file without formatting, include names on first row
# then copy to felek in '/home/ben/allSpp/' and 'var/lib/data/share'

# NOW IN ECOPIFFLE, so commented out:
#source('./addLaggedForward.r')     #includes addEnvironmentalData2 (w/ skew)


#
pheno$length[ pheno$length == "-9999" ] <- NA  #get rid of -9999s
pheno$weight[ pheno$weight == "-9999" ] <- NA  #get rid of -9999s
pheno$riverM[ pheno$riverM == "-9999" ] <- NA  #get rid of -9999s
pheno$mature[ pheno$mature == "-9999" ] <- NA  #get rid of -9999s
pheno$sex[ pheno$sex == "-9999" | pheno$sex == "9999" ] <- NA  #get rid of Us (unknowns)
pheno$sex <- factor(pheno$sex)
#      pheno$age <- factor(pheno$age)
pheno$cohort <- factor(pheno$cohort)
#      pheno$agesInt <- as.integer(pheno$ages)
pheno <- subset(pheno,pheno$season != "Summer")     #get rid of the two summer obs
pheno$season[pheno$season == "Fall_Winter"] <- "Fall"
pheno$season <- factor(pheno$season,levels=c("PreSmolt","PostSmolt","Fall","PreWinter"), ordered=T)
#MAKE SURE there are only 4 seasons - otherwise numbering gets screwed up  when merging on season



pheno$date2 <- as.POSIXct(strptime(pheno$date, format = "%m/%d/%Y"))
pheno$medianDate2 <- as.POSIXct(strptime(pheno$medianDate, format = "%m/%d/%Y"))

#pheno$medianSampleDate2 <- as.POSIXct(strptime(pheno$medianSampleDate, format = "%m/%d/%Y"))
pheno$year <- as.factor(format(pheno$date2, format="%Y"))
pheno$julian <- as.numeric(format(pheno$date2,"%j"))
pheno$age <- as.numeric(as.character(pheno$year)) - as.numeric(as.character(pheno$cohort))
pheno$daysOld <- pheno$julian + as.numeric(as.character(pheno$age)) * 365

pheno$mature01 <- ifelse(pheno$mature == 'P',1,0)
matureList <- subset(pheno[,c('tagNumber','mature')],pheno$mature01 == 1 )
pheno$everMat01 <- ifelse(pheno$tagNumber %in% matureList$tagNumber,1,0)
#pheno <- pheno[ order(pheno$tagNumber,pheno$sampleNum),]

#pheno$agesLumpOld <- as.factor(ifelse(pheno$age %in% c(0,1), as.character(pheno$ages), as.character("other")))

#get rid of samples without many observations
#pheno <- subset(pheno, pheno$ages %in% c("0-3", "0-4","1-1","1-2","1-3","1-4","2-1","2-2","2-3","2-4","3-1","3-2","3-3"))
#xyplot(weight~length | ages, data=pheno)

pheno$riverN<-as.numeric(pheno$river)  #Jimmy 1, Mitchell 2, WB 3

#need to use 'date2' - it's in POSIX format
#pheno <- addLagged(data = pheno, individual = "tagNumber", time = "date2",
#           lag = c("riverN","riverM","length","weight","date2", "sampleNum"))


# SEEMS LIKE THIS NEEDS TO BE CHANGED BECASUE IT'S BASED ON CONSECTUTIVE CAPTURES IN PHENO, BUT SHOULD BE BASED
# ON CONSEC SAmPLeNUMNBERS THAT GET DEFINED IN PHENOcONSEC
# doint this down below now (9/1/2011)
#pheno <- addLaggedForward(data = pheno, individual = "tagNumber", time = "date2",
#           lag = c("riverN","riverM","length","weight","date2"))#, "sampleNum"))


#pheno$fromTo <- ifelse(!is.na(pheno$lagRiverN), paste(pheno$lagRiverN,pheno$riverN, sep=""),NA)
#pheno$fromTo <- as.factor(pheno$fromTo)

#pheno$move <- ifelse(pheno$riverN == pheno$lagRiverN,0,1)

#pheno$interval <- difftime(pheno$date2,pheno$lagDate2)
#pheno$moveDistance <- pheno$riverM - pheno$lagRiverM
#pheno$moveDistanceRate <- (pheno$riverM - pheno$lagRiverM) / as.numeric(pheno$interval)
#pheno$everEmPerm <- ifelse(substr(as.character(fromTo),2,2) == '0',1,0)

#pheno$grLength <- (pheno$length - pheno$lagLength) / as.numeric(pheno$interval)
#pheno$grWeight <- (log(pheno$weight) - log(pheno$lagWeight) ) /  as.numeric(pheno$interval)
#b <- 0.31
#pheno$grWeightOst <- ((pheno$weight ^ b) - (pheno$lagWeight ^ b)) / (b * as.numeric(pheno$interval))

#pheno$grWeightOst <- (log(pheno$weight) - log(pheno$lagWeight) ) /  as.numeric(pheno$interval)

#get some  negative growth rates in length because get recaptures within a sample, but in a different
# river. for now, will set length and weight growth rates for these fish to 0 to avoid ugly graphs - this may introduce some bias, though
#pheno$grLength[pheno$grLength < 0 ] <- 0
#pheno$grWeight[pheno$grLength < 0 ] <- 0
#pheno$grWeightOst[pheno$grLength < 0 ] <- 0

# Line up sample numbers across drainages
# if need to change sample numbers, run "C:/PITTAGMAIN/CMR Analyses/adjustedDates/sampleDates.r")
adjDates <- read.csv('./adjSampNums.csv', header=T)
pheno <- merge(pheno,adjDates[,c("drainage","sampleNum","sampleNumAdj")], by = c("drainage", "sampleNum"), all.x = T)

#rename wB sample 41.8 to 42. 42 was a sampleComple==NO, but there are equivalent samples in the other rivers
pheno$sampleNumAdj[ pheno$sampleNumAdj == 41.8 ] <- 42
# Delete any sampleNumberAdj that gets an NA (e.g. Shorey fall2 samples)
pheno <- pheno[! is.na(pheno$sampleNumAdj),]
#rename all sampleNumAdj == 34 to 'postSmolt' -Catamaran is named 'presmolt
#which adds an extra row to sampleToSeason (below) and then matrix dims
#don't match up - would be best to change this in the input file...
dim(unique(pheno[,c('sampleNumAdj','season')]))
pheno$season[ pheno$sampleNumAdj == 34 ] <- 'PostSmolt'
dim(unique(pheno[,c('sampleNumAdj','season')]))

# get rid of duplicate tags on the same occasion
pheno <- pheno[ !duplicated(pheno[,c('tagNumber','sampleNum')]), ]
#pheno[pheno$tagNumber=='1BF0DDFFB6',]
#pheno[pheno$tagNumber=='1BF20EA597',]

str(pheno)
##########################
# End of initial data prep
##########################

#############################################################################
# subset species
pheno2 <- subset(pheno, species %in% speciesSubset)
##########################
#subset river
pheno2 <- subset(pheno2, river %in% riverSubset)
##########################
#subset area
pheno2 <- subset(pheno2, area %in% areaSubset)
#############################################################################

#pheno2 <- subset(pheno2, sampComplete == 'YES')
pheno2$age <- factor(pheno2$age)
#pheno2$ages <- factor(pheno2$ages)
pheno2$cohort <- factor(pheno2$cohort)
pheno2$tagNumber <- factor(pheno2$tagNumber)
pheno2$river <- factor(pheno2$river)

#############################################################################
# Generate data to set up long data format
#############################################################################
firstYear <- min(as.numeric(as.character(pheno2$year)))
firstSeason <- as.numeric( pheno2$season[min(pheno2$sampleNumAdj)])

lastYear <- max(as.numeric(as.character(pheno2$year)))
#lastSeason is wrong, but doesn't get used below
lastSeason <- as.numeric( pheno2$season[max(pheno2$sampleNumAdj)])

# set up template for all possible samples
yearSeasonList <- as.data.frame(matrix(NA,(lastYear-firstYear+1)*4,2))
i=0
for (y in 1:(lastYear-firstYear+1)){
  for (s in 1:4){
    i=i+1
    yearSeasonList[i,1] <- y + firstYear -1
    yearSeasonList[i,2] <- s
  }
}
names(yearSeasonList) <- c('year','season')

# get list of actual samples from data
sampleYearSeasonList <- unique(pheno2[,c('sampleNumAdj','season','year')])
sampleYearSeasonList$season <- as.numeric(sampleYearSeasonList$season)
# sample 47 has some fish in 2003 and 2004 (caught in Jan)
#need to make sure sample 51 is year 2004
sampleYearSeasonList$year[sampleYearSeasonList$sampleNumAdj == 51] <- 2004
sampleYearSeasonList$year[sampleYearSeasonList$sampleNumAdj == 47] <- 2003

#get rid of repeats for sample 47 and 51
sampleYearSeasonList <- unique(sampleYearSeasonList)

#don't do this, puts NA in for sampleNumAdj and keeps year wrong
#delete sample 47 and 2004 in this list - otherwise double up sample 47 in 2003 and 2004
#sampleYearSeasonList <- subset(sampleYearSeasonList,
#                              (sampleNumAdj==47 & year == 2004)==FALSE)


# add sampleNumAdj to yearSeasonList (template)
yearSeasonList2 <- merge(
  x = yearSeasonList, y = sampleYearSeasonList,
  by = c('year','season'), all.x = TRUE
)

# fill in NAs for sampleNumAdj. Happens when species=ATS for winter samples
#just hard-coding this for now. sample #s will always line up w year/season
yearSeasonList2$sampleNumAdj[#yearSeasonList2$drainage == 'WEST' &  # will need to add this back in if run for multiple drainages at once
  yearSeasonList2$year == 2005 &
    yearSeasonList2$season == 4 ] <- 55
yearSeasonList2$sampleNumAdj[#yearSeasonList2$drainage == 'WEST' &
  yearSeasonList2$year == 2002 &
    yearSeasonList2$season == 4 ] <- 42  # don't seem to need first and last (2005,2007), but do need this one (2002)'
yearSeasonList2$sampleNumAdj[#yearSeasonList2$drainage == 'WEST' &
  yearSeasonList2$year == 2007 &
    yearSeasonList2$season == 4 ] <- 63

#get rid of leading unsampled samples
#yearSeasonList2 <- yearSeasonList2[
#                    min(which(!is.na(yearSeasonList2$sampleNumAdj))):
#                    nrow(yearSeasonList2),]
#instead of previous 3 lines, get rid of ANY unsampled samples - if left in, trailing NA sampleNumAdj screw up evalRows [1/12/2012]
yearSeasonList2 <- yearSeasonList2[ !is.na( yearSeasonList2$sampleNumAdj ), ]

yearSeasonList2$sampleNumAdjConsec <- 1:nrow(yearSeasonList2)

##################################################
# add preceding sampleNums so old fish caught early can be augmented back before the study started
##################################################
firstYear <- (yearSeasonList2$year)[1]
firstSeason <-(yearSeasonList2$season)[1]
numOccAug <- subsetDMdataAgeInSamples - 1

augBack <- as.data.frame(matrix(NA,numOccAug,ncol(yearSeasonList2)))
names(augBack) <- names(yearSeasonList2)
augBack$sampleNumAdjConsec <- seq(-(numOccAug-1),0) # -1 to account for sampleNumConsec==0

augBack$season[nrow(augBack)] <- firstSeason - 1
if (augBack$season[nrow(augBack)] == 0) augBack$season[nrow(augBack)] <- 4

augBack$year[nrow(augBack)] <- firstYear - 0
if (augBack$season[nrow(augBack)] == 4) augBack$year[nrow(augBack)] <- firstYear - 1

for ( i in (numOccAug - 1) : 1 ){
  augBack$season[i] <- augBack$season[i+1] - 1
  if (augBack$season[i] == 0) augBack$season[i] <- 4

  augBack$year[i] <- augBack$year[i+1] - 0
  if (augBack$season[i] == 4) augBack$year[i] <- augBack$year[i+1] - 1

}

yearSeasonList3 <- rbind(augBack,yearSeasonList2)
##################################################
##################################################

# add sampleNumAdjConsec to pheno2
pheno2 <- merge(
  x = pheno2, y = yearSeasonList2[,c('sampleNumAdj','sampleNumAdjConsec')],
  by = 'sampleNumAdj', all.x = TRUE
)

nOccTemplate <- nrow(yearSeasonList3)
nInd <- length(unique(pheno2$tagNumber))

##############################################################################
#make sampleNumConsec template for long format data
##############################################################################
sampNumTemplate <- rep(yearSeasonList3$sampleNumAdjConsec,nInd)
tagNumTemplate <- sort(rep( as.character(unique(pheno2$tagNumber)), nOccTemplate ))

template <- as.data.frame( cbind( tagNumTemplate, sampNumTemplate) )
names(template) <- c('tagNumber', 'sampleNumAdjConsec')


##############################################################################
# merge pheno2 into template of tags and sampleNumAdjConsec
##############################################################################
#pheno2$sampleNumConsec) <- 'sampleNumAdjConsec'

pheno2Long <- merge(
  x = template, y = pheno2[ , pheno2LongList ],
  by = c('tagNumber','sampleNumAdjConsec'), all.x = TRUE
)

pheno2Long$enc <- ifelse(is.na(pheno2Long$sampleNumAdj),0,1)

pheno2Long$tagNumber <- as.character(pheno2Long$tagNumber)
pheno2Long$sampleNumAdjConsec <- as.numeric(as.character(pheno2Long$sampleNumAdjConsec))
pheno2Long$age <- as.numeric(as.character(pheno2Long$age))
##############################################################################
# season index for each sample
##############################################################################
sampleToSeason <- unique(yearSeasonList3[,c('sampleNumAdjConsec','season')])
sampleToSeason$seasonN <- as.integer(sampleToSeason$season)
#sort by sampleNumAdj - otherwise yearsMatrix is out of order
sampleToSeason <- sampleToSeason[ order(sampleToSeason$sampleNumAdjConsec), ]

##############################################################################
# year index for each sample
##############################################################################
sampleToYear <- unique(yearSeasonList3[,c('sampleNumAdjConsec','year')])
sampleToYear$yearN <- as.numeric(as.character(sampleToYear$year))
#sort by sampleNumAdj - otherwise yearsMatrix is out of order
sampleToYear <- sampleToYear[ order(sampleToYear$sampleNumAdjConsec), ]

sampleTo <- as.data.frame(
  cbind( sampleToSeason$sampleNumAdjConsec,
         sampleToSeason$seasonN,
         sampleToYear$yearN ) )
names(sampleTo) <- c( 'sampleNumAdjConsec','seasonConsec','yearConsec' )

pheno2Long <- merge(
  x = pheno2Long, y = sampleTo,
  by = c('sampleNumAdjConsec'), all.x = TRUE, sort=FALSE
)

##############################################################################
#get individual-specifc lists and prepare to merge into pheno2Long
##############################################################################
##############################################################################
# variables that don't change within fish
##############################################################################
first <- aggregate((pheno2[,c('sampleNumAdjConsec') ] ),
                   by=list(( pheno2$tagNumber) ),FUN=min, na.rm =TRUE )
names(first) <- c('tagNumber','firstConsec')

last <- aggregate((pheno2[,c('sampleNumAdjConsec') ] ),
                  by=list(( pheno2$tagNumber) ),FUN=max, na.rm =TRUE )
names(last) <- c('tagNumber2','lastConsec')

firstLast <- cbind(first,last)

pheno2Long <- merge(
  x = pheno2Long, y = firstLast[, c('tagNumber','firstConsec',
                                    'lastConsec')],
  by = 'tagNumber', all.x = TRUE
)

# first/last and cohort/drainage have different #'s of individauls. need ot merge separately
# make matix of ages for each fish
cohortList <- unique(pheno2[,c('tagNumber','cohort')])
cohortList$tagNumber <- as.character(cohortList$tagNumber)
cohortList$cohort <- as.numeric(as.character(cohortList$cohort))
names(cohortList) <- c('tagNumber','cohortConsec')

# make matix of drainage for each fish
drainageList <- unique(pheno2[,c('tagNumber','drainage')])
drainageList$tagNumber <- as.character(drainageList$tagNumber)
drainageList$drainage <- as.character(drainageList$drainage)
names(drainageList) <- c('tagNumber','drainageConsec')

# make matix of everMat01 for each fish
everMat01List <- unique(pheno2[,c('tagNumber','everMat01')])
everMat01List$tagNumber <- as.character(everMat01List$tagNumber)
everMat01List$everMat01 <- as.character(everMat01List$everMat01)
names(everMat01List) <- c('tagNumber','everMat01Consec')

# make matix of species for each fish
speciesList <- unique(pheno2[,c('tagNumber','species')])
speciesList$tagNumber <- as.character(speciesList$tagNumber)
speciesList$species <- as.character(speciesList$species)
names(speciesList) <- c('tagNumber','speciesConsec')


cohortDrainageEverMat01 <- cbind(cohortList,drainageList,everMat01List,speciesList)

pheno2Long <- merge(
  x = pheno2Long, y = cohortDrainageEverMat01[, c('tagNumber',
                                                  'cohortConsec','drainageConsec','everMat01Consec','speciesConsec')],
  by = 'tagNumber', all.x = TRUE
)

#create ageConsec column for each fish
pheno2Long$ageConsec <- pheno2Long$yearConsec - pheno2Long$cohortConsec
# + 2 so season 2 is ageInSamples
pheno2Long$ageInSamplesConsec <- 4*(pheno2Long$ageConsec-1) + 2 + pheno2Long$seasonConsec

# delete fish that were caught for the first time too old
#hist(pheno2Long[pheno2Long$sampleNumAdjConsec == pheno2Long$firstConsec, c('ageInSamplesConsec')] )
firstAIS <- (pheno2Long[ pheno2Long$sampleNumAdjConsec == pheno2Long$firstConsec , c('tagNumber','ageInSamplesConsec') ] )
firstNotTooOld <- firstAIS[firstAIS$ageInSamples <= maxAgeInSamplesFirstCapt, 'tagNumber']
pheno2Long <- subset(pheno2Long, pheno2Long$tagNumber %in% firstNotTooOld)

#need to do everEm Separately - tag nums don;t line up. maybe some values missing
# make matix of everEmMain for each fish
emList <- unique(pheno2[,c('tagNumber','everEmWest')])   #ok to use everEmWest across rivers. identical to everEmMain, except for WB
emList$tagNumber <- as.character(emList$tagNumber)
emList$everEmWest01 <- ifelse( emList$everEmWest == 'YES', 1, 0 )
#emList$everEmMain <- as.numeric(as.character(emList$everEmMain))
names(emList) <- c('tagNumber','everEmConsec','everEmConsec01')

pheno2Long <- merge(
  x = pheno2Long, y = emList[, c('tagNumber','everEmConsec01')],
  by = 'tagNumber', all.x = TRUE
)

pheno2Long <- pheno2Long[ order(pheno2Long$tagNumber,
                                pheno2Long$sampleNumAdjConsec), ]

#pheno2Long[pheno2Long$tagNumber=='1BF20EA597',]


##############################################################################
# prepare to create availability column
##############################################################################

#returns Inf if never gets to subsetDMdataAgeInSamples w/in ind
oldestSampConsec <- aggregate((pheno2Long[,c('ageInSamplesConsec') ] ),
                              by=list( pheno2Long$tagNumber ),
                              FUN = function(x){ pheno2Long$sampleNumAdjConsec[
                                min(which(x==(subsetDMdataAgeInSamples-1))) - 0
                                ] } )
# min(which( pheno2Long$ageInSamplesConsec[pheno2Long$tagNumber=='1BF0DDFFB6']==(subsetDMdataAgeInSamples-1)))
names(oldestSampConsec) <- c('tagNumber','lastPossibleConsec')

# set lastPossible to the last sample if the fish never got older than subsetDMdataAgeInSamples
nOcc <- max(yearSeasonList3$sampleNumAdjConsec)
oldestSampConsec$lastPossibleConsec <- ifelse(
  is.finite(oldestSampConsec$lastPossibleConsec)==TRUE,
  oldestSampConsec$lastPossibleConsec,
  nOcc)

pheno2Long <- merge(
  x = pheno2Long, y = oldestSampConsec,
  by = 'tagNumber', all.x = TRUE
)

# [used to] truncate rows for permanent emigrants
pheno2Long$lastPossible01 <- ifelse( pheno2Long$everEmConsec01 == 0,
                                     pheno2Long$lastPossibleConsec,
                                     pheno2Long$lastPossibleConsec)   # with this line on, not truncating for everEm yes fish. Using emigration model to do this now
#pheno2Long$lastConsec )         # this line truncates

# fish are not available when they are too old                        ##or if they have emigrated
# fish are not available before first capture
pheno2Long$available <- ifelse(
  ( pheno2Long$sampleNumAdjConsec <= pheno2Long$lastPossible01 ) &
    #     ( pheno2Long$ageInSamplesConsec >= 1 | pheno2Long$sampleNumAdjConsec >= pheno2Long$firstConsec ),
    #the previous line augments backwards to ageInSamples == 1 or smaller if the fish was captured in the summer 0+
    # the following lines starts at first capture
    ( pheno2Long$sampleNumAdjConsec >= pheno2Long$firstConsec ),
  1, 0
)

# delete fish that were seen for the first time after subsetDMdataAgeInSamples
tooOld <- unique(pheno2Long$tagNumber[ pheno2Long$lastPossibleConsec < pheno2Long$firstConsec ])
pheno2Long <- pheno2Long[ !(pheno2Long$tagNumber %in% tooOld), ]

# need to make fish unavailable '' and after the end    # not doing this any more-> 'before the beginning of a study'
minSample <- melt( tapply(pheno2Long$sampleNumAdjConsec[
  pheno2Long$ageConsec >= 0 & pheno2Long$seasonConsec > 2],
  INDEX=pheno2Long$drainageConsec[
    pheno2Long$ageConsec >= 0 & pheno2Long$seasonConsec > 2],
  FUN=min) )
names(minSample) <- c('drainageConsec', 'minSampleNumConsec')

#had '> 2' instead of '>0' in earlier version. this assumed the last sample was greater than season 2 which deosn't work when the last sample is in season 1 or 2.
#end up with no data for fish caught in the last sample
maxSample <-  melt( tapply(pheno2Long$sampleNumAdjConsec[
  pheno2Long$ageConsec >= 0 & pheno2Long$seasonConsec > 0],
  INDEX=pheno2Long$drainageConsec[
    pheno2Long$ageConsec >= 0 & pheno2Long$seasonConsec > 0],
  FUN=max) )
names(maxSample) <- c('drainageConsec', 'maxSampleNumConsec')
minMaxSampleNum <- cbind(minSample,maxSample)
minMaxSampleNum <- minMaxSampleNum[,-1]

pheno2Long <- merge(
  x = pheno2Long, y = minMaxSampleNum,
  by = 'drainageConsec', all.x = TRUE
)

pheno2Long$available <- ifelse(
  ( pheno2Long$sampleNumAdjConsec < pheno2Long$minSampleNumConsec  ) |
    ( pheno2Long$sampleNumAdjConsec > pheno2Long$maxSampleNumConsec  ),
  0, pheno2Long$available
)



# indicator for first sample
pheno2Long$gtFirstOcc <- ifelse(
  pheno2Long$sampleNumAdjConsec >
    pheno2Long$firstConsec,
  1, 0
)


# pheno2Long[pheno2Long$tagNumber=='1BF20EA597',]
pheno2Long <- pheno2Long[ order(pheno2Long$tagNumber, pheno2Long$sampleNumAdjConsec),]




##############################################################################
# subset for availability ==1. Final input dataset for analysis
##############################################################################

dMData <- subset( pheno2Long[ , dMDataList ] ,
                  pheno2Long$available == 1 )

#for now get rid of fish with ageInSamplesConsec == -1
minusOnes <- subset(dMData,dMData$ageInSamplesConsec == -1)
dMData <- dMData[ which((dMData$tagNumber %in% minusOnes$tagNumber) == FALSE ), ]

# dMDataList names defined up top
names(dMData) <- dMDataNames

# Subset dMData
#dMData <- subset(dMData, year >= subsetDMdataYearMin & year <= subsetDMdataYearMax)
dMData <- subset(dMData, cohort >= subsetDMdataCohortMin & cohort <= subsetDMdataCohortMax)
dMData <- subset(dMData, ageInSamples < subsetDMdataAgeInSamples)


dMData$drainage <- as.factor(dMData$drainage)
dMData$season <- as.factor(dMData$season)
dim(dMData)
#dMData <- subset(dMData,as.numeric(as.character(dMData$age)) < maxAgePlusOne )      #get rid of old fish. not best place for this
dim(dMData)
dMData$age <- factor(dMData$age)


##############################################################################
#add in proportion of sections sampled
#ftable(pheno2$drainage,pheno2$year,pheno2$season,pheno2$section)
numDrainages <- length(unique(dMData$drainage))
yearSeasonList4 <- yearSeasonList3
drainages <- levels(unique(dMData$drainage))
yearSeasonList4$drainage <- drainages[1]

if (numDrainages > 1){
  for (i in 2:numDrainages){
    yearSeasonListHold <- yearSeasonList3
    yearSeasonListHold$drainage <- drainages[i]
    yearSeasonList4 <- rbind(yearSeasonList4,yearSeasonListHold)
  }
}

yearSeasonList4$propSampled <- 1

# fix incomplete samples to the proportion of sections sampled

# WEST - 2002, no river was sampled in winter 2002
yearSeasonList4$propSampled[ yearSeasonList4$drainage == 'WEST' &
                               yearSeasonList4$year == 2002 &
                               yearSeasonList4$season == 4 ] <- 0

# West Brook samples were not complete in the following years in winter
yearSeasonList4$propSampled[ yearSeasonList4$river == 'WEST BROOK' &
                               yearSeasonList4$year == 2003 &
                               yearSeasonList4$season == 4 ] <- 30/47

yearSeasonList4$propSampled[ yearSeasonList4$river == 'WEST BROOK' &
                               yearSeasonList4$year == 2004 &
                               yearSeasonList4$season == 4 ] <- 3/47

yearSeasonList4$propSampled[ yearSeasonList4$river == 'WEST BROOK' &
                               yearSeasonList4$year == 2005 &
                               yearSeasonList4$season == 4 ] <- 0

yearSeasonList4$propSampled[ yearSeasonList4$river == 'WEST BROOK' &
                               yearSeasonList4$year == 2007 &
                               yearSeasonList4$season == 4 ] <- 0
# SHorey
# seems like all samples were complete



# merge yearSeasonList2$propSampledATS into dMData
yearSeasonList4$sampleNum <- yearSeasonList4$sampleNumAdjConsec
dMData <- merge(
  x = dMData, y = yearSeasonList4[,c('drainage','sampleNum','propSampled')],
  by = c('drainage','sampleNum'), all.x = TRUE
)
##############################################################################


##############################################################################
# add environmental data to each row of dMData
##############################################################################
envData <- as.data.frame(read.csv(file="./envData.txt", header=TRUE))
# file created using access file query  fourRiversAllSpp:envData query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
#export query as text file without formatting
# then copy to felek

names(envData) <- c('drainage','river','date','temperature','discharge','rain')
envData$date2 <- as.POSIXct(strptime(envData$date, format = "%m/%d/%Y"))
envData$year <-  format(envData$date2,'%Y')
envData$julian <- as.numeric(format(envData$date2,"%j"))

#need to add medianDate and lagMedianDate to row for which fish were not caught
# fill in median dates for samples with no data (e.g winters)

minDateListForMissing <- aggregate(dMData$date2[dMData$enc==1],
                                   by=list(dMData$drainage[dMData$enc==1],dMData$sampleNum[dMData$enc==1]),
                                   FUN=min,na.rm=TRUE)
names(minDateListForMissing) <- c('drainage','sampleNum','minDateSampleNum')
minDateListForMissing$minDateSampleNum <-
  as.POSIXct(minDateListForMissing$minDateSampleNum, origin="1970-01-01" )  #default origin


medianDateListForMissing <- aggregate(dMData$date2[dMData$enc==1],
                                      by=list(dMData$drainage[dMData$enc==1],dMData$sampleNum[dMData$enc==1]),
                                      FUN=median,na.rm=TRUE)
names(medianDateListForMissing) <- c('drainage','sampleNum','medianDateSampleNum')
medianDateListForMissing$medianDateSampleNum <-
  as.POSIXct(medianDateListForMissing$medianDate, origin="1970-01-01" )  #default origin

maxDateListForMissing <- aggregate(dMData$date2[dMData$enc==1],
                                   by=list(dMData$drainage[dMData$enc==1],dMData$sampleNum[dMData$enc==1]),
                                   FUN=max,na.rm=TRUE)
names(maxDateListForMissing) <- c('drainage','sampleNum','maxDateSampleNum')
maxDateListForMissing$maxDateSampleNum <-
  as.POSIXct(maxDateListForMissing$maxDateSampleNum, origin="1970-01-01" )  #default origin
maxDateListForMissing$julianMax <- as.numeric(format(maxDateListForMissing$maxDateSampleNum,"%j"))


#for dat and graphs below
medianDateListForMissing$year <- format(medianDateListForMissing$medianDateSampleNum, '%Y')
medianDateListForMissing$julianMedian <- as.numeric(format(medianDateListForMissing$medianDateSampleNum, '%j'))
julianMin<-as.numeric(format(minDateListForMissing$minDateSampleNum, '%j'))
julianMax<-as.numeric(format(maxDateListForMissing$maxDateSampleNum, '%j'))
dateMin<-minDateListForMissing$minDateSampleNum
dateMax<-maxDateListForMissing$maxDateSampleNum

medianDateListForMissing <- cbind(medianDateListForMissing,
                                  julianMin,
                                  julianMax,
                                  dateMin,
                                  dateMax)

# saving for antennaDuda.r, so we can make graphs of antenna data
save(medianDateListForMissing, file= './medianDateListForMissing.RData')

#calulate T and discharge means between dateMin and dataMax for each sample
# I'm sure these's a better way to do this, but this works...
meansDuringSample <- as.data.frame(matrix(NA,nrow(medianDateListForMissing),6))
names(meansDuringSample) <- c('drainage','sampleNum',
                              'tempDuringOcc','dischDuringOcc',
                              'tempDuringOccSampDays','dischDuringOccSampDays')

sampleDateList <- unique(dMData[,c('drainage','sampleNum','date2')])

sampleDateList <- sampleDateList[ order(sampleDateList$drainage,
                                        sampleDateList$sampleNum,
                                        sampleDateList$date2),]
sampleDateList <- sampleDateList[!is.na(sampleDateList$date2),]

for (i in 1:nrow(medianDateListForMissing)){
  meansDuringSample[i,1] <- as.character(medianDateListForMissing$drainage[i])
  meansDuringSample[i,2] <- medianDateListForMissing$sampleNum[i]
  # don't really need cols 3 and 4
  meansDuringSample[i,3] <- mean(envData$temperature[
    envData$drainage ==
      as.character(medianDateListForMissing$drainage[i]) &
      envData$date2 >=
      medianDateListForMissing$dateMin[i] &
      envData$date2 <=
      medianDateListForMissing$dateMax[i] ],
    na.rm=TRUE)
  meansDuringSample[i,4] <- mean(envData$discharge[
    envData$drainage ==
      as.character(medianDateListForMissing$drainage[i]) &
      envData$date2 >=
      medianDateListForMissing$dateMin[i] &
      envData$date2 <=
      medianDateListForMissing$dateMax[i] ],
    na.rm=TRUE)
  meansDuringSample[i,5] <- mean(envData$temperature[
    envData$date2 %in%

      sampleDateList$date2[sampleDateList$sampleNum==
                             medianDateListForMissing$sampleNum[i]] &

      envData$drainage ==
      as.character(medianDateListForMissing$drainage[i]) ]
    , na.rm=TRUE)
  meansDuringSample[i,6] <- mean(envData$discharge[
    envData$date2 %in%

      sampleDateList$date2[sampleDateList$sampleNum==
                             medianDateListForMissing$sampleNum[i]] &

      envData$drainage ==
      as.character(medianDateListForMissing$drainage[i]) ]
    , na.rm=TRUE)

}
# use tempDuringOccSampDays for analyses - is the mean of actual days samples, not over range like tempDuringOcc

#need to fill in NaN entires in meansDuringSample. 1st get means, then merge in when entry is NaN
sampleToSeason2 <- sampleToSeason
names(sampleToSeason2)<- c('sampleNum','season','season2')
meansDuringSample <- merge(
  x = meansDuringSample, y = sampleToSeason2[,c('sampleNum','season')],
  by = c('sampleNum'), all.x = TRUE
)

# add rows for samples not Attempted
missed <- (matrix(NA,nrow(meansDuringSample),1))
ii <- 0
for ( i in 1:( nrow(meansDuringSample) - 1 ) ) {
  if( meansDuringSample$sampleNum[ i+1 ] != meansDuringSample$sampleNum[ i ] + 1 ) {
    ii <- ii + 1
    missed[ ii ] <- meansDuringSample$sampleNum[ i ] + 1
    print(c(i,ii,missed[ii]))
  }
  #print(c(i,meansDuringSample$sampleNum[ i ],meansDuringSample$sampleNum[ i+1 ] , meansDuringSample$sampleNum[ i ] + 1))
}
missedSamples <- missed[!is.na(missed)]

temp <- as.data.frame(matrix(NA,length(missedSamples),ncol(meansDuringSample))) #3

temp[,1] <- missedSamples #c(23,35,43)           #sampleNum
temp[,2] <- rep( unique(dMData$drainage), length(missedSamples) ) #c('WEST','WEST','WEST')   #drainage
temp[,7] <- rep( 4, length(missedSamples) ) #c(4,4,4)             #season - assume all are in winter
names(temp) <- names(meansDuringSample)

meansDuringSample <- rbind(meansDuringSample,temp)

duplicated(meansDuringSample$sampleNum)

meansBySeason <-  aggregate(meansDuringSample,
                            by=list(meansDuringSample$drainage,meansDuringSample$season),
                            FUN=mean,na.rm=TRUE)
meansBySeason <- meansBySeason[,c('Group.1', 'Group.2', 'tempDuringOcc', 'dischDuringOcc', 'tempDuringOccSampDays', 'dischDuringOccSampDays')]
names(meansBySeason) <- c('drainage', 'season', 'tempDuringOcc', 'dischDuringOcc', 'tempDuringOccSampDays', 'dischDuringOccSampDays')

# replace NaNs in meansDuringSample. these come from incomplete data in envData
#should thoroughly check indexing and column names here
for ( i in 1:nrow(meansDuringSample) ) {
  for ( j in 3:6 ) {
    if ( is.na(meansDuringSample[ i,j ]) ) {
      holdI <- which(meansDuringSample$drainage[i] == meansBySeason$drainage  &
                       meansDuringSample$season[i] == meansBySeason$season )
      meansDuringSample[ i,j ] <- meansBySeason[ holdI,j ]
    }
  }
}

# for ATS, 33 is NaN becasue minDate=maxDate
#add T and flow to each row of data as covariates for p(capt)
# merging on drainage for now because data are incomplete for river
dMData <- merge(
  x = dMData, y = meansDuringSample[,c('drainage','sampleNum',
                                       'tempDuringOccSampDays', 'dischDuringOccSampDays')],
  by = c('drainage','sampleNum'), all.x = TRUE
)

#merge in actual temp and disch for each fish
dMData <- merge(
  x = dMData, y = envData[,c('river','date2',
                             'temperature', 'discharge')],
  by = c('river','date2'), all.x = TRUE
)



#keep actual env data if capture, otherwise keep mean during sample.
# with movement model, could estimate location and then estimate day of capture to get est env data
dMData$temperatureForP <- ifelse(is.na(dMData$date2)==TRUE |
                                   is.na(dMData$temperature)==TRUE,  # this adds in data when it's missing for the date of capture (SHOREY)
                                 dMData$tempDuringOccSampDays,
                                 dMData$temperature )
dMData$dischargeForP <- ifelse(is.na(dMData$date2)==TRUE |
                                 is.na(dMData$discharge)==TRUE,
                               dMData$dischDuringOccSampDays,
                               dMData$discharge )




#ggplot(envData,aes(date2,temperature)) + geom_point() +
#  facet_wrap(~drainage)
#ggplot(subset(envData,drainage=='WEST'),aes(date2,temperature)) + geom_point() +
#  facet_wrap(~river)
#ggplot(subset(envData,river=='WEST BROOK'),aes(julian,(discharge+0.0)))  +
#  geom_line() + ylim(c(0,15)) + theme_bw() +
#  facet_wrap(~year) +
#  geom_vline(aes(xintercept=julianMin),medianDateListForMissing, colour='blue') +
#  geom_vline(aes(xintercept=julianMedian),medianDateListForMissing, linetype=2) +
#  geom_vline(aes(xintercept=julianMax),medianDateListForMissing,colour='red')
#
dMData <- merge(
  x = dMData, y = medianDateListForMissing[,c('drainage','sampleNum','medianDateSampleNum')],
  by = c('drainage','sampleNum'), all.x = TRUE
)


# set up template for all possible samples for antenna efficiency
firstYearAnt <- 1997; lastYearAnt <- 2015 # need ot update this for samnles > 2015
yearSeasonListAnt <- as.data.frame(matrix(NA,(lastYearAnt-firstYearAnt+1)*4,2))
i=0
for (y in 1:(lastYearAnt-firstYearAnt+1)){
  for (s in 1:4){
    i=i+1
    yearSeasonListAnt[i,1] <- y + firstYearAnt -1
    yearSeasonListAnt[i,2] <- s
  }
}
names(yearSeasonListAnt) <- c('year','season')
yearSeasonListAnt$propDaysOn <- 0.1    # set default to 0.1 so fish below have a chance of detection when antennas were not working
#numbers based on graph in antennaDuda.r. Fill in rest for other year/season combos (e.g. when run trout)
yearSeasonListAnt$propDaysOn[19:43] <- (c(10,80,100,95,100,65,100,80,90,100,100,
                                          80,100,100,100,90,100,50,100,100,100,
                                          100,95,85,90) ) / 100
yearSeasonListAnt$drainage <- 'WEST'
#yearSeasonListAnt$antEfficiency <- 1 #0.9 #make interval-specific
#yearSeasonListAnt$propDaysOnAntEff <- yearSeasonListAnt$propDaysOn *
#                                  yearSeasonListAnt$antEfficiency

dMData <- merge(
  x = dMData, y = yearSeasonListAnt[,c('year','season','drainage','propDaysOn')],
  by = c('drainage','year','season'), all.x = TRUE
)

dMData <- dMData[ order(dMData$tagNumber,dMData$sampleNum),]

#hard code missing winter samples if they are missing
dMData$medianDateSampleNum[ is.na( dMData$medianDateSampleNum ) &
                              dMData$drainage == 'WEST' &
                              dMData$year == 2002 &
                              dMData$season == 4 ] <- '2002-12-04'
dMData$medianDateSampleNum[ is.na( dMData$medianDateSampleNum ) &
                              dMData$drainage == 'WEST' &
                              dMData$year == 2005 &
                              dMData$season == 4 ] <- '2005-12-04'
dMData$medianDateSampleNum[ is.na( dMData$medianDateSampleNum ) &
                              dMData$drainage == 'WEST' &
                              dMData$year == 2007 &
                              dMData$season == 4 ] <- '2007-12-04'

dMData$dateForEnv <- ifelse(is.na(dMData$date2)==TRUE,
                            dMData$medianDateSampleNum,
                            dMData$date2 )
dMData$dateForEnv <- as.POSIXct(dMData$dateForEnv, origin="1970-01-01" )  #default origin

dMData <- dMData[ order(dMData$tagNumber,dMData$sampleNum),]

#table(dMData$sampleNum,dMData$river)
#yearSeasonList3
dMData <- addLaggedForward(data = dMData, individual = "tagNumber", time = "dateForEnv",
                           lag = c("dateForEnv"))

#names(dMData[['lagDateForEnv']]) <- 'lagDateForEnvForward'
dMData$tagNumberCH <- as.character(dMData$tagNumber)
############################################################
# add in smolt trap data for salmon
############################################################
if (speciesSubset == 'ATS') {

  #setwd("C://PITTAGMAIN//CMR Analyses//Hierach_Bugs//allSpp//")
  smoltTrap <- as.data.frame(read.csv(file="./smolttrapFourRiversAllSpp.txt", header=TRUE))
  # file created using access file query  fourRiversAts:smoltTrapfourRiversAts query
  # in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
  #export query as text file without formatting
  smoltTrap$date2 <- as.POSIXct(strptime(smoltTrap$date, format = "%m/%d/%Y"))
  smoltTrap$tagNumberCH <- as.character(smoltTrap$tagNumber)

  smoltTrap$year <- as.factor(format(smoltTrap$date2, format="%Y"))
  smoltTrap$age <- as.numeric(as.character(smoltTrap$year)) - as.numeric(as.character(smoltTrap$cohort))
  smoltTrap$ais <- ifelse(smoltTrap$age == 2, 7,ifelse(smoltTrap$age == 3, 11, NA))


  # get rid of duplicate tagNumbers
  # ddply(dMData[,c('tagNumber','ageInSamples')],.(tagNumber), summarize, max=max(dMData$ageInSamples))
  smoltTrap <- smoltTrap[!duplicated(smoltTrap[,c('tagNumber')]),]

  #hold2 <- ddply(dMData[,c('tagNumber','ageInSamples')],
  #               .(tagNumber),
  #               function(x) (x$tagNumber %in% smoltTrap$tagNumber &
  #                            x$ageInSamples %in% smoltTrap$ais) )

  #smoltTrap[smoltTrap$tagNumber=='1BF0DFFDC7',]

  dMData$smoltTrap <- NA
  smoltTrap$smoltTrap <- 1
  smoltTrap$ageInSamples <- smoltTrap$ais

  dMData <- merge(
    x = dMData, y = smoltTrap[,c('tagNumberCH','ageInSamples','smoltTrap')],
    by = c('tagNumberCH','ageInSamples'), all.x = TRUE
  )

}

#dMData[dMData$tagNumber=='1BF0DFFDC7',]
#dMData[dMData$tagNumber=='1BF0DFFF26',]

############################################################
# separates first occ from others for likelihood
############################################################

# use below if there is only one row of gtFirstOcc per fish
# not true when augmenting back to ageInSamples==1
#firstObsRows <- subset(1:nrow(dMData),dMData$gtFirstOcc==0)

lastAIS <- aggregate((dMData[,c('ageInSamples') ] ),
                     by=list(( dMData$tagNumber) ),FUN=max, na.rm =TRUE )
names(lastAIS) <- c('tagNumber','lastAIS')

dMData <- merge(
  x = dMData, y = lastAIS[, c('tagNumber','lastAIS')],
  by = 'tagNumber', all.x = TRUE
)

#if limit years so that early samples for an older cohort are cut off, can get fish with 0
#observations. This messes up the evalRows etc
#minSAmpleNum <- min(dMData$sampleNum)
#dMData2  <- dMData[ dMData$last > minSAmpleNum, ]

dMData <- dMData[ order(dMData$tagNumber,dMData$sampleNum),]  # MUST order before doing eval rows.

firstObsRows= which( dMData$first == dMData$sampleNum )
nFirstObsRows <- length( firstObsRows )
lastObsRows <- subset(1:nrow(dMData),is.na(dMData$lagDateForEnv))
nLastObsRows <- length(lastObsRows)

#firstObsJSRows= which( dMData$ageInSamples == 1 ) # make sure all fish have an ageInSample==1. 1996 cohort fish won't -not augmented back currently
firstObsJSRows <- which(dMData$sampleNum == dMData$first ) # make sure sampleNums cover all obs of first
nFirstObsJSRows <- length( firstObsJSRows )

evalRows <- subset( 1:nrow(dMData),
                    is.na(dMData$lagDateForEnv)==FALSE &
                      dMData$sampleNum >= dMData$first )
nEvalRows <- length(evalRows)

evalJSRows <- subset(1:nrow(dMData), dMData$ageInSamples != dMData$lastAIS )
nEvalJSRows <- length(evalJSRows)
############################################################

# variables for estimating pMat
summerObsRows <- subset(1:nrow(dMData),dMData$season==2)
nSummerObsRows <- length(summerObsRows)

nonSummerObsRows <- subset(1:nrow(dMData),dMData$season!=2)
nNonSummerObsRows <- length(nonSummerObsRows)

occasions <- unique(dMData[,c('ageInSamples','season')])
summerAIS <- occasions$ageInSamples[occasions$season == 2]
nonSummerAIS <- occasions$ageInSamples[occasions$season != 2]

# 1= matured in that maturity year (suummer-spring)
# 0 = observed immature in fall
# NA = not observed mature and not captured in fall

#age1
dMDataTemp <- subset(dMData, dMData$ageInSamples %in% 4:7)
age1Mat <- aggregate((dMDataTemp[,c('mature01') ] ),
                     by=list( dMDataTemp$tagNumber),FUN=sum, na.rm =TRUE )
age1Mat$mat <- ifelse( age1Mat$x > 0, 1,0 )
age1Mat$ageInSamples <- 4
names(age1Mat)[1:2] <- c('tagNumber','sum')

dMDataTemp2 <- subset(dMData, dMData$ageInSamples %in% 5)
dMDataTemp2$ageInSamples <- dMDataTemp2$ageInSamples - 1
age1Mat <- merge(
  x = age1Mat, y = dMDataTemp2[,c('tagNumber','ageInSamples','mature01')],
  by = c('tagNumber','ageInSamples'), all.x = TRUE
)

# age2
dMDataTemp <- subset(dMData, dMData$ageInSamples %in% 8:11)
age2Mat <- aggregate((dMDataTemp[,c('mature01') ] ),
                     by=list( dMDataTemp$tagNumber ),FUN=sum, na.rm =TRUE )
age2Mat$mat <- ifelse( age2Mat$x > 0, 1,0 )
age2Mat$ageInSamples <- 8
names(age2Mat)[1:2] <- c('tagNumber','sum')

dMDataTemp2 <- subset(dMData, dMData$ageInSamples %in% 9)
dMDataTemp2$ageInSamples <- dMDataTemp2$ageInSamples - 1
age2Mat <- merge(
  x = age2Mat, y = dMDataTemp2[,c('tagNumber','ageInSamples','mature01')],
  by = c('tagNumber','ageInSamples'), all.x = TRUE
)

#age3
dMDataTemp <- subset(dMData, dMData$ageInSamples %in% 12:14)
age3Mat <- aggregate((dMDataTemp[,c('mature01') ] ),
                     by=list( dMDataTemp$tagNumber ),FUN=sum, na.rm =TRUE )
age3Mat$mat <- ifelse( age3Mat$x > 0, 1,0 )
age3Mat$ageInSamples <- 12
names(age3Mat)[1:2] <- c('tagNumber','sum')

dMDataTemp2 <- subset(dMData, dMData$ageInSamples %in% 13)
dMDataTemp2$ageInSamples <- dMDataTemp2$ageInSamples - 1
age3Mat <- merge(
  x = age3Mat, y = dMDataTemp2[,c('tagNumber','ageInSamples','mature01')],
  by = c('tagNumber','ageInSamples'), all.x = TRUE
)



ageAllMat <- rbind(age1Mat,age2Mat,age3Mat)
ageAllMat$matNA <- ageAllMat$mat
# set 0s to NA to start
ageAllMat$matNA[ageAllMat$matNA == 0] <- NA
# set seen in fall, but not mature, to 0
ageAllMat$matNA[ is.na(ageAllMat$matNA) & ageAllMat$mature01 == 0] <- 0

dMData <- merge(
  x = dMData, y = ageAllMat[,c('tagNumber','ageInSamples','matNA')],
  by = c('tagNumber','ageInSamples'), all.x = TRUE
)

dMData$intervalLength <- difftime(dMData$lagDateForEnv,dMData$dateForEnv, unit="days")

#########################################################################
### AntennaDuda
#########################################################################
# make col of rows lastseen when fish were coded as permanent emigrants
# fish acutally emigrated after this (antenna or outside)

antennaDuda <- as.data.frame(read.csv(file="./antennaDuda.txt", header=TRUE))
# file created using access file query  fourRiversAts:antennaDuda query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
#export query as text file without formatting
# then copy to felek

antennaDuda$date2 <- as.POSIXct(strptime(antennaDuda$date, format = "%m/%d/%Y"))
antennaDuda$tagNumberCH <- as.character(antennaDuda$tagNumber)

# maxDateAntenna <- ddply( antennaDuda, .(tagNumber), summarize, maxAntennaDate2=max(antennaDuda$date2) )

# get last (max antenna date for each fish
maxDateAntenna <- unique( ddply( antennaDuda[,c('tagNumberCH','date2')], .(tagNumberCH), subset, date2==max(date2) ) )
names(maxDateAntenna) <- c('tagNumberCH','maxAntennaDate2')

# get last sampled date for each fish
lastDate <- dMData[ dMData$sampleNum == dMData$last, c('tagNumberCH','date2') ]
names(lastDate) <- c('tagNumberCH','lastDate2')

#antennaDuda[antennaDuda$tagNumber=='1748C762DD',]
#maxDateAntenna[maxDateAntenna$tagNumber=='1748C762DD',]


dMData <- merge(
  x = dMData, y = maxDateAntenna,
  by = c('tagNumberCH'), all.x = TRUE
)
dMData <- merge(
  x = dMData, y = lastDate,
  by = c('tagNumberCH'), all.x = TRUE
)
dMData <- dMData[ order(dMData$tagNumber,dMData$ageInSamples),]

dMData$emPerm <- 0

#dMData[dMData$tagNumber=='1748C762DD',]

# this puts a 1 at the last obervation
#dMData$emPerm[  dMData$everEmWest == 1 &
#  			dMData$sampleNum == dMData$last
#             ] <- 1

###
#  want to put 1 for emPerm on the line before the date of maxAntenna
###

dMData$emPerm <- ifelse( !is.na(dMData$maxAntennaDate2) &
                           dMData$maxAntennaDate2 > dMData$lastDate2 &  #antenna hit after last capture ?

                           dMData$maxAntennaDate2 > dMData$dateForEnv &
                           dMData$maxAntennaDate2 <= dMData$lagDateForEnv,
                         1,
                         0
)


dMData$emPermNA <- dMData$emPerm
#dMData$emPermNA <- ifelse( dMData$everEmWest != 1 &
#                           dMData$sampleNum >= dMData$last,
#                           NA, dMData$emPermNA)
dMData$emPermNA <- ifelse( dMData$everEmWest == 0 &
                             dMData$sampleNum > dMData$last,   #keep 0 for last obs. this works better with 'censored' in bugs. rather have phi kill fish than censored (when we're not modeling PIT ant det)
                           NA, dMData$emPermNA)

# indicator for permanent emigration for models where emigration is not estiamted
dMData$available01 <- ifelse( dMData$everEmWest == 1 &
                                dMData$sampleNum > dMData$last,
                              0, 1
)

#
dMData$antennaCountDudaGT1<- ifelse( dMData$maxAntennaDate2 > dMData$dateForEnv &
                                       dMData$maxAntennaDate2 <= dMData$lagDateForEnv,
                                     1, 0 )


#########################################################################
#########################################################################

#########################################################################
### AntennaTrib
#########################################################################
#
#

antennaTribs <- as.data.frame(read.csv(file="./antennaTribs.txt", header=TRUE))
# file created using access file query  fourRiversAts:antennaDuda query
# in C:\PITTAGMAIN\CMR Analyses\Hierach_Bugs\allSpp
#export query as text file without formatting
# then copy to felek

antennaTribs$date2 <- as.POSIXct(strptime(antennaTribs$date, format = "%m/%d/%Y"))
antennaTribs$tagNumberCH <- as.character(antennaTribs$tagNumber)

#countAntennaTribs <- ddply( antennaTribs[,c('tagNumber','river','date2')], .(tagNumber, river),
#							summarise,  ct=length(date2)
#		   				  )
#countAntennaTribs <- countAntennaTribs[ order(countAntennaTribs$tagNumber,countAntennaTribs$river),]


# count the number of time a fiSh waS oberved on trib antenna
# for each interval between dateForEnv and lagDataForEnv
# takeS a LONG time to run - need to improve on loopS
# AlSo need to makeSure dateForEnv iS acting ok.
#if (runAntennaCount) {

#dMData$antennaCountJimmy <- 0
#dMData$antennaCountMitchell <- 0

#for ( i in 1 : nrow( dMData ) ) {
#	hold <- antennaTribs[ antennaTribs$tagNumberCH == dMData$tagNumberCH[i], ]

#	cJimmy <- 0
#	cMitchell <- 0

#	if ( nrow(hold) > 0 &
#	     !is.na(dMData$river[i]) &
#	     !is.na(dMData$dateForEnv[ i ]) &
#	     !is.na(dMData$lagDateForEnv[ i ] ) ) {

#		for ( j in 1 : nrow(hold) ) {
#			if( hold$river[j] == 'WB JIMMY' &
#				hold$date2[j] > dMData$dateForEnv[ i ]  &
#				hold$date2[j] <= dMData$lagDateForEnv[ i ] ) cJimmy <- cJimmy + 1
#
#			if( hold$river[j] == 'WB MITCHELL' &
#				hold$date2[j] > dMData$dateForEnv[ i ]  &
#				hold$date2[j] <= dMData$lagDateForEnv[ i ] ) cMitchell <- cMitchell + 1
#     }

#	}

#	dMData$antennaCountJimmy[ i ] <- cJimmy
#	dMData$antennaCountMitchell[ i ] <- cMitchell
#	if ( i %% 50 == 0 ) print( c( i,nrow(dMData),Sys.time() ) )
#}

#}

#thiS work faster
m <- merge ( x = dMData, y=antennaTribs[,c('river','tagNumberCH','date2')], by = c('tagNumberCH'), all.x=TRUE)
#m <- m[order(m$tagNumberCH,m$sampleNum),]
m$antennaHit <- ifelse( m$date2.y > m$dateForEnv &
                          m$date2.y <= m$lagDateForEnv,
                        1, 0 )

m$antennaHit[ is.na( m$antennaHit) ] <- 0

# too Slow
#mHit <- ddply( m, .(river.y, tagNumberCH, sampleNum), summarise, sum( antennaHit ), .progress = "text" )

#mHit <- tapply(m$antennaHit,list(m$river.y, m$tagNumberCH, m$sampleNum), sum)#, na.rm=TRUE)
mHit <- aggregate( m$antennaHit , by=list( m$river.y, m$tagNumberCH, m$sampleNum ),FUN=sum, na.rm =TRUE )

names(mHit) <- c('river','tagNumberCH','sampleNum','antennaCount')

rm(m)

mHitJimmy <- subset(mHit, river == 'WB JIMMY')
mHitMitchell <- subset(mHit, river == 'WB MITCHELL')
mHitOBear <- subset(mHit, river == 'WB OBEAR')

names(mHitJimmy) <- c('antennaTrib','tagNumberCH', 'sampleNum', 'antennaCountJimmy')
names(mHitMitchell) <- c('antennaTrib','tagNumberCH', 'sampleNum', 'antennaCountMitchell')
names(mHitOBear) <- c('antennaTrib','tagNumberCH', 'sampleNum', 'antennaCountOBear')


dMData <- merge ( x = dMData, y=mHitJimmy[,2:4], by = c('tagNumberCH', 'sampleNum'), all.x=TRUE)
dMData <- merge ( x = dMData, y=mHitMitchell[,2:4], by = c('tagNumberCH', 'sampleNum'), all.x=TRUE)
dMData <- merge ( x = dMData, y=mHitOBear[,2:4], by = c('tagNumberCH', 'sampleNum'), all.x=TRUE)


dMData$antennaCountJimmy[ is.na( dMData$antennaCountJimmy ) ] <- 0
dMData$antennaCountMitchell[ is.na( dMData$antennaCountMitchell ) ] <- 0
dMData$antennaCountOBear[ is.na( dMData$antennaCountOBear ) ] <- 0

dMData$antennaCountJimmyGT1 <- ifelse(dMData$antennaCountJimmy >0, 1, 0)
dMData$antennaCountMitchellGT1 <- ifelse(dMData$antennaCountMitchell >0, 1, 0)
dMData$antennaCountOBearGT1 <- ifelse(dMData$antennaCountOBear >0, 1, 0)

# get list of fish that moved from OBear to WB - this is better than fromTo becasue is doesn't rely on consec samples
#could do the same thing for other movements...
wBAndOB <- ddply(dMData[!is.na(dMData$river),], .(tagNumberCH), function(x) {any(x$river == c('WEST BROOK')) * any(x$river == c('WB OBEAR'))})
wBAndOBTN <- wBAndOB[wBAndOB$V1 == T,1]

dMData$wBAndOB <- ifelse( dMData$tagNumberCH %in% wBAndOBTN, 1,0 )


#######################################################
# Fill in river observations for missing occasions
# observations per individual in each river

countRivsInd  <-  (ftable( dMData$tagNumberCH, dMData$river ))
countRivsInd2 <- melt( countRivsInd ) # / rowSums( countRivsInd ) )
names(countRivsInd2) <- c('tagNumber','river','count')

# get river with most observations per fish - use as default before capture
r=ddply(countRivsInd2, .(tagNumber), function(x) { c( max(x$count), which.max(x$count), maxRiver = x$river[which.max(x$count)] ) }   )

dMData$riverConsec <- NA
# loops would take 10 hours
# for ( i in 1:nrow(dMData) ) {
#   if ( dMData$sampleNum[ i ] < dMData$first[ i ] ) {
#     dMData$riverConsec[ i ] <- r$maxRiver[ dMData$tagNumber[ i ] ]
#     print(c(i,'first', r$maxRiver[ dMData$tagNumber[ i ] ]))
#   }
#
#   else if ( !is.na( dMData$river[ i ] ) ) {
#     dMData$riverConsec[ i ] <-  dMData$river[ i ]
#     print(c(i,'notFirstRiver Yes',dMData$riverConsec[ i ]))
#   }
#
#   else if ( is.na( dMData$river[ i ] ) ) {
#     dMData$riverConsec[ i ] <-  dMData$riverConsec[ i - 1 ]
#     print(c(i,'notFirst River No',dMData$riverConsec[ i ]))
#   }
#
# }

# NOTE - this lags backwards
dMData <- addLagged(data = dMData, individual = "tagNumberCH", time = "sampleNum",
                    lag = c("river"), pf = 'lagBack')
#names(dMData)[names(dMData)=="lagRiver"] <- 'lagBackRiver'
# Fills in most common river per tagNumber (r$maxRiver[ dMData$tagNumber ]) before first sample
# Keeps river for obscerved occasions and grabs lagged (backwards) river for first uncaptured occasion
dMData$riverConsec <- ifelse( dMData$sampleNum < dMData$first,
                              r$maxRiver[ dMData$tagNumber ],
                              ifelse( !is.na( dMData$river ),
                                      dMData$river,
                                      dMData$lagBackRiver
                              )
)


#########################################################################
#########################################################################

# do thiS here becauSe need conec Sample numS and will give orig Sample numS when done in pheno
dMData <- dMData[ order(dMData$tagNumberCH,dMData$sampleNum),]

dMData <- addLaggedForward(data = dMData, individual = "tagNumberCH", time = "sampleNum",
                           lag = c("sampleNum","riverN","riverM","length","weight","date2"))
# flowing 2 lines were above in pheno - moved out of phenoe because they were based on consec dates and not sample nums
#pheno <- addLaggedForward(data = pheno, individual = "tagNumber", time = "date2",
#           lag = c("riverN","riverM","length","weight","date2"))#, "sampleNum"))

#########################################################################


#'dateTime...' are the date cols in the pheno data
(before <- Sys.time())
dMData <- addEnvironmentalData2(
  dMData, envData, lagDateTimeCol = "dateForEnv",
  dateTimeCol = "lagDateForEnv",
  dateCol = "date2", temperatureColumn = "temperature",
  dischargeColumn = "discharge" )
(after <-Sys.time())
difftime(after,before)

# fill in 'samples' before the study starts with means

envMeans <- aggregate((dMData[,c('intervalLength','maxT','meanT','medianT','minT','sdT','skewT',
                                 'maxD','meanD','medianD','minD','sdD','skewD') ] ),
                      by=list( dMData$drainage,dMData$season ),FUN=mean, na.rm =TRUE )
names(envMeans)[1:2] <- c('drainage','season')

##############################################################################
##############################################################################
#### Temporary fix to fill in mean data for augmented back rows ##############
##############################################################################
dMData$seasonMeanIntLen<-envMeans$intervalLength[match(dMData$season,envMeans$season)]
dMData$seasonMeanT <- envMeans$meanT[ match(dMData$season,envMeans$season) ]
dMData$seasonMeanD<-envMeans$meanD[match(dMData$season,envMeans$season)]
dMData$fullMeanIntLen<-dMData$intervalLength
dMData$fullMeanT<-dMData$meanT
dMData$fullMeanD<-dMData$meanD
dMData$fullMeanIntLen[which(is.na(dMData$intervalLength))]<-dMData$seasonMeanIntLen[which(is.na(dMData$intervalLength))]
dMData$fullMeanT[which(is.na(dMData$meanT))]<-dMData$seasonMeanT[which(is.na(dMData$meanT))]
dMData$fullMeanD[which(is.na(dMData$meanD))]<-dMData$seasonMeanD[which(is.na(dMData$meanT))]
##############################################################################
##############################################################################
##############################################################################
##############################################################################

dMData <- dMData[ order(dMData$tagNumber,dMData$sampleNum),]

##############################################################################

evalList <- list(firstObsRows=firstObsRows,nFirstObsRows=nFirstObsRows,
                 lastObsRows=lastObsRows,nLastObsRows=nLastObsRows,
                 evalRows=evalRows,nEvalRows=nEvalRows,
                 evalJSRows=evalJSRows,nEvalJSRows=nEvalJSRows,
                 summerObsRows=summerObsRows,nSummerObsRows=nSummerObsRows,
                 nonSummerObsRows=nonSummerObsRows, nNonSummerObsRows=nNonSummerObsRows,
                 summerAIS=summerAIS, nonSummerAIS=nonSummerAIS)

# means for standardizing

lengthStd <- tapply(dMData$length,dMData$ageInSamples,mean, na.rm=TRUE)
lengthStd <- ifelse( is.na(lengthStd), as.numeric(dimnames(lengthStd)[[1]])*15, lengthStd )     # need to do this for the SR, so fish in AIS 10,13,14

stdList <-
  list(
    lengthStd = lengthStd,
    tempStd = tapply(dMData$fullMeanT,dMData$season,mean),
    flowStd = tapply(dMData$fullMeanD,dMData$season,mean)
  )

# for cohort by ais means
nAgeInSamples = length(table(dMData$ageInSamples))
nCohorts = length(table(dMData$cohort))
cohortIndex<-as.numeric(as.factor(dMData$cohort))

aisMeans <- tapply(dMData$length,dMData$ageInSamples, mean, na.rm=TRUE)
aisCohortMeans <- tapply(dMData$length,list(dMData$cohort,dMData$ageInSamples), mean, na.rm=TRUE)

temp=array(NA,c(8,14))
for(c in 1:nCohorts){
  for(a in 1:nAgeInSamples){
    if (is.na(aisCohortMeans[c,a])) aisCohortMeans[c,a] <- aisMeans[a]
  }
}

stdList_cohort <-
  list(
    lengthStd = aisCohortMeans,
    tempStd = tapply(dMData$fullMeanT,dMData$season,mean),
    flowStd = tapply(dMData$fullMeanD,dMData$season,mean)
  )

##########################################################
# do counts for stdN in the growthSurvivalMove model. need counts of fish from all cohorts for estimate of N
pheno2Long$riverOrdered <- factor(pheno2Long$river,levels=c('WEST BROOK','WB JIMMY','WB MITCHELL','WB OBEAR'), ordered=T)
pheno2Long$season2 = as.numeric(as.character(pheno2Long$seasonConsec))

count2 <- melt( ftable(pheno2Long[pheno2Long$enc == 1 ,c('seasonConsec','yearConsec','riverOrdered')]) )
names(count2) <- c('season','year','riverName','count')
count2$river <- as.numeric(count2$riverName) + 1
count2$year <- as.numeric(as.character(count2$year))
count2$season <- as.numeric(as.character(count2$season))
fillIn5 <- cbind(unique(count2[,c('season','year')]),'dead',0,1)
names(fillIn5) <- names(count2)
count2a <- rbind(count2,fillIn5)
count2a$count<- ifelse( count2a$river == 1, rbinom(  sum(( count2a$river == 1)*1),20,0.75 ), count2a$count )  # put in dummy numbers so don't get a sd of 0 calculated in bugs


meanCount <- ddply( count2a, .( season,river ), summarise, mean=mean(count), sd=sd(count))
count3 <- merge(
  x = count2a, y = meanCount,
  #by = c('time','river'),
  all.x = TRUE
)

count3$countAdj <- (count3$count-count3$mean)/count3$sd
count3$countAdj[count3$river == 1]  <- 0
countForN2 <- melt(count3, id.vars=c('season','year','river'), measure.vars=c('count','mean','sd'))
#  countForN2$countAdj <- countForN2$value
countForN2$year2 <- as.numeric(as.character(countForN2$year)) - min(countForN2$year) +1
#  countForN2$season <- as.numeric(as.character(countForN2$season))
countForN2 <- countForN2[order(countForN2$variable,countForN2$river,countForN2$year,countForN2$season),]

nRivers <- length(unique(countForN2 $river)) -1 # -1 for NA
nYears = max(countForN2 $year)-min(countForN2 $year)+1

countForN <- array(countForN2$value[countForN2$variable == 'count'],c(4,nYears,nRivers+1))
meanForN <- array(countForN2$value[countForN2$variable == 'mean'],c(4,nRivers+1))
sdForN <- array(countForN2$value[countForN2$variable == 'sd'],c(4,nRivers+1))
sdForN[sdForN == 0 ]  <- 1 # so we don't divide by 0

statsForN <- list(
  countForN = countForN,
  meanForN = meanForN,
  sdForN = sdForN,
  nYears = nYears,
  minYear = min(countForN2$year),
  maxYear = max(countForN2$year)

)

##########################################################
#stop('temp')
##############################################################################
# output dMData  to species_DMData_river.RData
##############################################################################

directory <- '/home/ben/allSpp/'
directory <- tempfile( pattern="output-", tmpdir ='.', fileext='-dMData')
dir.create(directory)
file.copy(from='./makeDMData.R', to=paste(directory,'makeDMData.R',sep='/'))
writeLines(text=directory, con='./latest_directory')
print(directory)

fileName <- paste('dMDataOut',speciesSubset,subsetDMdataCohortMin,'_', subsetDMdataCohortMax,'.RData', sep='')

#save.image(paste(directory,fileName, sep=''))
save(dMData, evalList, stdList, stdList_cohort, statsForN,
     file = paste(directory,fileName, sep='/'))

##############################################################################
##############################################################################
# end of data prep
##############################################################################
##############################################################################

stop('Done') # stop when this file is 'sourced'

# data checks
table(pheno2$river)
length(unique(pheno2$tagNumber))

table(pheno2Long$river[pheno2Long$enc==1])
table(pheno2Long$river[pheno2Long$available==1])
length(unique(pheno2Long$tagNumber))

table(dMData$river[dMData$enc==1])
length(unique(dMData$tagNumber))

table(dMData$river[dMData$enc==1],dMData$ageInSamples[dMData$enc==1])
ftable(dMData$river[dMData$enc==1],
       dMData$cohort[dMData$enc==1],
       dMData$ageInSamples[dMData$enc==1]
)


table(dMData$river[dMData$enc==1],
      dMData$sampleNum[dMData$enc==1]
)


table(dMData$ageInSamples[dMData$enc==1],
      dMData$cohort[dMData$enc==1]
)

tapply(dMData$length[dMData$enc==1],
       list(dMData$ageInSamples[dMData$enc==1],
            dMData$cohort[dMData$enc==1]),
       function(x){mean(x,na.rm=TRUE)})




################################################################################
# below sets up design matrix - not using for now because it's slow- may return to....
################################################################################

# set intercept levels
#dMData$river <- relevel(dMData$river, ref='WEST')
#dMData$season <- relevel(dMData$season, ref='3')
#dMData$age <- relevel(dMData$age, ref='1')

# get indices for random effects
dMData$fullRandom <- paste(dMData$drainage,dMData$season, dMData$age, dMData$cohort,sep=":")
dMData$fullRandomIndex <- as.numeric(as.factor(dMData$fullRandom))

randomIndex <- unique(cbind(dMData$fullRandomIndex,dMData$drainage,dMData$season, dMData$age, dMData$cohort))
randomIndex <- randomIndex[order(randomIndex[,1]),]

#if put length in dM, it'll exclude all the unobserved length rows
dM1 <- model.matrix(enc ~ ( season * age )^2 , dMData )

dM <- cbind(dM1,dMData$fullRandomIndex)
dimnames(dM)[[2]][ncol(dM)] <- "randomIndex"

#exclude colSums ==0
dM <- dM[ , colSums(dM) > 0 ]

#get length vector for input - this is ordered same as dM
lengthIn <- as.numeric(dMData$length)
#replace NAs with 0s. this will mult sizeBetas by 0 when not observed
lengthIn[ is.na( lengthIn )] <- 0

table(dMData$age)
