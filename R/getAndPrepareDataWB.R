# data end up here: file='/home/ben/dataForModels/cdForMovementModelWB.RData')
# to avoid problems with .gitignore not ignoring .RData files



# to rebuild library:
# install.packages("devtools")
# devtools::install_github('Conte-Ecology/westBrookData/getWBData')

library(getWBData)
library(dplyr)
library(dbplyr)
library(lubridate)


#'Extract data from the PIT tag database
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

getCoreData <- function(drainage = "west"){

  cdWB <- createCoreData(sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
                         whichDrainage = drainage,
                         columnsToAdd=c("sampleNumber","river","riverMeter","survey",'observedLength','observedWeight')) %>%
    addTagProperties( columnsToAdd=c("cohort","species","dateEmigrated","sex","species")) %>%
    dplyr::filter( !is.na(tag), area %in% c("trib","inside","below","above") ) %>%
    createCmrData( maxAgeInSamples=20, inside=F) %>%
    addSampleProperties() %>%
    addEnvironmental() %>%
    addKnownZ() %>%
    fillSizeLocation(size = F) #assumes fish stay in same location until observed elsewhere
}


#'Get data from sites table
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

getSites <- function(drainageIn = "west"){
  # get sites table
  sitesIn <- data.frame(tbl(conDplyr,"data_sites") )
  sites <- sitesIn %>% filter(is.na(quarter) & !is.na(quarter_length) & drainage == drainageIn) %>% select(-quarter)
  sites$section <- as.numeric(sites$section)
  return(sites)
}

############# 2_prepare data
# install dev version to fix NA problem with lag()
#if (packageVersion("devtools") < 1.6) {
#  install.packages("devtools")
#}
#devtools::install_github("hadley/lazyeval")
#devtools::install_github("hadley/dplyr")

# sample 2.5 = fyke net. not sure if we should keep it - prob with obs model
#cdWB <-  filter(cdWB, sampleNumber != 2.5 & sampleNumber != 10.1)

#'Clean data from the PIT tag database
#'
#'@param d dataframe created with getCoreData()
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

cleanData <- function(d,drainageIn){
  # some formatting fixes
  d$sectionOriginal <- d$section
  d$section <- as.numeric( d$section )

  if(drainageIn == "west") {
    maxSectionNum <- 47
    d$riverOrdered <- factor(d$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)
  }
  else if(drainageIn == "stanley"){
    maxSectionNum <- 52

    # need to fix
    d$riverOrdered <- factor(d$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)
  }

  d$inside <- ifelse( d$section %in% 1:maxSectionNum | d$survey == "stationaryAntenna", T, F )

  d$year <- year(d$detectionDate)
  d$yday <- yday(d$detectionDate)

  d <- d %>%
    group_by(tag) %>%
    # arrange(tag,sampleNumber) %>%
    mutate( lagSection = lead(section),
            distMoved = section - lagSection,
            lagObservedWeight = lead(observedWeight),
            lagObservedLength = lead(observedLength),
            grWeight = exp(lagObservedWeight - observedWeight)/as.numeric((lagDetectionDate - detectionDate)),
            grLength = (lagObservedLength - observedLength)/as.numeric((lagDetectionDate - detectionDate)),
            minSample = min(sampleNumber),
            maxSample = max(sampleNumber)) %>%
    ungroup()

  d$moveDir <- ifelse( d$section == d$lagSection, 0, ifelse( d$section > d$lagSection, 1,-1 ) )
  return(d)
}



