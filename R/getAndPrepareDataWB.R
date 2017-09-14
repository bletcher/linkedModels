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

#'Clean data from the PIT tag database
#'
#'@param d dataframe created with getCoreData()
#'@param drainageIn Which drainage, "west" or "stanley"
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
    d$riverOrdered <- factor(d$river,levels=c('mainstem', 'west', 'east'),labels = c('mainstem', 'west', 'east'), ordered=T)
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


#'Merge sites table
#'
#'@param d dataframe created with getCoreData()
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

mergeSites <- function(d,drainageIn){
  sites <- getSites(drainageIn)
  # merge in riverMeter for sections
  d <- left_join(d, sites, by = c("river","section","area"))
  d$riverMeter <- ifelse( d$survey == "shock" | d$survey == "portableAntenna", d$river_meter, d$riverMeter )
  return(d)
}

#'Show minimal data
#'
#'@param d dataframe created with getCoreData()
#'@return a data frame with key selected columns
#'@export

minimalData <- function(d){
  d %>% select(tag,detectionDate,sampleNumber,riverOrdered,observedLength,survey,enc,knownZ,grLength)
}
############# 2_prepare data
# install dev version to fix NA problem with lag()
#if (packageVersion("devtools") < 1.6) {
#  install.packages("devtools")
#}
#devtools::install_github("hadley/lazyeval")
#devtools::install_github("hadley/dplyr")
