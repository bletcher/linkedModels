#'Extract data from the PIT tag database
#'
#'@param drainage Which drainage, "west" or "stanley"#'Extract data from the PIT tag database
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

getCoreData <- function(drainage = "west"){

  cdWB <- createCoreData(sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
                         whichDrainage = drainage,
                         columnsToAdd=c("sampleNumber","river","riverMeter","survey",'observedLength','observedWeight')) %>%
    addTagProperties( columnsToAdd = c("cohort","species","dateEmigrated","sex","species")) %>%
    dplyr::filter( !is.na(tag), area %in% c("trib","inside","below","above") ) %>%
    createCmrData( maxAgeInSamples = 20, inside = F, censorDead = F, censorEmigrated = T) %>%
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
    minYear = 2002
  }
  else if(drainageIn == "stanley"){
    maxSectionNum <- 50
    d$riverOrdered <- factor(d$river,levels=c('mainstem', 'west', 'east'),labels = c('mainstem', 'west', 'east'), ordered=T)
    minYear = 2006
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
            maxSample = max(sampleNumber),
            minYear = minYear) %>%
    ungroup()

  d$moveDir <- ifelse( d$section == d$lagSection, 0, ifelse( d$section > d$lagSection, 1,-1 ) )
  d$sampleInterval <- as.numeric(d$lagDetectionDate - d$detectionDate)

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
  d %>% select(tag,detectionDate,sampleNumber,riverOrdered,observedLength,
               survey,enc,knownZ,grLength,lagDetectionDate,lagObservedLength)
}


#'Get propoortion of sections sampled by river,sample
#'
#'@param nSeasons,nRivers,nYears counts
#'@return a list of propSampleDATA, and zeroSectionsDATA
#'@export

getPropSampled <- function(nSeasons,nRivers,nYears){

  #check data

#  tmp <- dddD %>%
#           group_by(riverOrdered,year,season,section) %>%
#           summarize( s = sum(enc) ) %>%
#           filter( s == 0 )
#  table(tmp$riverOrdered,tmp$year,tmp$season)

  # propSamp - proportion of each season,river,year combo that got sampled (proportion of sctions sampled)
  propSampledDATA <- array( 1, c(nSeasons,nRivers,nYears) ) #season, river year

  propSampledDATA[ c(1,4),1:4,1 ] <- 0     #all spring and winter samples in 2002
  propSampledDATA[ 2,2:3,1 ] <- 0          #J and M summer samples in 2002
  propSampledDATA[ 4,1,2 ] <- 30/47        #WB winter sample in 2003
  propSampledDATA[ 4,1,3 ] <- 3/47         #WB winter sample in 2004
  propSampledDATA[ 4,1,4 ] <- 0            #WB winter sample in 2005
  propSampledDATA[ 4,1,6 ] <- 0            #WB winter sample in 2007
  propSampledDATA[ 4,1,11 ] <- 0           #WB winter sample in 2012
  propSampledDATA[ 1,1,14 ] <- 0           #WB spring sample in 2015
  propSampledDATA[ c(3,4),c(1,4),14 ] <- 0 #all fall and winter samples in 2015

  # zeroSectionsDATA - completely unsampled winter samples. not including samples before season 4,year 1 because we didn't want to rewrite the meanPhiS34 indexing [mostly noise ni these estimates anyway]
  zeroSectionsDATA <- array( 0, c(nSeasons,nRivers,nYears) ) #season, river year

  zeroSectionsDATA[ 3:4,1:4,1 ] <- 1          #all winter samples in 2002
  zeroSectionsDATA[ 3:4,1,4 ] <- 1            #WB winter sample in 2005
  zeroSectionsDATA[ 3:4,1,6 ] <- 1            #WB winter sample in 2007

  return(list(propSampledDATA = propSampledDATA,zeroSectionsDATA = zeroSectionsDATA))
}


#'Get cutoFFYoy data
#'
#'@param d dataframe created with getCoreData()
#'@param dr Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

getYOYCutoffs <- function(d,dr){

  ####################################
  # create a data frame of max lengths for YOYs from Matt,..

  yoySizeLimits <- read.csv(file='./data/yoySizeLimits.csv', header=T)
  yoySizeLimits$year <- yoySizeLimits$YOS
  yoySizeLimits$river <- yoySizeLimits$River

  snOrigSeason <- data.frame(unique(cbind(d$sampleNumber,d$season)))
  names(snOrigSeason) <- c('Sample', 'season')

  y2 <- merge(x=yoySizeLimits, y=snOrigSeason, by='Sample', all.x=T)
  y2$year <- y2$YOS
  y2$river <- tolower(y2$river)

  yTemplate <- data.frame( year=rep(2002:max(d$year), each=5*4), river= rep(c('DEAD',"wb obear","west brook","wb jimmy","wb mitchell"),each=4), season=1:4 )
  y <- merge( x=yTemplate, y=y2, by=c('river','year','season'), all.x=T)

  y$maxLength <- ifelse( is.na(y$Max.Length), 90, y$Max.Length) # fill in missing obs with 90
  y$maxLength <- ifelse( y$season == 1, 50, y$maxLength )

  # some visual fixes
  y$maxLength[ y$year==2003 & y$season==3 & y$river == 'wb jimmy' ] <- 85
  y$maxLength[ y$year==2011 & y$season==2 & y$river == 'wb obear' ] <- 70
  y$maxLength[ y$year==2010 & y$season==3 & y$river == 'wb obear' ] <- 72
  y$maxLength[ y$year==2011 & y$season==3 & y$river == 'west brook' ] <- 100
  y$maxLength[ y$year==2011 & y$season==4 & y$river == 'west brook' ] <- 110


  y$riverOrdered <- factor(y$river, levels=c('DEAD','west brook','wb jimmy','wb mitchell','wb obear'), ordered=T)
  y <- y[ order(y$year,y$riverOrdered,y$season),]
  cutoffYOYDATA <- array( y$maxLength, c(4,5,11) )

  save(cutoffYOYDATA,file = './data/cutoffYOYDATA.RData')
  #######################################

  # river <- 'wb mitchell'#,
  # river <- 'west brook'#'wb obear'##'#''#'wb obear' #
  # river <- 'wb obear'##'#''#'wb obear' #
  # river <- 'wb jimmy'
  #
  # ggplot( cd %>% filter(drainage == dr), aes(observedLength) ) +
  #   geom_histogram( binwidth=3 )+
  #   geom_vline(aes(xintercept=maxLength), y[y$year>=2002 & y$river==river,])+
  #   facet_grid(season~year ) +
  #   ggtitle(river)

  # incorporating the spring sample 1+ fish into the yoy category
  yy <- y

  yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'west brook' ] <- 102
  yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'west brook' ] <- 105
  yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'west brook' ] <- 115
  yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'west brook' ] <- 100
  yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'west brook' ] <- 109
  yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'west brook' ] <- 125
  yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'west brook' ] <- 127
  yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'west brook' ] <- 114
  yy$maxLength[ yy$year==2013 & yy$season==1 & yy$river == 'west brook' ] <- 118
  yy$maxLength[ yy$year==2014 & yy$season==1 & yy$river == 'west brook' ] <- 116
  yy$maxLength[ yy$year==2015 & yy$season==1 & yy$river == 'west brook' ] <- 114

  yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'wb jimmy' ] <- 82
  yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'wb jimmy' ] <- 87
  yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'wb jimmy' ] <- 93
  yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'wb jimmy' ] <- 90
  yy$maxLength[ yy$year==2013 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$year==2014 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$year==2015 & yy$season==1 & yy$river == 'wb jimmy' ] <- 93

  yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'wb mitchell' ] <- 95
  yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'wb mitchell' ] <- 110
  yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'wb mitchell' ] <- 109
  yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'wb mitchell' ] <- 107
  yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'wb mitchell' ] <- 88
  yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'wb mitchell' ] <- 83
  yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'wb mitchell' ] <- 102
  yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'wb mitchell' ] <- 114
  yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'wb mitchell' ] <- 132
  yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$year==2013 & yy$season==1 & yy$river == 'wb mitchell' ] <- 90
  yy$maxLength[ yy$year==2014 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$year==2015 & yy$season==1 & yy$river == 'wb mitchell' ] <- 78

  yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'wb obear' ] <- 95
  yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'wb obear' ] <- 90
  yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'wb obear' ] <- 83
  yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'wb obear' ] <- 106
  yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$year==2013 & yy$season==1 & yy$river == 'wb obear' ] <- 100
  yy$maxLength[ yy$year==2014 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$year==2015 & yy$season==1 & yy$river == 'wb obear' ] <- 84

  yy$riverOrdered <- factor(yy$river, levels=c('DEAD','west brook','wb jimmy','wb mitchell','wb obear'), ordered=T)
  yy <- yy[ order(yy$year,yy$riverOrdered,yy$season),]
  cutoffYOYInclSpring1DATA <- array( yy$maxLength, c(4,5,16) )

  save(yy,cutoffYOYInclSpring1DATA,file='./data/cutoffYOYInclSpring1DATA.RData')

  #check cutOffs
  #  yy[  yy$season==4 & yy$river == 'wb obear' ,c('year','maxLength')]

  # river <- 'wb mitchell'#,
  # river <- 'west brook'#'wb obear'##'#''#'wb obear' #
  # river <- 'wb obear'##'#''#'wb obear' #
  # river <- 'wb jimmy'
  #
  # ggplot( cd[cd$river==river & !is.na(cd$season),], aes(observedLength) ) +
  #   geom_histogram( binwidth=3 )+
  #   geom_vline(aes(xintercept=maxLength), yy[yy$year>=2002 & yy$river==river,])+
  #   facet_grid(season~year ) +
  #   ggtitle(river)

}

############# 2_prepare data
# install dev version to fix NA problem with lag()
#if (packageVersion("devtools") < 1.6) {
#  install.packages("devtools")
#}
#devtools::install_github("hadley/lazyeval")
#devtools::install_github("hadley/dplyr")
