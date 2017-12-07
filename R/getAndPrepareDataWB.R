#'Extract data from the PIT tag database
#'
#'@param drainage Which drainage, "west" or "stanley"#'Extract data from the PIT tag database
#'@return a data frame
#'@export

getCoreData <- function(drainage = "west"){

  cdWB <- createCoreData(sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
                         whichDrainage = drainage,
                         columnsToAdd=c("sampleNumber","river","riverMeter","survey",'observedLength','observedWeight')) %>%
    addTagProperties( columnsToAdd = c("cohort","species","dateEmigrated","sex","species")) %>%
    dplyr::filter( !is.na(tag), area %in% c("trib","inside","below","above"), !is.na(sampleNumber) ) %>%
    createCmrData( maxAgeInSamples = 20, inside = F, censorDead = F, censorEmigrated = T) %>%
    addSampleProperties() %>%
    addEnvironmental() %>%
    addKnownZ() %>%
    fillSizeLocation(size = F) #assumes fish stay in same location until observed elsewhere
}

#'Get counts and summed mass of all fish, including untagged
#'
#'@param drainage Which drainage, "west" or "stanley"#'Extract data from the PIT tag database
#'@return a data frame
#'@export

getCountsAllFish <- function(drainage = "west", filteredAreas = c("inside","trib")){

  cdWBAll <- createCoreData(sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
                            whichDrainage = drainage,
                            columnsToAdd=c("sampleNumber","river","riverMeter","survey",'observedLength','observedWeight'),
                            includeUntagged = T) %>%
    addTagProperties( columnsToAdd = c("cohort","species","dateEmigrated","sex","species")) %>%
    dplyr::filter( area %in% c("trib","inside","below","above"), !is.na(sampleNumber) ) %>%
    # createCmrData( maxAgeInSamples = 20, inside = F, censorDead = F, censorEmigrated = T) %>%
    addSampleProperties() %>%
    # addEnvironmental() %>%
    # addKnownZ() %>%
    # fillSizeLocation(size = F) #assumes fish stay in same location until observed elsewhere
    filter( species %in% c('bkt','bnt','ats'))

    cdWBAll$riverOrdered <- factor(cdWBAll$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)
    cdWBAll$riverN <- as.numeric(cdWBAll$riverOrdered)

    counts <- cdWBAll %>%
      filter( area %in% filteredAreas ) %>%
      group_by( species,river,season,year ) %>%
      summarize( nAllFishBySpecies = n(),
                 massAllFishBySpecies = sum(observedWeight,na.rm=T))

    #ggplot(counts, aes(year,nAllFishBySpecies,color=species)) +
    # ggplot(counts, aes(year,massAllFishBySpecies,color=species)) +
    #   geom_point() +
    #   geom_line() +
    #   labs( x = "Year", y = "Estimated count") +
    #   facet_grid(season ~ riverN, scales = "free")

    return(counts)
}

#'Add counts for all fish to cd
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame with vars nAllFishBySpecies, massAllFishBySpecies and nAllFish, massAllFish
#'@export
#'
addCounts <- function(cd,drainage){
  allFishBySpecies <- getCountsAllFish(drainage)
  cd <- cd %>% left_join(allFishBySpecies)

  allFish <- allFishBySpecies %>%
    group_by(river,season,year) %>%
    summarize( nAllFish = sum(nAllFishBySpecies, na.rm=T),
               massAllFish = sum(massAllFishBySpecies, na.rm=T))
  cd <- cd %>% left_join(allFish)

  save(allFishBySpecies,allFish, file = './data/out/countsAndMasses.RData')

  return(cd)
}

#'Get pass data from raw data table
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

addNPasses <- function(cd,dr){

  reconnect()

  nPasses <- tbl(conDplyr,"data_tagged_captures") %>%
    filter(drainage == dr) %>%
    dplyr::select(river,sample_number,pass) %>%
    distinct() %>%
    collect() %>%
    arrange(sample_number,river) %>%

    group_by(river,sample_number) %>%
    summarize( nPasses = max(pass,na.rm=T) ) %>%
    rename( sampleNumber = sample_number )


  cd <- left_join( cd,nPasses, by = c('river',"sampleNumber") )
  cd$nPasses <- ifelse( is.na(cd$nPasses) & cd$proportionSampled == 0, 1, cd$nPasses )
  #coreData[is.na(nPasses)&proportionSampled==0,nPasses:=1] #when proportionSampled==0 no pass info
  return(cd)
}


#'Get data from sites table
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

getSites <- function(drainageIn = "west"){
  # get sites table
  sitesIn <- data.frame(tbl(conDplyr,"data_sites") )
  sites <- sitesIn %>% filter(is.na(quarter) & !is.na(quarter_length) & drainage == drainageIn) %>% dplyr::select(-quarter)
  sites$section <- as.numeric(sites$section)
  return(sites)
}

#'Get counts of data from untagged fish
#'
#'@param d a dataframe
#'@return a data frame
#'@export

getCountOfUntagged <- function(drainage = 'west'){
  d1 <- getCoreDataAllFish(drainage) %>%
    filter(species %in% c('bkt','bnt','ats'),
           observedLength > 61,
           area %in% c('inside','trib') ) %>%
    mutate( year = year(detectionDate) )

  # d <- d1 %>%
  #   mutate( isTagged = ifelse( is.na(tag),0,1 ),
  #           season = seasonNumber ) %>%
  #   filter( observedLength > 61, area %in% c('inside','trib') ) %>%
  #   group_by( species,season,river,year,isTagged ) %>%
  #   summarise( n = n() )
  #
  # ggplot(d, aes(year,n,color=isTagged)) + geom_point()  + facet_grid(river~season+species)


  return(d1)
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
    minYear = min(d$year) #1997
  }
  else if(drainageIn == "stanley"){
    maxSectionNum <- 50
    d$riverOrdered <- factor(d$river,levels=c('mainstem', 'west', 'east'),labels = c('mainstem', 'west', 'east'), ordered=T)
    minYear = min(d$year) #2006
  }

  d$inside <- ifelse( d$section %in% 1:maxSectionNum | d$survey == "stationaryAntenna", T, F )

  d$year <- year(d$detectionDate)
  d$yday <- yday(d$detectionDate)

  dUntagged <- d %>%
    filter( is.na(tag) ) %>%
    mutate( minSample = min(sampleNumber),
            maxSample = max(sampleNumber),
            minYear = minYear,
            moveDir = 0,
            sampleInterval = 0)

  d <- d %>%
    filter( !is.na(tag) ) %>%
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

  d <- bind_rows( d,dUntagged )

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
  d %>% dplyr::select(tag,detectionDate,sampleNumber,riverOrdered,observedLength,
                      survey,enc,knownZ,grLength,lagDetectionDate,lagObservedLength)
}


#'Get proportion of sections sampled by river,sample
#'
#'@param nSeasons,nRivers,nYears counts
#'@return a list of propSampleDATA, and zeroSectionsDATA
#'@export

getPropSampled <- function(nSeasons,nRivers,nYears,minYear){

  #check data

   # tmp <- dddD[[2]] %>%
   #          group_by(riverOrdered,year,season,section) %>%
   #          summarize( s = sum(enc) ) %>%
   #          filter( s == 0 )
   # table(tmp$riverOrdered,tmp$year,tmp$season)
   # Season 4
   #              1998 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015
   # west brook     0    1    0   50   22   49   52    0   53    1    0    0   12   53    6   33    0
   # wb jimmy       0    0    0   14    0    0    0    0    0    0    2    0    0    0    0    0    0
   # wb mitchell    0    0    0   15    2    1    1    4    4    4    0    4    4    1    2    1    0
   # wb obear       0    0    0   15    6   10    0    0    1    0    0    0    2    0    0    0    0

  # propSamp - proportion of each season,river,year combo that got sampled (proportion of sctions sampled)
  propSampledDATA <- array( 1, c(nSeasons,nRivers,nYears) ) #season, river year

  #propSampledDATA[ c(1,4),1:4,2002 - minYear + 1 ] <- 0     #all spring and winter samples in 2002
  propSampledDATA[ c(1,4),2:4,2002 - minYear + 1 ] <- 0     #all spring and winter samples in 2002


  propSampledDATA[ 2,2:3,2002 - minYear + 1 ] <- 0          #J and M summer samples in 2002
  propSampledDATA[ 4,1,2003 - minYear + 1 ] <- 30/47        #WB winter sample in 2003
  propSampledDATA[ 4,1,2004 - minYear + 1 ] <- 3/47         #WB winter sample in 2004
  propSampledDATA[ 4,1,2005 - minYear + 1 ] <- 0            #WB winter sample in 2005
  propSampledDATA[ 4,1,2007 - minYear + 1 ] <- 0            #WB winter sample in 2007
  propSampledDATA[ 4,1,2012 - minYear + 1 ] <- 0           #WB winter sample in 2012
  propSampledDATA[ 1,1,2015 - minYear + 1 ] <- 0           #WB spring sample in 2015
  propSampledDATA[ c(3,4),c(1,4),2015 - minYear + 1 ] <- 0 #all fall and winter samples in 2015

  # zeroSectionsDATA - completely unsampled winter samples. not including samples before season 4,year 1 because we didn't want to rewrite the meanPhiS34 indexing [mostly noise ni these estimates anyway]

  # tmp2 <- dddD[[2]] %>%
  #          group_by(riverOrdered,year,season) %>%
  #          summarize( s = sum(enc) ) %>%
  #          filter( s == 0 )
  # # show unsampled river/year/season combos
  # table(tmp2$riverOrdered,tmp2$year,tmp2$season)
  #
  zeroSectionsDATA <- array( 0, c(nSeasons,nRivers,nYears) ) #season, river year

  zeroSectionsDATA[ 4,2:4,2002 - minYear + 1 ] <- 1          #trib winter samples in 2002
  zeroSectionsDATA[ 4,1,2005 - minYear + 1 ] <- 1            #WB winter sample in 2005
  zeroSectionsDATA[ 4,1,2007 - minYear + 1 ] <- 1            #WB winter sample in 2007
  zeroSectionsDATA[ 4,1,2012 - minYear + 1 ] <- 1            #WB winter sample in 2007
  zeroSectionsDATA[ 1,1,2015 - minYear + 1 ] <- 1            #WB spring sample in 2015
  zeroSectionsDATA[ c(3,4),1:4,2015 - minYear + 1 ] <- 1  #all fall and winter samples in 2015

  array2df(zeroSectionsDATA, label.x = "zero") %>%
    rename(est=est,season=d1,river=d2,year=d3)

  return(list(propSampledDATA = propSampledDATA,zeroSectionsDATA = zeroSectionsDATA))
}


#'List of unsampled samples for filtering out enc == 0 in the det model
#'
#'@param d dataframe created with getCoreData()
#'@return a data frame
#'@export

# zeroSectionsDATA - completely unsampled winter samples. not including samples before season 4,year 1 because we didn't want to rewrite the meanPhiS34 indexing [mostly noise ni these estimates anyway]

# tmp2 <- dddD[[2]] %>%
#          group_by(riverOrdered,year,season) %>%
#          summarize( s = sum(enc) ) %>%
#          filter( s == 0 )
# # show unsampled river/year/season combos
# table(tmp2$riverOrdered,tmp2$year,tmp2$season)
#
removeUnsampledRows <- function(d,drainage,removeIncomplete = F){

  nSeasons <- n_distinct(d$season, na.rm = T)
  nRivers <- n_distinct(d$riverOrdered, na.rm = T)
  nYears <- n_distinct(d$year, na.rm = T)
  minYear <- min(d$year)

  zeroSectionsDATA <- array( 0, c(nSeasons,nRivers,nYears) ) #season, river year

  if (drainage == 'west'){
    zeroSectionsDATA[ 4,2:4,2002 - minYear + 1 ] <- 1          #trib winter samples in 2002
    zeroSectionsDATA[ 4,1,2005 - minYear + 1 ] <- 1            #WB winter sample in 2005
    zeroSectionsDATA[ 4,1,2007 - minYear + 1 ] <- 1            #WB winter sample in 2007
    zeroSectionsDATA[ 4,1,2012 - minYear + 1 ] <- 1            #WB winter sample in 2012
    zeroSectionsDATA[ 1,1,2015 - minYear + 1 ] <- 1            #WB spring sample in 2015
    zeroSectionsDATA[ c(3,4),1:4,2015 - minYear + 1 ] <- 1  #all fall and winter samples in 2015

    if ( removeIncomplete ){
      # these had only a few sections done - remove for the detection model
      zeroSectionsDATA[ 4,1,2002 - minYear + 1 ] <- 1          #WB winter sample in 2002
      zeroSectionsDATA[ 1,1,2004 - minYear + 1 ] <- 1          #WB spring sample in 2004
    }
  }

  if (drainage == 'stanley'){

  }

  unsampledSamples <- array2df(zeroSectionsDATA, label.x = "unsampled01") %>%
    rename( season = d1, riverN = d2, yearN = d3 )

  d$yearN = d$year - minYear + 1
  d$riverN <- as.numeric(d$riverOrdered)

  d <- left_join(d,unsampledSamples) %>%
    filter( unsampled01 == 0 )

}

#'Get cutoFFYoy data
#'
#'@param d dataframe created with getCoreData()
#'@param dr Which drainage, "west" or "stanley"
#'@return a data frame and an array
#'@export

getYOYCutoffs <- function(d,dr = 'west'){

  ####################################
  # create a data frame of max lengths for YOYs from Matt, available here:
  # /home/projects/westbrook/dataIn/originalData/yoy_bins.csv

  # need to get riverLists for different watersheds...ie.e riverOrdered Not hardcoded

  y3 <- read.csv(file='./data/yoy_bins.csv', header = T)

  y2 <- y3 %>%
    rename( maxLength = Max.Length,
            minLength = Min.Length,
            drainage = Drainage,
            year = YOS,
            river = River,
            species = Species,
            sampleName = Sample,
            age = Age) %>%
    mutate( species = factor(tolower(species), levels = c('bkt','bnt','ats'), ordered = T),
            drainage = tolower(drainage),
            river = tolower(river),
            riverOrdered = factor(river, levels = c('west brook','wb jimmy','wb mitchell','wb obear'), ordered = T) ) %>%
    filter( drainage == dr,
            age == 0,
            species %in% speciesIn ) %>%
    arrange( species,riverOrdered,sampleName ) %>%
    dplyr::select( -river )

  # get season for each original sample
  snOrigSeason <- d %>% distinct(sampleName,sampleNumber,season) %>% arrange(sampleName) %>% mutate( sampleName = as.numeric(sampleName) ) #data.frame(unique(cbind(d$sampleName,d$sampleNumber,d$season)))
  #  snOrigSeason <- snOrigSeason %>% mutate( origSample = as.numeric(origSample),sample = as.numeric(sample),season = as.numeric(season)) %>%
  #                    arrange( origSample )

  y1 <- left_join(y2, snOrigSeason, by = 'sampleName')

  riverList <- unique(d$riverOrdered)
  nRivers <- n_distinct(d$riverOrdered, na.rm = T)
  nSeasons <- n_distinct(d$season, na.rm = T)
  nYears <- n_distinct(d$year, na.rm = T)
  nSpecies <- n_distinct(d$species, na.rm = T)

  yTemplate <- data.frame( species = rep(speciesIn,     each = nYears*nRivers*nSeasons),
                           year = rep(min(d$year):max(d$year), each = nRivers*nSeasons),
                           river = rep(riverList,                      each = nSeasons),
                           season =                                         1:nSeasons
  ) %>%
    mutate( riverOrdered = factor(river, levels = c('west brook','wb jimmy','wb mitchell','wb obear'), ordered = T) )

  y <- left_join( yTemplate, y1, by = c('riverOrdered','year','species','season') )

  y$maxLength <- ifelse( is.na(y$maxLength), 90, y$maxLength) # fill in missing obs with 90
  y$maxLength <- ifelse( y$season == 1, 50, y$maxLength )

  # some visual fixes
  y$maxLength[ y$species=='bkt' & y$year == 2003 & y$season == 3 & y$river == 'wb jimmy' ] <- 85
  y$maxLength[ y$species=='bkt' & y$year == 2011 & y$season == 2 & y$river == 'wb obear' ] <- 70
  y$maxLength[ y$species=='bkt' & y$year == 2010 & y$season == 3 & y$river == 'wb obear' ] <- 72
  y$maxLength[ y$species=='bkt' & y$year == 2011 & y$season == 3 & y$river == 'west brook' ] <- 100
  y$maxLength[ y$species=='bkt' & y$year == 2011 & y$season == 4 & y$river == 'west brook' ] <- 110

  y <- y[ order(y$species,y$year,y$riverOrdered,y$season),]
  cutoffYOYDATA <- array( y$maxLength, c(nSeasons,nRivers,nYears,nSpecies) )

  #  save(cutoffYOYDATA,file = './data/cutoffYOYDATA.RData')
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

  yy$maxLength[ yy$species=='bkt' & yy$year==2002 & yy$season==1 & yy$river == 'west brook' ] <- 102
  yy$maxLength[ yy$species=='bkt' & yy$year==2003 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$species=='bkt' & yy$year==2004 & yy$season==1 & yy$river == 'west brook' ] <- 105
  yy$maxLength[ yy$species=='bkt' & yy$year==2005 & yy$season==1 & yy$river == 'west brook' ] <- 115
  yy$maxLength[ yy$species=='bkt' & yy$year==2006 & yy$season==1 & yy$river == 'west brook' ] <- 100
  yy$maxLength[ yy$species=='bkt' & yy$year==2007 & yy$season==1 & yy$river == 'west brook' ] <- 109
  yy$maxLength[ yy$species=='bkt' & yy$year==2008 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$species=='bkt' & yy$year==2009 & yy$season==1 & yy$river == 'west brook' ] <- 125
  yy$maxLength[ yy$species=='bkt' & yy$year==2010 & yy$season==1 & yy$river == 'west brook' ] <- 127
  yy$maxLength[ yy$species=='bkt' & yy$year==2011 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$species=='bkt' & yy$year==2012 & yy$season==1 & yy$river == 'west brook' ] <- 114
  yy$maxLength[ yy$species=='bkt' & yy$year==2013 & yy$season==1 & yy$river == 'west brook' ] <- 118
  yy$maxLength[ yy$species=='bkt' & yy$year==2014 & yy$season==1 & yy$river == 'west brook' ] <- 116
  yy$maxLength[ yy$species=='bkt' & yy$year==2015 & yy$season==1 & yy$river == 'west brook' ] <- 114

  yy$maxLength[ yy$species=='bkt' & yy$year==2002 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bkt' & yy$year==2003 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bkt' & yy$year==2004 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bkt' & yy$year==2005 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2006 & yy$season==1 & yy$river == 'wb jimmy' ] <- 82
  yy$maxLength[ yy$species=='bkt' & yy$year==2007 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2008 & yy$season==1 & yy$river == 'wb jimmy' ] <- 87
  yy$maxLength[ yy$species=='bkt' & yy$year==2009 & yy$season==1 & yy$river == 'wb jimmy' ] <- 93
  yy$maxLength[ yy$species=='bkt' & yy$year==2010 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2011 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2012 & yy$season==1 & yy$river == 'wb jimmy' ] <- 90
  yy$maxLength[ yy$species=='bkt' & yy$year==2013 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$species=='bkt' & yy$year==2014 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$species=='bkt' & yy$year==2015 & yy$season==1 & yy$river == 'wb jimmy' ] <- 93

  yy$maxLength[ yy$species=='bkt' & yy$year==2002 & yy$season==1 & yy$river == 'wb mitchell' ] <- 95
  yy$maxLength[ yy$species=='bkt' & yy$year==2003 & yy$season==1 & yy$river == 'wb mitchell' ] <- 110
  yy$maxLength[ yy$species=='bkt' & yy$year==2004 & yy$season==1 & yy$river == 'wb mitchell' ] <- 109
  yy$maxLength[ yy$species=='bkt' & yy$year==2005 & yy$season==1 & yy$river == 'wb mitchell' ] <- 107
  yy$maxLength[ yy$species=='bkt' & yy$year==2006 & yy$season==1 & yy$river == 'wb mitchell' ] <- 88
  yy$maxLength[ yy$species=='bkt' & yy$year==2007 & yy$season==1 & yy$river == 'wb mitchell' ] <- 83
  yy$maxLength[ yy$species=='bkt' & yy$year==2008 & yy$season==1 & yy$river == 'wb mitchell' ] <- 102
  yy$maxLength[ yy$species=='bkt' & yy$year==2009 & yy$season==1 & yy$river == 'wb mitchell' ] <- 114
  yy$maxLength[ yy$species=='bkt' & yy$year==2010 & yy$season==1 & yy$river == 'wb mitchell' ] <- 132
  yy$maxLength[ yy$species=='bkt' & yy$year==2011 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2012 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2013 & yy$season==1 & yy$river == 'wb mitchell' ] <- 90
  yy$maxLength[ yy$species=='bkt' & yy$year==2014 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2015 & yy$season==1 & yy$river == 'wb mitchell' ] <- 78

  yy$maxLength[ yy$species=='bkt' & yy$year==2002 & yy$season==1 & yy$river == 'wb obear' ] <- 95
  yy$maxLength[ yy$species=='bkt' & yy$year==2003 & yy$season==1 & yy$river == 'wb obear' ] <- 90
  yy$maxLength[ yy$species=='bkt' & yy$year==2004 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$species=='bkt' & yy$year==2005 & yy$season==1 & yy$river == 'wb obear' ] <- 83
  yy$maxLength[ yy$species=='bkt' & yy$year==2006 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2007 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$species=='bkt' & yy$year==2008 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$species=='bkt' & yy$year==2009 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$species=='bkt' & yy$year==2010 & yy$season==1 & yy$river == 'wb obear' ] <- 106
  yy$maxLength[ yy$species=='bkt' & yy$year==2011 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2012 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2013 & yy$season==1 & yy$river == 'wb obear' ] <- 100
  yy$maxLength[ yy$species=='bkt' & yy$year==2014 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bkt' & yy$year==2015 & yy$season==1 & yy$river == 'wb obear' ] <- 84

  # bnt ###############################################################################################

  yy$maxLength[ yy$species=='bnt' & yy$year==2002 & yy$season==1 & yy$river == 'west brook' ] <- 100
  yy$maxLength[ yy$species=='bnt' & yy$year==2003 & yy$season==1 & yy$river == 'west brook' ] <- 105
  yy$maxLength[ yy$species=='bnt' & yy$year==2004 & yy$season==1 & yy$river == 'west brook' ] <- 105
  yy$maxLength[ yy$species=='bnt' & yy$year==2005 & yy$season==1 & yy$river == 'west brook' ] <- 115
  yy$maxLength[ yy$species=='bnt' & yy$year==2006 & yy$season==1 & yy$river == 'west brook' ] <- 100
  yy$maxLength[ yy$species=='bnt' & yy$year==2007 & yy$season==1 & yy$river == 'west brook' ] <- 111
  yy$maxLength[ yy$species=='bnt' & yy$year==2008 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$species=='bnt' & yy$year==2009 & yy$season==1 & yy$river == 'west brook' ] <- 125
  yy$maxLength[ yy$species=='bnt' & yy$year==2010 & yy$season==1 & yy$river == 'west brook' ] <- 132
  yy$maxLength[ yy$species=='bnt' & yy$year==2011 & yy$season==1 & yy$river == 'west brook' ] <- 110
  yy$maxLength[ yy$species=='bnt' & yy$year==2012 & yy$season==1 & yy$river == 'west brook' ] <- 122
  yy$maxLength[ yy$species=='bnt' & yy$year==2013 & yy$season==1 & yy$river == 'west brook' ] <- 124
  yy$maxLength[ yy$species=='bnt' & yy$year==2014 & yy$season==1 & yy$river == 'west brook' ] <- 120
  yy$maxLength[ yy$species=='bnt' & yy$year==2015 & yy$season==1 & yy$river == 'west brook' ] <- 114

  yy$maxLength[ yy$species=='bnt' & yy$year==2002 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2003 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2004 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2005 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2006 & yy$season==1 & yy$river == 'wb jimmy' ] <- 82
  yy$maxLength[ yy$species=='bnt' & yy$year==2007 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2008 & yy$season==1 & yy$river == 'wb jimmy' ] <- 87
  yy$maxLength[ yy$species=='bnt' & yy$year==2009 & yy$season==1 & yy$river == 'wb jimmy' ] <- 97
  yy$maxLength[ yy$species=='bnt' & yy$year==2010 & yy$season==1 & yy$river == 'wb jimmy' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2011 & yy$season==1 & yy$river == 'wb jimmy' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2012 & yy$season==1 & yy$river == 'wb jimmy' ] <- 100
  yy$maxLength[ yy$species=='bnt' & yy$year==2013 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$species=='bnt' & yy$year==2014 & yy$season==1 & yy$river == 'wb jimmy' ] <- 94
  yy$maxLength[ yy$species=='bnt' & yy$year==2015 & yy$season==1 & yy$river == 'wb jimmy' ] <- 93

  # so few fish, just kept what we had for bkt
  yy$maxLength[ yy$species=='bnt' & yy$year==2002 & yy$season==1 & yy$river == 'wb mitchell' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2003 & yy$season==1 & yy$river == 'wb mitchell' ] <- 110
  yy$maxLength[ yy$species=='bnt' & yy$year==2004 & yy$season==1 & yy$river == 'wb mitchell' ] <- 109
  yy$maxLength[ yy$species=='bnt' & yy$year==2005 & yy$season==1 & yy$river == 'wb mitchell' ] <- 107
  yy$maxLength[ yy$species=='bnt' & yy$year==2006 & yy$season==1 & yy$river == 'wb mitchell' ] <- 88
  yy$maxLength[ yy$species=='bnt' & yy$year==2007 & yy$season==1 & yy$river == 'wb mitchell' ] <- 83
  yy$maxLength[ yy$species=='bnt' & yy$year==2008 & yy$season==1 & yy$river == 'wb mitchell' ] <- 102
  yy$maxLength[ yy$species=='bnt' & yy$year==2009 & yy$season==1 & yy$river == 'wb mitchell' ] <- 114
  yy$maxLength[ yy$species=='bnt' & yy$year==2010 & yy$season==1 & yy$river == 'wb mitchell' ] <- 132
  yy$maxLength[ yy$species=='bnt' & yy$year==2011 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2012 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2013 & yy$season==1 & yy$river == 'wb mitchell' ] <- 90
  yy$maxLength[ yy$species=='bnt' & yy$year==2014 & yy$season==1 & yy$river == 'wb mitchell' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2015 & yy$season==1 & yy$river == 'wb mitchell' ] <- 78

  # no fish, just kept what we had for bkt
  yy$maxLength[ yy$species=='bnt' & yy$year==2002 & yy$season==1 & yy$river == 'wb obear' ] <- 95
  yy$maxLength[ yy$species=='bnt' & yy$year==2003 & yy$season==1 & yy$river == 'wb obear' ] <- 90
  yy$maxLength[ yy$species=='bnt' & yy$year==2004 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$species=='bnt' & yy$year==2005 & yy$season==1 & yy$river == 'wb obear' ] <- 83
  yy$maxLength[ yy$species=='bnt' & yy$year==2006 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2007 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$species=='bnt' & yy$year==2008 & yy$season==1 & yy$river == 'wb obear' ] <- 94
  yy$maxLength[ yy$species=='bnt' & yy$year==2009 & yy$season==1 & yy$river == 'wb obear' ] <- 103
  yy$maxLength[ yy$species=='bnt' & yy$year==2010 & yy$season==1 & yy$river == 'wb obear' ] <- 106
  yy$maxLength[ yy$species=='bnt' & yy$year==2011 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2012 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2013 & yy$season==1 & yy$river == 'wb obear' ] <- 100
  yy$maxLength[ yy$species=='bnt' & yy$year==2014 & yy$season==1 & yy$river == 'wb obear' ] <- 92
  yy$maxLength[ yy$species=='bnt' & yy$year==2015 & yy$season==1 & yy$river == 'wb obear' ] <- 84

  yy <- yy[ order(yy$species,yy$year,yy$riverOrdered,yy$season),]
  cutoffYOYInclSpring1DATA <- array( yy$maxLength, c(nSeasons,nRivers,nYears,nSpecies) )

  save(yy,cutoffYOYInclSpring1DATA,file = paste0('./data/cutoffYOYInclSpring1DATA_',dr,'.RData'))

  #check cutOffs
  #  yy[  yy$season==4 & yy$river == 'wb obear' ,c('year','maxLength')]

  # river <- 'wb mitchell'#,
  #  river <- 'west brook'#'wb obear'##'#''#'wb obear' #
  river <- 'wb obear'##'#''#'wb obear' #
  # river <- 'wb jimmy'
  #
  spp <- "bkt"

  #   ggplot( cd[cd$river==river & !is.na(cd$season) & cd$species == spp,], aes(observedLength) ) +
  #     geom_histogram( binwidth=3 )+
  #     geom_vline(aes(xintercept=maxLength), yy[yy$year >= 2002 & yy$river == river & yy$species == spp,])+
  #     facet_grid(season~year ) +
  #     ggtitle(paste(river,spp))

}

#'Get biomass deltas
#'
#'@param d dataframe created with getCoreData()
#'@return a data frame including meanBiomassAllSppStdDelta (changes in biomass with bioomass by sample standardized by all species) and meanBiomassStdDelta (changes in biomass with bioomass by sample standardized by each species)
#'@export

# addBiomassDeltas <- function(d){
#
#   # first get biomass means and deltas by cohort and ageInSamples
#   meanBiomassCohort <- d %>%
#     group_by(species, cohort, ageInSamples,year,season) %>%
#     summarize( meanBiomassAllSppStd = mean(biomassAllSppStd, na.rm=T ),
#                sdBiomassAllSppStd = sd(biomassAllSppStd, na.rm=T ),
#                meanBiomassStd = mean(biomassStd, na.rm=T ),
#                sdBiomassStd = sd(biomassStd, na.rm=T ), n=n()) %>%
#     group_by(species,cohort) %>%
#     mutate( meanBiomassAllSppStdLag = lead(meanBiomassAllSppStd),
#             meanBiomassAllSppStdDelta = meanBiomassAllSppStd - meanBiomassAllSppStdLag,
#             meanBiomassStdLag = lead(meanBiomassStd),
#             meanBiomassStdDelta = meanBiomassStd - meanBiomassStdLag )
#
#   #ggplot( filter(meanBiomassCohort, species == "bkt", n>25), aes(ageInSamples,meanBiomassStdDelta, color = factor(cohort))) + geom_point() + geom_line()# + facet_wrap(~year)
#   #ggplot( filter(meanBiomassCohort, species == "bkt", n>20), aes(season,meanBiomassAllSppStdDelta, color = factor(cohort))) + geom_point() + geom_line() + facet_wrap(~year)
#
#   # get means across cohorts
#   biomassDeltaMeans <- meanBiomassCohort %>%
#     group_by(species, year,season) %>%
#     summarize( meanBiomassAllSppStdDelta = mean(meanBiomassAllSppStdDelta, na.rm=T),
#                meanBiomassStdDelta = mean(meanBiomassStdDelta, na.rm=T))
#
#   #ggplot( filter(biomassDeltaMeans), aes(season,meanBiomassAllSppStdDelta, color=species)) + geom_point() + geom_line() + facet_wrap(~year)
#   #ggplot( filter(biomassDeltaMeans), aes(season,meanBiomassStdDelta, color=species)) + geom_point() + geom_line() + facet_wrap(~year)
#
#   d <- left_join(d,biomassDeltaMeans)
#   return( d )
# }
#
