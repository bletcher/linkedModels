library(ggplot2)
load('D:/PITTAGMAIN/CMR Analyses/Hierach_Bugs/allSpp/dMDataOutBKT2002_2012.RData')

####################################
# create a data frame of max lengths for YOYs from Matt,..

yoySizeLimits <- read.csv(file='D:/PITTAGMAIN/CMR Analyses/Hierach_Bugs/allSpp/yoySizeLimits.csv', header=T)
yoySizeLimits$year <- yoySizeLimits$YOS
yoySizeLimits$river <- yoySizeLimits$River

snOrigSeason <- data.frame(unique(cbind(dMData$sampleNumOrig,dMData$season)))
names(snOrigSeason) <- c('Sample', 'season')

y2 <- merge(x=yoySizeLimits, y=snOrigSeason, by='Sample', all.x=T)
y2$year <- y2$YOS


yTemplate <- data.frame( year=rep(2002:2012, each=5*4), river= rep(c('DEAD', 'WEST BROOK','WB JIMMY','WB MITCHELL','WB OBEAR'),each=4), season=1:4 )
y <- merge( x=yTemplate, y=y2, by=c('river','year','season'), all.x=T)

y$maxLength <- ifelse( is.na(y$Max.Length), 90, y$Max.Length) # fill in missing obs with 90
y$maxLength <- ifelse( y$season == 1, 50, y$maxLength )

# some visual fixes
y$maxLength[ y$year==2003 & y$season==3 & y$river == 'WB JIMMY' ] <- 85
y$maxLength[ y$year==2011 & y$season==2 & y$river == 'WB OBEAR' ] <- 70
y$maxLength[ y$year==2010 & y$season==3 & y$river == 'WB OBEAR' ] <- 72
y$maxLength[ y$year==2011 & y$season==3 & y$river == 'WEST BROOK' ] <- 100
y$maxLength[ y$year==2011 & y$season==4 & y$river == 'WEST BROOK' ] <- 110


y$riverOrdered <- factor(y$river, levels=c('DEAD','WEST BROOK','WB JIMMY','WB MITCHELL','WB OBEAR'), ordered=T)
y <- y[ order(y$year,y$riverOrdered,y$season),]
cutoffYOYDATA <- array( y$maxLength, c(4,5,11) )

setwd('D:/PITTAGMAIN/CMR Analyses/Hierach_Bugs/allSpp')
save(cutoffYOYDATA,file='cutoffYOYDATA.RData')
#######################################

river <- 'WB MITCHELL'#,
river <- 'WEST BROOK'#'WB OBEAR'##'#''#'WB OBEAR' #
river <- 'WB OBEAR'##'#''#'WB OBEAR' #
river <- 'WB JIMMY'

ggplot( dMData[dMData$river==river & !is.na(dMData$season),], aes(length) ) +
  geom_histogram( binwidth=3 )+
  geom_vline(aes(xintercept=maxLength), y[y$year>=2002 & y$river==river,])+
  facet_grid(season~year ) +
  ggtitle(river)



# incorporating the spring sample 1+ fish into the yoy category
yy <- y

yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 110
yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 110
yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 105
yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 115
yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 110
yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 112
yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 110
yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 125
yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 125
yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 110
yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'WEST BROOK' ] <- 112

yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 95
yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 95
yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 95
yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 92
yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 82
yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 92
yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 87
yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 93
yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 92
yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 92
yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'WB JIMMY' ] <- 90

yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 95
yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 110
yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 109
yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 107
yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 88
yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 83
yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 102
yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 114
yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 132
yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 92
yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'WB MITCHELL' ] <- 92

yy$maxLength[ yy$year==2002 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 95
yy$maxLength[ yy$year==2003 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 95
yy$maxLength[ yy$year==2004 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 94
yy$maxLength[ yy$year==2005 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 83
yy$maxLength[ yy$year==2006 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 92
yy$maxLength[ yy$year==2007 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 103
yy$maxLength[ yy$year==2008 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 94
yy$maxLength[ yy$year==2009 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 103
yy$maxLength[ yy$year==2010 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 106
yy$maxLength[ yy$year==2011 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 92
yy$maxLength[ yy$year==2012 & yy$season==1 & yy$river == 'WB OBEAR' ] <- 92

yy$riverOrdered <- factor(yy$river, levels=c('DEAD','WEST BROOK','WB JIMMY','WB MITCHELL','WB OBEAR'), ordered=T)
yy <- yy[ order(yy$year,yy$riverOrdered,yy$season),]
cutoffYOYInclSpring1DATA <- array( yy$maxLength, c(4,5,11) )

setwd('D:/PITTAGMAIN/CMR Analyses/Hierach_Bugs/allSpp')
save(yy,cutoffYOYInclSpring1DATA,file='cutoffYOYInclSpring1DATA.RData')

#check cutOffs
yy[  yy$season==4 & yy$river == 'WB OBEAR' ,c('year','maxLength')]

river <- 'WB MITCHELL'#,
river <- 'WEST BROOK'#'WB OBEAR'##'#''#'WB OBEAR' #
river <- 'WB OBEAR'##'#''#'WB OBEAR' #
river <- 'WB JIMMY'

ggplot( dMData[dMData$river==river & !is.na(dMData$season),], aes(length) ) +
  geom_histogram( binwidth=3 )+
  geom_vline(aes(xintercept=maxLength), yy[yy$year>=2002 & yy$river==river,])+
  facet_grid(season~year ) +
  ggtitle(river)
